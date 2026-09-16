rm(list = ls())
gc()

options(stringsAsFactors = FALSE, contrasts = c("contr.treatment", "contr.poly"), scipen = 999, width = 180)
ANALYSIS_VERSION <- "2026-09-16-publication-v2.0"

paper2_analysis_file <- function() {
  for (frame_index in rev(seq_len(sys.nframe()))) {
    source_file <- tryCatch(sys.frame(frame_index)$ofile, error = function(e) NULL)
    if (!is.null(source_file) && length(source_file) == 1L && nzchar(source_file)) return(normalizePath(source_file, winslash = "/", mustWork = FALSE))
  }
  args <- commandArgs(FALSE)
  hit <- grep("^--file=", args)
  if (length(hit)) return(normalizePath(sub("^--file=", "", args[hit[1]]), winslash = "/", mustWork = FALSE))
  normalizePath(file.path(getwd(), "Paper2_Final_Analysis_and_Figures.R"), winslash = "/", mustWork = FALSE)
}

ANALYSIS_FILE <- paper2_analysis_file()
PROJECT_DIR <- normalizePath(Sys.getenv("PAPER2_PROJECT_DIR", unset = dirname(ANALYSIS_FILE)), winslash = "/", mustWork = FALSE)
INPUT_DIR <- file.path(PROJECT_DIR, "data_consolidated")
OUTPUT_ROOT <- file.path(PROJECT_DIR, "final_results_output")
PERIOD_FILE <- file.path(INPUT_DIR, "PERIOD_GRAPH_LEVEL.csv")
BASIN_FILE <- file.path(INPUT_DIR, "BASIN_SYSTEM_LEVEL.csv")
PERTURBATION_FILE <- file.path(INPUT_DIR, "PERTURBATION_TRIAL_LEVEL.csv")
HAMMING_FILE <- file.path(INPUT_DIR, "HAMMING_SYSTEM_TIMESTEP.csv")

analysis_packages <- c("data.table", "mgcv", "ggplot2", "patchwork", "scales", "knitr", "svglite", "survival")
missing_packages <- analysis_packages[!vapply(analysis_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) install.packages(missing_packages, repos = "https://cloud.r-project.org", dependencies = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(mgcv)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(knitr)
  library(svglite)
  library(survival)
})

CI_LEVEL <- 0.95
RANDOM_SEED <- 20260814L
CV_FOLDS <- 5L
RUN_CROSS_VALIDATION <- tolower(Sys.getenv("PAPER2_RUN_CV", unset = "true")) %in% c("true", "1", "yes")
MODEL_THREADS <- suppressWarnings(as.integer(Sys.getenv("PAPER2_MODEL_THREADS", unset = max(1L, parallel::detectCores(logical = TRUE) - 1L))))
if (!is.finite(MODEL_THREADS) || MODEL_THREADS < 1L) MODEL_THREADS <- 1L

RUN_ID <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(OUTPUT_ROOT, paste0("FINAL_RESULTS_", RUN_ID))
TABLE_DIR <- file.path(OUTPUT_DIR, "tables_csv")
LATEX_DIR <- file.path(OUTPUT_DIR, "tables_latex")
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
MODEL_DIR <- file.path(OUTPUT_DIR, "models")
AUDIT_DIR <- file.path(OUTPUT_DIR, "audit")
for (path in c(OUTPUT_ROOT, OUTPUT_DIR, TABLE_DIR, LATEX_DIR, FIGURE_DIR, MODEL_DIR, AUDIT_DIR)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
writeLines(normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE), file.path(OUTPUT_ROOT, "LATEST_FINAL_RESULTS_OUTPUT.txt"))
required_input_files <- c(PERIOD_FILE, BASIN_FILE, PERTURBATION_FILE, HAMMING_FILE)
if (any(!file.exists(required_input_files))) stop("Missing required input files: ", paste(required_input_files[!file.exists(required_input_files)], collapse = ", "), call. = FALSE)

canonical_name <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))
safe_numeric <- function(x) suppressWarnings(as.numeric(as.character(x)))
safe_logical <- function(x) {
  if (is.logical(x)) return(x)
  value <- tolower(trimws(as.character(x)))
  output <- rep(NA, length(value))
  output[value %in% c("1", "true", "t", "yes", "y", "returned", "found")] <- TRUE
  output[value %in% c("0", "false", "f", "no", "n", "notreturned", "notfound")] <- FALSE
  numeric_value <- suppressWarnings(as.numeric(value))
  output[is.na(output) & is.finite(numeric_value)] <- numeric_value[is.na(output) & is.finite(numeric_value)] != 0
  output
}

find_column <- function(DT, aliases, required = FALSE, label = "column") {
  keys <- canonical_name(names(DT))
  hits <- match(canonical_name(aliases), keys, nomatch = 0L)
  hits <- hits[hits > 0L]
  if (length(hits)) return(names(DT)[hits[1L]])
  if (required) stop("Missing ", label, ". Expected one of: ", paste(aliases, collapse = ", "), ". Available columns: ", paste(names(DT), collapse = ", "), call. = FALSE)
  NA_character_
}

rename_standard <- function(DT, target, aliases, required = FALSE) {
  source <- find_column(DT, unique(c(target, aliases)), required, target)
  if (!is.na(source) && source != target) setnames(DT, source, target)
  invisible(DT)
}

read_input <- function(path) {
  if (!file.exists(path)) stop("Required input file not found: ", path, call. = FALSE)
  extension <- tolower(tools::file_ext(path))
  if (extension == "rds") return(as.data.table(readRDS(path)))
  if (extension %in% c("csv", "gz", "txt", "tsv")) return(fread(path, showProgress = FALSE))
  stop("Unsupported input file: ", path, call. = FALSE)
}

write_csv <- function(x, filename) {
  path <- file.path(TABLE_DIR, filename)
  fwrite(as.data.table(x), path, na = "")
  invisible(path)
}

write_latex <- function(x, filename, digits = 4) {
  path <- file.path(LATEX_DIR, filename)
  text <- knitr::kable(as.data.frame(x), format = "latex", booktabs = TRUE, digits = digits, escape = TRUE, row.names = FALSE)
  writeLines(text, path, useBytes = TRUE)
  invisible(path)
}

rmse <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  sqrt(mean((observed[keep] - predicted[keep])^2))
}

mae <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  mean(abs(observed[keep] - predicted[keep]))
}

r_squared <- function(observed, predicted) {
  keep <- is.finite(observed) & is.finite(predicted)
  observed <- observed[keep]
  predicted <- predicted[keep]
  denominator <- sum((observed - mean(observed))^2)
  if (!is.finite(denominator) || denominator <= 0) return(NA_real_)
  1 - sum((observed - predicted)^2) / denominator
}

mean_summary <- function(DT, value_column, group_columns) {
  DT[
    is.finite(get(value_column)),
    {
      values <- get(value_column)
      n <- length(values)
      estimate <- mean(values)
      standard_deviation <- if (n > 1L) sd(values) else NA_real_
      standard_error <- if (n > 1L) standard_deviation / sqrt(n) else NA_real_
      critical <- if (n > 1L) qt(1 - (1 - CI_LEVEL) / 2, n - 1L) else NA_real_
      list(N = n, Estimate = estimate, SD = standard_deviation, SE = standard_error, Lower = estimate - critical * standard_error, Upper = estimate + critical * standard_error)
    },
    by = group_columns
  ]
}

endpoint_difference <- function(DT, value_column, treatment_column, reference_value, target_value, group_columns) {
  reference <- mean_summary(DT[abs(get(treatment_column) - reference_value) < 1e-10], value_column, group_columns)
  target <- mean_summary(DT[abs(get(treatment_column) - target_value) < 1e-10], value_column, group_columns)
  setnames(reference, c("N", "Estimate", "SD", "SE", "Lower", "Upper"), paste0(c("N", "Estimate", "SD", "SE", "Lower", "Upper"), "Reference"))
  setnames(target, c("N", "Estimate", "SD", "SE", "Lower", "Upper"), paste0(c("N", "Estimate", "SD", "SE", "Lower", "Upper"), "Target"))
  output <- merge(reference, target, by = group_columns, all = FALSE)
  output[, `:=`(
    ReferenceValue = reference_value,
    TargetValue = target_value,
    Difference = EstimateTarget - EstimateReference,
    DifferenceSE = sqrt((SDReference^2 / NReference) + (SDTarget^2 / NTarget)),
    DifferenceDF = ((SDReference^2 / NReference) + (SDTarget^2 / NTarget))^2 / (((SDReference^2 / NReference)^2 / pmax(NReference - 1, 1)) + ((SDTarget^2 / NTarget)^2 / pmax(NTarget - 1, 1)))
  )]
  output[, `:=`(
    DifferenceLower = Difference - qt(1 - (1 - CI_LEVEL) / 2, DifferenceDF) * DifferenceSE,
    DifferenceUpper = Difference + qt(1 - (1 - CI_LEVEL) / 2, DifferenceDF) * DifferenceSE,
    PValue = 2 * pt(-abs(Difference / DifferenceSE), DifferenceDF)
  )]
  output[]
}

paper2_theme <- function(base_size = 10) {
  theme_classic(base_size = base_size, base_family = "sans") +
    theme(
      axis.title = element_text(size = base_size, colour = "black"),
      axis.text = element_text(size = base_size - 0.5, colour = "black"),
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks.length = grid::unit(1.5, "mm"),
      legend.position = "right",
      legend.direction = "vertical",
      legend.box = "vertical",
      legend.justification = "center",
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 0.5),
      legend.key.width = grid::unit(9, "mm"),
      legend.key.height = grid::unit(5.5, "mm"),
      legend.spacing.y = grid::unit(1.2, "mm"),
      legend.margin = margin(0, 0, 0, 5),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold", margin = margin(2, 2, 3, 2)),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_text(size = base_size + 2, face = "bold"),
      plot.tag.position = c(0.01, 0.99),
      panel.spacing = grid::unit(5, "mm"),
      plot.margin = margin(6, 8, 6, 7)
    )
}

prep_palette <- c("0.10" = "#0072B2", "0.24" = "#E69F00", "0.30" = "#009E73", "0.41" = "#CC79A7", "0.50" = "#D55E00")
qc_palette <- c("0.00" = "#000000", "0.25" = "#0072B2", "0.50" = "#009E73", "0.75" = "#E69F00", "1.00" = "#D55E00")
qc_shapes <- c("0.00" = 16, "0.25" = 17, "0.50" = 15, "0.75" = 18, "1.00" = 8)
prep_linetypes <- c("0.10" = "solid", "0.24" = "dashed", "0.30" = "dotdash", "0.41" = "longdash", "0.50" = "twodash")
series_linetypes <- c("solid", "dashed", "dotdash", "longdash", "twodash", "dotted")
series_shapes <- c(16, 17, 15, 18, 8, 3)
format_qc_axis <- function(x) {
  output <- sprintf("%.2f", x)
  output[abs(x) < 1e-10] <- "0"
  output[abs(x - 1) < 1e-10] <- "1"
  output
}

save_figure <- function(plot_object, filename, height, width = 7.2) {
  png_path <- file.path(FIGURE_DIR, paste0(filename, ".png"))
  pdf_path <- file.path(FIGURE_DIR, paste0(filename, ".pdf"))
  svg_path <- file.path(FIGURE_DIR, paste0(filename, ".svg"))
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggsave(png_path, plot_object, width = width, height = height, units = "in", dpi = 600, bg = "white")
  ggsave(pdf_path, plot_object, width = width, height = height, units = "in", device = pdf_device, bg = "white")
  ggsave(svg_path, plot_object, width = width, height = height, units = "in", device = svglite::svglite, bg = "white")
  output_paths <- c(png = png_path, pdf = pdf_path, svg = svg_path)
  if (any(!file.exists(output_paths)) || any(file.info(output_paths)$size <= 0)) stop("Figure export failed for: ", filename, call. = FALSE)
  invisible(output_paths)
}

period_graph <- read_input(PERIOD_FILE)
rename_standard(period_graph, "GraphID", c("SystemID", "graph_id", "GraphRuleID", "NetworkID", "ReplicateID"), TRUE)
rename_standard(period_graph, "N", c("NumberOfNodes", "Nodes", "n_nodes"), TRUE)
rename_standard(period_graph, "d_input", c("d", "DegreeInput", "AvgDegreeInput", "degree_control"), TRUE)
rename_standard(period_graph, "kappa", c("Kappa", "WSKappa", "k", "RingLatticeParameter"), FALSE)
rename_standard(period_graph, "p", c("p_rewire", "RewireProb", "RewireProbability", "beta"), TRUE)
rename_standard(period_graph, "qc", c("q_c", "CanalizingProbability", "CanalisingProbability"), TRUE)
rename_standard(period_graph, "prep", c("p_rep", "InhibitionFraction", "RepressionProbability"), TRUE)
rename_standard(period_graph, "Y", c("MeanLogAttractorPeriod", "AvgLogPeriod", "MeanLogPeriod", "mean_log_period"), TRUE)
rename_standard(period_graph, "NFound", c("FoundRuns", "n_found", "FoundCount", "ValidTrajectories"), TRUE)
rename_standard(period_graph, "NAttempted", c("Runs", "n_attempted", "AttemptedTrajectories", "NumberInitialStates"), TRUE)
rename_standard(period_graph, "FoundRate", c("found_rate", "DetectionRate"), TRUE)
rename_standard(period_graph, "FixedPointRate", c("fixed_point_rate", "FixedRate", "PeriodOneRate"), TRUE)

if (!"kappa" %in% names(period_graph)) period_graph[, kappa := pmax(floor(safe_numeric(d_input) / 2), 1)]
period_graph[, `:=`(
  N = as.integer(round(safe_numeric(N))),
  d_input = as.integer(round(safe_numeric(d_input))),
  kappa = as.integer(round(safe_numeric(kappa))),
  p = round(safe_numeric(p), 4),
  qc = round(safe_numeric(qc), 2),
  prep = round(safe_numeric(prep), 2),
  Y = safe_numeric(Y),
  NFound = safe_numeric(NFound),
  NAttempted = safe_numeric(NAttempted),
  FoundRate = safe_numeric(FoundRate),
  FixedPointRate = safe_numeric(FixedPointRate)
)]
period_graph <- period_graph[is.finite(N) & is.finite(d_input) & is.finite(kappa) & is.finite(p) & is.finite(qc) & is.finite(prep) & is.finite(Y) & is.finite(FixedPointRate)]
if (!nrow(period_graph)) stop("No valid graph-level period rows remain.", call. = FALSE)
if (anyDuplicated(period_graph$GraphID)) stop("GraphID is not unique in the period file.", call. = FALSE)

degree_audit <- period_graph[, .(GraphRows = .N, KappaValues = paste(sort(unique(kappa)), collapse = ",")), by = d_input]
degree_audit[, ExpectedKappa := pmax(floor(d_input / 2), 1)]
degree_audit[, EffectiveUndirectedLatticeDegree := 2L * ExpectedKappa]
degree_audit[, MappingValid := KappaValues == as.character(ExpectedKappa)]
write_csv(degree_audit, "01_DEGREE_MAPPING_AUDIT.csv")
if (any(!degree_audit$MappingValid)) stop("The d-to-kappa mapping is inconsistent. Inspect 01_DEGREE_MAPPING_AUDIT.csv.", call. = FALSE)

design_cells <- period_graph[, .(
  Y = mean(Y),
  FixedPointRate = mean(FixedPointRate),
  Graphs = .N,
  DInputValues = paste(sort(unique(d_input)), collapse = ","),
  MeanFoundRate = mean(FoundRate, na.rm = TRUE),
  TotalFound = sum(NFound, na.rm = TRUE),
  TotalAttempted = sum(NAttempted, na.rm = TRUE)
), by = .(N, kappa, p, qc, prep)]
setorder(design_cells, N, kappa, p, qc, prep)
if (min(design_cells$Graphs) < 2L) stop("At least one design cell has fewer than two graph replicates.", call. = FALSE)

data_audit <- data.table(
  Quantity = c("Graph-level rows", "Kappa design cells", "Unique graph IDs", "N levels", "Original d input levels", "Actual kappa levels", "p levels", "q_c levels", "p_rep levels", "Minimum graphs per kappa cell", "Maximum graphs per kappa cell", "Total found trajectories", "Total attempted trajectories", "Overall detection rate"),
  Value = c(nrow(period_graph), nrow(design_cells), uniqueN(period_graph$GraphID), paste(sort(unique(period_graph$N)), collapse = ", "), paste(sort(unique(period_graph$d_input)), collapse = ", "), paste(sort(unique(period_graph$kappa)), collapse = ", "), paste(sort(unique(period_graph$p)), collapse = ", "), paste(sort(unique(period_graph$qc)), collapse = ", "), paste(sort(unique(period_graph$prep)), collapse = ", "), min(design_cells$Graphs), max(design_cells$Graphs), sum(period_graph$NFound), sum(period_graph$NAttempted), sum(period_graph$NFound) / sum(period_graph$NAttempted))
)
write_csv(data_audit, "00_DATA_AND_DESIGN_AUDIT.csv")
write_csv(design_cells, "02_MODEL_DESIGN_CELLS.csv")

observed_period <- mean_summary(period_graph, "Y", c("prep", "qc"))
setnames(observed_period, c("N", "Estimate", "SD", "SE", "Lower", "Upper"), c("Graphs", "MeanLogPeriod", "SDLogPeriod", "SELogPeriod", "LowerLogPeriod", "UpperLogPeriod"))
observed_period[, `:=`(GeometricMeanPeriod = exp(MeanLogPeriod), LowerGeometricMeanPeriod = exp(LowerLogPeriod), UpperGeometricMeanPeriod = exp(UpperLogPeriod))]
observed_fixed <- mean_summary(period_graph, "FixedPointRate", c("prep", "qc"))
setnames(observed_fixed, c("N", "Estimate", "SD", "SE", "Lower", "Upper"), c("GraphsFixed", "FixedPointRate", "SDFixedPointRate", "SEFixedPointRate", "LowerFixedPointRate", "UpperFixedPointRate"))
observed <- merge(observed_period, observed_fixed, by = c("prep", "qc"), all = TRUE)
setorder(observed, prep, qc)
write_csv(observed, "02_OBSERVED_PERIOD_FIXED_POINT.csv")

observed_period_endpoint <- endpoint_difference(period_graph, "Y", "qc", min(period_graph$qc), max(period_graph$qc), c("prep"))
observed_period_endpoint[, `:=`(PeriodRatio = exp(Difference), LowerRatio = exp(DifferenceLower), UpperRatio = exp(DifferenceUpper), PercentChange = 100 * (exp(Difference) - 1))]
observed_fixed_endpoint <- endpoint_difference(period_graph, "FixedPointRate", "qc", min(period_graph$qc), max(period_graph$qc), c("prep"))
write_csv(observed_period_endpoint, "03_RULE_ENDPOINT_PERIOD_CONTRASTS.csv")
write_csv(observed_fixed_endpoint, "03_RULE_ENDPOINT_FIXED_POINT_CONTRASTS.csv")

p_levels <- sort(unique(period_graph$p))
qc_levels <- sort(unique(period_graph$qc))
prep_levels <- sort(unique(period_graph$prep))
N_levels <- sort(unique(period_graph$N))
kappa_levels <- sort(unique(period_graph$kappa))
centres <- list(N = mean(period_graph$N), kappa = mean(period_graph$kappa), p = mean(period_graph$p), qc = mean(period_graph$qc), prep = mean(period_graph$prep))
p_text <- sprintf("%.4f", p_levels)
qc_text <- sprintf("%.2f", qc_levels)
prep_text <- sprintf("%.2f", prep_levels)

build_model_data <- function(DT) {
  output <- copy(as.data.table(DT))
  output[, `:=`(
    Nc = safe_numeric(N) - centres$N,
    kappac = safe_numeric(kappa) - centres$kappa,
    p_num = safe_numeric(p) - centres$p,
    qc_num = safe_numeric(qc) - centres$qc,
    prep_num = safe_numeric(prep) - centres$prep,
    p_fac = factor(sprintf("%.4f", safe_numeric(p)), levels = p_text),
    qc_fac = factor(sprintf("%.2f", safe_numeric(qc)), levels = qc_text),
    prep_fac = factor(sprintf("%.2f", safe_numeric(prep)), levels = prep_text)
  )]
  output
}

safe_k <- function(x, requested) {
  available <- uniqueN(round(safe_numeric(x), 10))
  max(3L, min(as.integer(requested), as.integer(available - 1L)))
}

model_graph_data <- build_model_data(period_graph)
model_selection_data <- build_model_data(design_cells)
k_N <- safe_k(model_graph_data$Nc, 9L)
k_kappa <- safe_k(model_graph_data$kappac, 6L)
k_p <- safe_k(model_graph_data$p_num, 5L)
k_qc <- safe_k(model_graph_data$qc_num, 5L)
k_prep <- safe_k(model_graph_data$prep_num, 5L)
surface_Nkappa <- sprintf("te(Nc,kappac,k=c(%d,%d),bs=c('cr','cr'))", k_N, k_kappa)
ti_N_p <- sprintf("ti(Nc,p_num,k=c(%d,%d),bs=c('cr','cr'))", min(k_N, 6L), min(k_p, 4L))
ti_kappa_p <- sprintf("ti(kappac,p_num,k=c(%d,%d),bs=c('cr','cr'))", min(k_kappa, 5L), min(k_p, 4L))
ti_N_qc <- sprintf("ti(Nc,qc_num,k=c(%d,%d),bs=c('cr','cr'))", min(k_N, 6L), min(k_qc, 4L))
ti_kappa_qc <- sprintf("ti(kappac,qc_num,k=c(%d,%d),bs=c('cr','cr'))", min(k_kappa, 5L), min(k_qc, 4L))
ti_N_prep <- sprintf("ti(Nc,prep_num,k=c(%d,%d),bs=c('cr','cr'))", min(k_N, 6L), min(k_prep, 4L))
ti_kappa_prep <- sprintf("ti(kappac,prep_num,k=c(%d,%d),bs=c('cr','cr'))", min(k_kappa, 5L), min(k_prep, 4L))
ti_kappa_qc_prep <- sprintf("ti(kappac,qc_num,prep_num,k=c(%d,%d,%d),bs=c('cr','cr','cr'))", min(k_kappa, 5L), min(k_qc, 4L), min(k_prep, 4L))

make_formula <- function(terms) as.formula(paste("Y ~", paste(unique(terms), collapse = " + ")))
model_formulas <- list(
  K1_rule_main = make_formula(c(surface_Nkappa, "p_fac", "qc_fac * prep_fac")),
  K2_rewiring_rule = make_formula(c(surface_Nkappa, "p_fac * qc_fac * prep_fac")),
  K3_all_pairwise = make_formula(c(surface_Nkappa, "p_fac * qc_fac * prep_fac", ti_N_p, ti_kappa_p, ti_N_qc, ti_kappa_qc, ti_N_prep, ti_kappa_prep)),
  K4_exploratory_kappa_qc_prep = make_formula(c(surface_Nkappa, "p_fac * qc_fac * prep_fac", ti_N_p, ti_kappa_p, ti_N_qc, ti_kappa_qc, ti_N_prep, ti_kappa_prep, ti_kappa_qc_prep))
)

formula_table <- rbindlist(lapply(names(model_formulas), function(name) data.table(Model = name, Formula = paste(deparse(model_formulas[[name]]), collapse = " "))))
write_csv(formula_table, "04_PRIMARY_MODEL_FORMULAS.csv")

fit_selection_model <- function(formula_object) {
  gam(formula_object, data = model_selection_data, weights = Graphs, family = gaussian(), method = "ML", select = TRUE, na.action = na.fail)
}

selection_models <- lapply(model_formulas, fit_selection_model)
model_comparison <- rbindlist(lapply(names(selection_models), function(name) {
  model <- selection_models[[name]]
  prediction <- as.numeric(predict(model, model_selection_data, type = "response"))
  model_summary <- summary(model)
  data.table(Model = name, AIC = AIC(model), BIC = BIC(model), EDF = sum(model$edf), AdjustedR2 = unname(model_summary$r.sq), DevianceExplained = unname(model_summary$dev.expl), RMSE = rmse(model_selection_data$Y, prediction), MAE = mae(model_selection_data$Y, prediction), R2 = r_squared(model_selection_data$Y, prediction))
}))

three_way_test <- anova(selection_models$K3_all_pairwise, selection_models$K4_exploratory_kappa_qc_prep, test = "F")
p_column <- grep("Pr", names(three_way_test), value = TRUE)
three_way_p <- if (length(p_column)) as.numeric(three_way_test[[p_column[1L]]][nrow(three_way_test)]) else NA_real_
three_way_comparison <- data.table(
  BaseModel = "K3_all_pairwise",
  ExtendedModel = "K4_exploratory_kappa_qc_prep",
  BaseAIC = AIC(selection_models$K3_all_pairwise),
  ExtendedAIC = AIC(selection_models$K4_exploratory_kappa_qc_prep),
  DeltaAIC = AIC(selection_models$K4_exploratory_kappa_qc_prep) - AIC(selection_models$K3_all_pairwise),
  BaseBIC = BIC(selection_models$K3_all_pairwise),
  ExtendedBIC = BIC(selection_models$K4_exploratory_kappa_qc_prep),
  DeltaBIC = BIC(selection_models$K4_exploratory_kappa_qc_prep) - BIC(selection_models$K3_all_pairwise),
  LikelihoodRatioP = three_way_p
)

selection_cv_prediction <- function(formula_object, fold_ids) {
  prediction <- rep(NA_real_, nrow(model_selection_data))
  for (fold in seq_len(CV_FOLDS)) {
    training <- model_selection_data[fold_ids != fold]
    testing <- model_selection_data[fold_ids == fold]
    model <- gam(formula_object, data = training, weights = Graphs, family = gaussian(), method = "ML", select = TRUE, na.action = na.fail)
    prediction[fold_ids == fold] <- as.numeric(predict(model, testing, type = "response"))
  }
  prediction
}

if (RUN_CROSS_VALIDATION) {
  set.seed(RANDOM_SEED)
  selection_folds <- sample(rep(seq_len(CV_FOLDS), length.out = nrow(model_selection_data)))
  base_cv_prediction <- selection_cv_prediction(model_formulas$K3_all_pairwise, selection_folds)
  extended_cv_prediction <- selection_cv_prediction(model_formulas$K4_exploratory_kappa_qc_prep, selection_folds)
  three_way_comparison[, `:=`(
    BaseCVRMSE = rmse(model_selection_data$Y, base_cv_prediction),
    ExtendedCVRMSE = rmse(model_selection_data$Y, extended_cv_prediction),
    DeltaCVRMSE = rmse(model_selection_data$Y, extended_cv_prediction) - rmse(model_selection_data$Y, base_cv_prediction)
  )]
} else {
  three_way_comparison[, `:=`(BaseCVRMSE = NA_real_, ExtendedCVRMSE = NA_real_, DeltaCVRMSE = NA_real_)]
}

use_extended_primary <- is.finite(three_way_comparison$LikelihoodRatioP) && three_way_comparison$LikelihoodRatioP < 0.05 && three_way_comparison$DeltaAIC <= -2
primary_model_name <- if (use_extended_primary) "K4_exploratory_kappa_qc_prep" else "K3_all_pairwise"
primary_formula <- model_formulas[[primary_model_name]]
three_way_comparison[, `:=`(SelectedAsPrimary = use_extended_primary, PrimaryModelName = primary_model_name)]

write_csv(model_comparison, "05_PRIMARY_MODEL_COMPARISON.csv")
write_csv(three_way_comparison, "15_KAPPA_QC_PREP_MODEL_COMPARISON.csv")
writeLines(capture.output(three_way_test), file.path(MODEL_DIR, "KAPPA_QC_PREP_LIKELIHOOD_RATIO_TEST.txt"))
saveRDS(selection_models, file.path(MODEL_DIR, "MODEL_SELECTION_ML.rds"))

fit_bam <- function(formula_object, data_object) {
  bam(formula_object, data = data_object, family = gaussian(), method = "fREML", select = TRUE, discrete = TRUE, nthreads = MODEL_THREADS, gc.level = 1, na.action = na.fail)
}

primary_model <- fit_bam(primary_formula, model_graph_data)
exploratory_model <- if (primary_model_name == "K4_exploratory_kappa_qc_prep") primary_model else fit_bam(model_formulas$K4_exploratory_kappa_qc_prep, model_graph_data)
primary_prediction <- as.numeric(predict(primary_model, model_graph_data, type = "response"))
primary_summary <- summary(primary_model)
saveRDS(primary_model, file.path(MODEL_DIR, "PRIMARY_CONTROLLED_PARAMETER_BAM_FREML.rds"))
saveRDS(exploratory_model, file.path(MODEL_DIR, "EXPLORATORY_KAPPA_QC_PREP_BAM_FREML.rds"))
writeLines(paste(deparse(formula(primary_model)), collapse = " "), file.path(MODEL_DIR, "PRIMARY_MODEL_FORMULA.txt"))
writeLines(capture.output(primary_summary), file.path(MODEL_DIR, "PRIMARY_MODEL_SUMMARY.txt"))

forbidden_variables <- c("C", "ClusteringCoefficient", "MSP", "MeanShortestPath", "AverageDegree", "average_degree", "kbar")
if (any(tolower(all.vars(formula(primary_model))) %in% tolower(forbidden_variables))) stop("The primary model contains an uncontrolled graph metric.", call. = FALSE)

cv_prediction <- rep(NA_real_, nrow(model_graph_data))
cv_fold <- rep(NA_integer_, nrow(model_graph_data))
if (RUN_CROSS_VALIDATION) {
  model_graph_data[, DesignCellID := sprintf("N%d_K%d_P%.4f_Q%.2f_R%.2f", N, kappa, p, qc, prep)]
  cells <- unique(model_graph_data$DesignCellID)
  set.seed(RANDOM_SEED)
  cell_folds <- data.table(DesignCellID = sample(cells))
  cell_folds[, Fold := rep(seq_len(CV_FOLDS), length.out = .N)]
  model_graph_data[cell_folds, Fold := i.Fold, on = "DesignCellID"]
  cv_fold <- model_graph_data$Fold
  for (fold in seq_len(CV_FOLDS)) {
    training <- model_graph_data[Fold != fold]
    testing <- model_graph_data[Fold == fold]
    fold_model <- fit_bam(primary_formula, training)
    cv_prediction[model_graph_data$Fold == fold] <- as.numeric(predict(fold_model, testing, type = "response"))
  }
  write_csv(model_graph_data[, .(GraphID, DesignCellID, Fold)], "06_GROUPED_CV_ASSIGNMENTS.csv")
}

model_performance <- data.table(
  Quantity = c("Graph rows", "Adjusted R2", "Deviance explained", "In-sample R2", "In-sample RMSE", "In-sample MAE", "Effective degrees of freedom", "Cross-validation folds", "Cross-validation R2", "Cross-validation RMSE", "Cross-validation MAE"),
  Value = c(nrow(model_graph_data), unname(primary_summary$r.sq), unname(primary_summary$dev.expl), r_squared(model_graph_data$Y, primary_prediction), rmse(model_graph_data$Y, primary_prediction), mae(model_graph_data$Y, primary_prediction), sum(primary_model$edf), if (RUN_CROSS_VALIDATION) CV_FOLDS else NA_real_, if (RUN_CROSS_VALIDATION) r_squared(model_graph_data$Y, cv_prediction) else NA_real_, if (RUN_CROSS_VALIDATION) rmse(model_graph_data$Y, cv_prediction) else NA_real_, if (RUN_CROSS_VALIDATION) mae(model_graph_data$Y, cv_prediction) else NA_real_)
)
write_csv(model_performance, "06_PRIMARY_MODEL_DIAGNOSTICS.csv")
if (!is.null(primary_summary$p.table)) write_csv(as.data.table(primary_summary$p.table, keep.rownames = "Term"), "06_PRIMARY_PARAMETRIC_TERMS.csv")
if (!is.null(primary_summary$s.table)) write_csv(as.data.table(primary_summary$s.table, keep.rownames = "Term"), "06_PRIMARY_SMOOTH_TERMS.csv")
k_check <- tryCatch(as.data.table(mgcv::k.check(primary_model), keep.rownames = "Smooth"), error = function(e) data.table(Error = conditionMessage(e)))
concurvity_check <- tryCatch(as.data.table(t(mgcv::concurvity(primary_model, full = TRUE)), keep.rownames = "Term"), error = function(e) data.table(Error = conditionMessage(e)))
model_residuals <- residuals(primary_model, type = "deviance")
model_hat <- tryCatch(hatvalues(primary_model), error = function(e) rep(NA_real_, length(model_residuals)))
residual_diagnostics <- data.table(
  Quantity = c("Residual minimum", "Residual first quartile", "Residual median", "Residual third quartile", "Residual maximum", "Residual standard deviation", "Maximum leverage", "Leverage 99th percentile"),
  Value = c(min(model_residuals), quantile(model_residuals, 0.25), median(model_residuals), quantile(model_residuals, 0.75), max(model_residuals), sd(model_residuals), if (any(is.finite(model_hat))) max(model_hat, na.rm = TRUE) else NA_real_, if (any(is.finite(model_hat))) quantile(model_hat, 0.99, na.rm = TRUE) else NA_real_)
)
write_csv(k_check, "06_PRIMARY_BASIS_DIMENSION_CHECK.csv")
write_csv(concurvity_check, "06_PRIMARY_CONCURVITY_CHECK.csv")
write_csv(residual_diagnostics, "06_PRIMARY_RESIDUAL_DIAGNOSTICS.csv")

scenario_data <- function(base, replacements = list()) {
  output <- copy(as.data.table(base))
  for (name in names(replacements)) output[, (name) := replacements[[name]]]
  unique(build_model_data(output[, .(N, kappa, p, qc, prep)]))
}

average_prediction <- function(model, newdata) {
  matrix <- predict(model, newdata = newdata, type = "lpmatrix")
  vector <- colMeans(matrix)
  coefficients <- coef(model)
  covariance <- vcov(model, unconditional = TRUE)
  estimate <- as.numeric(sum(vector * coefficients))
  standard_error <- as.numeric(sqrt(drop(vector %*% covariance %*% vector)))
  critical <- qnorm(1 - (1 - CI_LEVEL) / 2)
  data.table(EstimateY = estimate, SE = standard_error, LowerY = estimate - critical * standard_error, UpperY = estimate + critical * standard_error, GeometricMeanPeriod = exp(estimate), LowerPeriod = exp(estimate - critical * standard_error), UpperPeriod = exp(estimate + critical * standard_error))
}

average_contrast <- function(model, target, reference) {
  target_matrix <- predict(model, newdata = target, type = "lpmatrix")
  reference_matrix <- predict(model, newdata = reference, type = "lpmatrix")
  vector <- colMeans(target_matrix - reference_matrix)
  coefficients <- coef(model)
  covariance <- vcov(model, unconditional = TRUE)
  estimate <- as.numeric(sum(vector * coefficients))
  standard_error <- as.numeric(sqrt(drop(vector %*% covariance %*% vector)))
  critical <- qnorm(1 - (1 - CI_LEVEL) / 2)
  data.table(DeltaY = estimate, SE = standard_error, LowerDeltaY = estimate - critical * standard_error, UpperDeltaY = estimate + critical * standard_error, PeriodRatio = exp(estimate), LowerRatio = exp(estimate - critical * standard_error), UpperRatio = exp(estimate + critical * standard_error), PercentChange = 100 * (exp(estimate) - 1), Z = estimate / standard_error, PValue = 2 * pnorm(-abs(estimate / standard_error)))
}

context_filter <- function(DT, context_row, context_columns) {
  keep <- rep(TRUE, nrow(DT))
  for (name in context_columns) keep <- keep & DT[[name]] == context_row[[name]][1L]
  DT[keep]
}

contrasts_for_levels <- function(model, support, treatment, reference_value, target_values, context_columns) {
  contexts <- if (length(context_columns)) unique(support[, ..context_columns]) else data.table(Context = 1)
  rows <- list()
  index <- 0L
  for (i in seq_len(nrow(contexts))) {
    base <- if (length(context_columns)) context_filter(support, contexts[i], context_columns) else support
    base <- unique(base[, .(N, kappa, p, qc, prep)])
    for (target_value in target_values) {
      index <- index + 1L
      heading <- if (length(context_columns)) copy(contexts[i]) else data.table()
      heading[, `:=`(Treatment = treatment, Reference = reference_value, Target = target_value, StandardisationRows = nrow(base))]
      rows[[index]] <- cbind(heading, average_contrast(model, scenario_data(base, setNames(list(target_value), treatment)), scenario_data(base, setNames(list(reference_value), treatment))))
    }
  }
  rbindlist(rows, fill = TRUE)
}

contrast_between_levels <- function(model, support, treatment, reference_value, target_value, context_columns) {
  contrasts_for_levels(model, support, treatment, reference_value, target_value, context_columns)
}

emm_for_variables <- function(model, support, variables) {
  combinations <- unique(support[, ..variables])
  rows <- vector("list", nrow(combinations))
  for (i in seq_len(nrow(combinations))) {
    base <- context_filter(support, combinations[i], variables)
    base <- unique(base[, .(N, kappa, p, qc, prep)])
    rows[[i]] <- cbind(copy(combinations[i]), data.table(StandardisationRows = nrow(base)), average_prediction(model, build_model_data(base)))
  }
  rbindlist(rows, fill = TRUE)
}

support <- unique(design_cells[, .(N, kappa, p, qc, prep)])
support_at_common_kappa <- support[kappa == min(kappa_levels)]
support_at_maximum_N <- support[N == max(N_levels)]
qc_contrasts <- contrasts_for_levels(primary_model, support, "qc", min(qc_levels), qc_levels[qc_levels != min(qc_levels)], c("prep"))
qc_contrasts[, AdjustedPValue := p.adjust(PValue, method = "holm"), by = prep]
prep_reference <- prep_levels[which.min(abs(prep_levels - 0.50))]
prep_contrasts <- contrasts_for_levels(primary_model, support, "prep", prep_reference, prep_levels[prep_levels != prep_reference], c("qc"))
prep_contrasts[, AdjustedPValue := p.adjust(PValue, method = "holm"), by = qc]
qc_prep_emm <- emm_for_variables(primary_model, support, c("qc", "prep"))
qc_endpoint_by_prep <- contrast_between_levels(primary_model, support, "qc", min(qc_levels), max(qc_levels), c("prep"))
p_qc_prep_emm <- emm_for_variables(primary_model, support, c("p", "qc", "prep"))
p_endpoints <- contrast_between_levels(primary_model, support, "p", min(p_levels), max(p_levels), c("qc", "prep"))
qc_endpoint_by_p_prep <- contrast_between_levels(primary_model, support, "qc", min(qc_levels), max(qc_levels), c("p", "prep"))
N_kappa_emm <- emm_for_variables(primary_model, support, c("N", "kappa"))
N_kappa_qc_emm <- emm_for_variables(primary_model, support, c("N", "kappa", "qc"))
N_kappa_prep_emm <- emm_for_variables(primary_model, support, c("N", "kappa", "prep"))
kappa_prep_contrasts <- contrasts_for_levels(primary_model, support_at_maximum_N, "prep", prep_reference, prep_levels[prep_levels != prep_reference], c("kappa"))
kappa_prep_contrasts[, AdjustedPValue := p.adjust(PValue, method = "holm"), by = kappa]
N_p_emm <- emm_for_variables(primary_model, support_at_common_kappa, c("N", "p"))
N_qc_emm <- emm_for_variables(primary_model, support_at_common_kappa, c("N", "qc"))
N_prep_emm <- emm_for_variables(primary_model, support_at_common_kappa, c("N", "prep"))
kappa_p_emm <- emm_for_variables(primary_model, support_at_maximum_N, c("kappa", "p"))
kappa_qc_emm <- emm_for_variables(primary_model, support_at_maximum_N, c("kappa", "qc"))
kappa_prep_emm <- emm_for_variables(primary_model, support_at_maximum_N, c("kappa", "prep"))
kappa_qc_prep_emm <- emm_for_variables(exploratory_model, support, c("kappa", "qc", "prep"))
kappa_qc_prep_endpoints <- contrast_between_levels(exploratory_model, support, "qc", min(qc_levels), max(qc_levels), c("N", "kappa", "prep"))

write_csv(qc_contrasts, "07_ADJUSTED_QC_CONTRASTS.csv")
write_csv(prep_contrasts, "08_ADJUSTED_PREP_CONTRASTS.csv")
write_csv(qc_prep_emm, "09_QC_PREP_ADJUSTED_ESTIMATES.csv")
write_csv(qc_endpoint_by_prep, "09_QC_ENDPOINTS_BY_PREP.csv")
write_csv(p_qc_prep_emm, "09_P_QC_PREP_ESTIMATES.csv")
write_csv(qc_endpoint_by_p_prep, "09_QC_ENDPOINTS_AT_EACH_P_PREP.csv")
write_csv(p_endpoints, "10_REWIRING_CONTRASTS.csv")
write_csv(N_kappa_emm, "11_N_KAPPA_ESTIMATES.csv")
write_csv(kappa_qc_prep_endpoints, "12_KAPPA_QC_CONTRASTS.csv")
write_csv(kappa_prep_emm, "13_KAPPA_PREP_ESTIMATES.csv")
write_csv(kappa_prep_contrasts, "13_KAPPA_PREP_CONTRASTS.csv")
write_csv(kappa_qc_prep_emm, "15_EXPLORATORY_KAPPA_QC_PREP_ESTIMATES.csv")

all_pairwise <- rbindlist(list(
  copy(N_kappa_emm)[, View := "N_by_kappa"],
  copy(N_p_emm)[, View := "N_by_p"],
  copy(N_qc_emm)[, View := "N_by_qc"],
  copy(N_prep_emm)[, View := "N_by_prep"],
  copy(kappa_p_emm)[, View := "kappa_by_p"],
  copy(kappa_qc_emm)[, View := "kappa_by_qc"],
  copy(kappa_prep_emm)[, View := "kappa_by_prep"]
), fill = TRUE)
setcolorder(all_pairwise, c("View", setdiff(names(all_pairwise), "View")))
write_csv(all_pairwise, "14_ALL_SEVEN_PAIRWISE_CONTROLLED_PARAMETER_VIEWS.csv")
pairwise_standardisation <- data.table(
  View = c("N_by_kappa", "N_by_p", "N_by_qc", "N_by_prep", "kappa_by_p", "kappa_by_qc", "kappa_by_prep"),
  StructuralRestriction = c("Observed valid N-kappa support", paste0("kappa = ", min(kappa_levels)), paste0("kappa = ", min(kappa_levels)), paste0("kappa = ", min(kappa_levels)), paste0("N = ", max(N_levels)), paste0("N = ", max(N_levels)), paste0("N = ", max(N_levels))),
  Reason = c("Shows the triangular graph-construction design without extrapolation", "Compares N on the common kappa support", "Compares N on the common kappa support", "Compares N on the common kappa support", "Compares kappa where every kappa value was simulated", "Compares kappa where every kappa value was simulated", "Compares kappa where every kappa value was simulated")
)
write_csv(pairwise_standardisation, "14_PAIRWISE_STANDARDISATION_AUDIT.csv")

basin_system <- read_input(BASIN_FILE)
rename_standard(basin_system, "SystemID", c("GraphID", "NetworkID", "BooleanSystemID"), TRUE)
rename_standard(basin_system, "N", c("NumberOfNodes", "Nodes"), TRUE)
rename_standard(basin_system, "d", c("d_input", "AvgDegreeInput", "DegreeInput"), TRUE)
rename_standard(basin_system, "p", c("RewireProb", "p_rewire", "RewireProbability"), TRUE)
rename_standard(basin_system, "qc", c("q_c", "CanalizingProbability", "CanalisingProbability"), TRUE)
rename_standard(basin_system, "prep", c("p_rep", "InhibitionFraction"), TRUE)
rename_standard(basin_system, "DistinctAttractors", c("DistinctAttractorsReached", "NumberDistinctAttractors"), TRUE)
rename_standard(basin_system, "LargestBasinFraction", c("LargestBasin", "MaxBasinFraction", "LargestObservedBasin"), TRUE)
rename_standard(basin_system, "BasinEntropy", c("ObservedBasinEntropy", "EmpiricalBasinEntropy", "Entropy"), TRUE)
basin_system[, `:=`(
  N = as.integer(round(safe_numeric(N))),
  d = as.integer(round(safe_numeric(d))),
  p = round(safe_numeric(p), 4),
  qc = round(safe_numeric(qc), 2),
  prep = round(safe_numeric(prep), 2),
  DistinctAttractors = safe_numeric(DistinctAttractors),
  LargestBasinFraction = safe_numeric(LargestBasinFraction),
  BasinEntropy = safe_numeric(BasinEntropy)
)]
basin_system <- basin_system[is.finite(qc) & is.finite(prep)]
if (anyDuplicated(basin_system$SystemID)) stop("SystemID is not unique in the basin file.", call. = FALSE)

basin_metrics <- c("DistinctAttractors", "LargestBasinFraction", "BasinEntropy")
basin_long <- melt(basin_system, id.vars = setdiff(names(basin_system), basin_metrics), measure.vars = basin_metrics, variable.name = "Metric", value.name = "Value")
basin_cell_estimates <- mean_summary(basin_long, "Value", c("Metric", "prep", "qc"))
basin_endpoints <- endpoint_difference(basin_long, "Value", "qc", min(basin_long$qc), max(basin_long$qc), c("Metric", "prep"))
write_csv(basin_cell_estimates, "16_BASIN_CELL_ESTIMATES.csv")
write_csv(basin_endpoints, "17_BASIN_ENDPOINT_CONTRASTS.csv")

perturbation_trial <- read_input(PERTURBATION_FILE)
rename_standard(perturbation_trial, "SystemID", c("GraphID", "NetworkID", "BooleanSystemID"), TRUE)
rename_standard(perturbation_trial, "N", c("NumberOfNodes", "Nodes"), TRUE)
rename_standard(perturbation_trial, "d", c("d_input", "AvgDegreeInput", "DegreeInput"), TRUE)
rename_standard(perturbation_trial, "p", c("RewireProb", "p_rewire", "RewireProbability"), TRUE)
rename_standard(perturbation_trial, "qc", c("q_c", "CanalizingProbability", "CanalisingProbability"), TRUE)
rename_standard(perturbation_trial, "prep", c("p_rep", "InhibitionFraction"), TRUE)
rename_standard(perturbation_trial, "Returned", c("ReturnedToSameAttractor", "ReturnIndicator", "ReturnToSameAttractor"), TRUE)
rename_standard(perturbation_trial, "ReturnTime", c("TimeToReturn", "TauReturn"), TRUE)
rename_standard(perturbation_trial, "TimeObserved", c("ObservedReturnTime", "CensoredReturnTime"), FALSE)
rename_standard(perturbation_trial, "Event", c("ReturnEvent", "Status"), FALSE)
rename_standard(perturbation_trial, "ReturnCap", c("ReturnWindow", "CensoringTime"), TRUE)
rename_standard(perturbation_trial, "H1", c("HammingFractionT1", "MeanHammingFractionT1"), TRUE)
rename_standard(perturbation_trial, "H10", c("HammingFractionT10", "MeanHammingFractionT10"), TRUE)
rename_standard(perturbation_trial, "H25", c("HammingFractionT25", "MeanHammingFractionT25"), TRUE)
perturbation_trial[, `:=`(
  N = as.integer(round(safe_numeric(N))),
  d = as.integer(round(safe_numeric(d))),
  p = round(safe_numeric(p), 4),
  qc = round(safe_numeric(qc), 2),
  prep = round(safe_numeric(prep), 2),
  Returned = safe_logical(Returned),
  ReturnTime = safe_numeric(ReturnTime),
  ReturnCap = safe_numeric(ReturnCap),
  H1 = safe_numeric(H1),
  H10 = safe_numeric(H10),
  H25 = safe_numeric(H25)
)]
if (!"TimeObserved" %in% names(perturbation_trial)) perturbation_trial[, TimeObserved := ifelse(Returned, ReturnTime, ReturnCap)]
if (!"Event" %in% names(perturbation_trial)) perturbation_trial[, Event := as.integer(Returned)]
perturbation_trial[, `:=`(TimeObserved = safe_numeric(TimeObserved), Event = as.integer(safe_numeric(Event)))]
perturbation_trial <- perturbation_trial[!is.na(Returned) & is.finite(TimeObserved) & is.finite(qc) & is.finite(prep)]

perturbation_system <- perturbation_trial[, .(
  Trials = .N,
  ReturnedTrials = sum(Returned),
  ReturnProbability = mean(Returned),
  SuccessfulMeanReturnTime = if (any(Returned & is.finite(ReturnTime))) mean(ReturnTime[Returned & is.finite(ReturnTime)]) else NA_real_,
  SuccessfulMedianReturnTime = if (any(Returned & is.finite(ReturnTime))) median(ReturnTime[Returned & is.finite(ReturnTime)]) else NA_real_,
  RestrictedMeanReturnTime = mean(TimeObserved),
  HammingFractionT1 = mean(H1, na.rm = TRUE),
  HammingFractionT10 = mean(H10, na.rm = TRUE),
  HammingFractionT25 = mean(H25, na.rm = TRUE)
), by = .(SystemID, N, d, p, qc, prep)]

perturbation_metrics <- c("ReturnProbability", "SuccessfulMeanReturnTime", "SuccessfulMedianReturnTime", "RestrictedMeanReturnTime", "HammingFractionT1", "HammingFractionT10", "HammingFractionT25")
perturbation_long <- melt(perturbation_system, id.vars = setdiff(names(perturbation_system), perturbation_metrics), measure.vars = perturbation_metrics, variable.name = "Metric", value.name = "Value")
perturbation_cell_estimates <- mean_summary(perturbation_long, "Value", c("Metric", "prep", "qc"))
perturbation_endpoints <- endpoint_difference(perturbation_long, "Value", "qc", min(perturbation_long$qc), max(perturbation_long$qc), c("Metric", "prep"))

perturbation_trial[, `:=`(qc_factor = factor(sprintf("%.2f", qc), levels = sprintf("%.2f", sort(unique(qc)))), prep_factor = factor(sprintf("%.2f", prep), levels = sprintf("%.2f", sort(unique(prep)))))]
km_fit <- survival::survfit(survival::Surv(TimeObserved, Event) ~ qc_factor + prep_factor, data = perturbation_trial, conf.int = CI_LEVEL)
km_table <- as.data.table(summary(km_fit)$table, keep.rownames = "Stratum")
write_csv(perturbation_cell_estimates, "18_PERTURBATION_CELL_ESTIMATES.csv")
write_csv(km_table, "18_PERTURBATION_KAPLAN_MEIER_ESTIMATES.csv")
write_csv(perturbation_endpoints, "19_PERTURBATION_ENDPOINT_CONTRASTS.csv")

hamming_system <- read_input(HAMMING_FILE)
rename_standard(hamming_system, "SystemID", c("GraphID", "NetworkID"), TRUE)
rename_standard(hamming_system, "qc", c("q_c", "CanalizingProbability", "CanalisingProbability"), TRUE)
rename_standard(hamming_system, "prep", c("p_rep", "InhibitionFraction"), TRUE)
rename_standard(hamming_system, "Timestep", c("Time", "Step"), TRUE)
rename_standard(hamming_system, "MeanHammingFraction", c("HammingFraction", "MeanFraction"), TRUE)
hamming_system[, `:=`(qc = round(safe_numeric(qc), 2), prep = round(safe_numeric(prep), 2), Timestep = as.integer(safe_numeric(Timestep)), MeanHammingFraction = safe_numeric(MeanHammingFraction))]
hamming_audit <- hamming_system[, .(Systems = uniqueN(SystemID), Rows = .N, MeanHammingFraction = mean(MeanHammingFraction, na.rm = TRUE)), by = .(prep, qc, Timestep)]
write_csv(hamming_audit, "00_HAMMING_AUDIT.csv")

basin_audit <- basin_system[, .(Systems = .N, MinimumSampledInitialStates = if ("SampledInitialStates" %in% names(basin_system) && any(is.finite(SampledInitialStates))) min(SampledInitialStates[is.finite(SampledInitialStates)]) else NA_real_, MinimumFoundRate = if ("FoundRate" %in% names(basin_system) && any(is.finite(FoundRate))) min(FoundRate[is.finite(FoundRate)]) else NA_real_), by = .(prep, qc)]
perturbation_counts <- perturbation_trial[, .(Trials = .N, Returned = sum(Returned), Censored = sum(!Returned)), by = .(SystemID, prep, qc)]
perturbation_audit <- perturbation_counts[, .(Systems = .N, Trials = sum(Trials), Returned = sum(Returned), Censored = sum(Censored), MinimumTrialsPerSystem = min(Trials), MaximumTrialsPerSystem = max(Trials)), by = .(prep, qc)]
write_csv(basin_audit, "00_BASIN_AUDIT.csv")
write_csv(perturbation_audit, "00_PERTURBATION_AUDIT.csv")

key_numbers <- rbindlist(list(
  observed_period_endpoint[, .(Section = "5.1", SentenceKey = "Observed_qc_endpoint_period", Metric = "MeanLogPeriod", prep, Context = NA_character_, Estimate = Difference, SE = DifferenceSE, Lower = DifferenceLower, Upper = DifferenceUpper, Ratio = PeriodRatio, PercentChange, PValue, AdjustedPValue = NA_real_, Scale = "mean_log_period")],
  observed_fixed_endpoint[, .(Section = "5.1", SentenceKey = "Observed_qc_endpoint_fixed_point", Metric = "FixedPointRate", prep, Context = NA_character_, Estimate = Difference, SE = DifferenceSE, Lower = DifferenceLower, Upper = DifferenceUpper, Ratio = NA_real_, PercentChange = NA_real_, PValue, AdjustedPValue = NA_real_, Scale = "rate_difference")],
  qc_contrasts[, .(Section = "5.1", SentenceKey = "Adjusted_qc_contrast", Metric = paste0("qc_", Target, "_vs_", Reference), prep, Context = NA_character_, Estimate = DeltaY, SE, Lower = LowerDeltaY, Upper = UpperDeltaY, Ratio = PeriodRatio, PercentChange, PValue, AdjustedPValue, Scale = "mean_log_period")],
  prep_contrasts[, .(Section = "5.1", SentenceKey = "Adjusted_prep_contrast", Metric = paste0("prep_", Target, "_vs_", Reference), prep = Target, Context = paste0("qc=", qc), Estimate = DeltaY, SE, Lower = LowerDeltaY, Upper = UpperDeltaY, Ratio = PeriodRatio, PercentChange, PValue, AdjustedPValue, Scale = "mean_log_period")],
  qc_endpoint_by_p_prep[, .(Section = "5.2", SentenceKey = "Qc_endpoint_at_each_p", Metric = "qc_max_vs_min", prep, Context = paste0("p=", p), Estimate = DeltaY, SE, Lower = LowerDeltaY, Upper = UpperDeltaY, Ratio = PeriodRatio, PercentChange, PValue, AdjustedPValue = NA_real_, Scale = "mean_log_period")],
  p_endpoints[, .(Section = "5.2", SentenceKey = "Rewiring_endpoint", Metric = "p_max_vs_min", prep, Context = paste0("qc=", qc), Estimate = DeltaY, SE, Lower = LowerDeltaY, Upper = UpperDeltaY, Ratio = PeriodRatio, PercentChange, PValue, AdjustedPValue = NA_real_, Scale = "mean_log_period")],
  kappa_qc_prep_endpoints[N == max(N_levels), .(Section = "5.3", SentenceKey = "Kappa_specific_qc_endpoint_at_N100", Metric = "qc_max_vs_min", prep, Context = paste0("kappa=", kappa), Estimate = DeltaY, SE, Lower = LowerDeltaY, Upper = UpperDeltaY, Ratio = PeriodRatio, PercentChange, PValue, AdjustedPValue = NA_real_, Scale = "mean_log_period")],
  kappa_prep_contrasts[, .(Section = "5.3", SentenceKey = "Kappa_specific_prep_contrast_at_N100", Metric = paste0("prep_", Target, "_vs_", Reference), prep = Target, Context = paste0("N=", max(N_levels), ";kappa=", kappa), Estimate = DeltaY, SE, Lower = LowerDeltaY, Upper = UpperDeltaY, Ratio = PeriodRatio, PercentChange, PValue, AdjustedPValue, Scale = "mean_log_period")],
  basin_endpoints[, .(Section = "5.4", SentenceKey = "Basin_qc_endpoint", Metric = as.character(Metric), prep, Context = NA_character_, Estimate = Difference, SE = DifferenceSE, Lower = DifferenceLower, Upper = DifferenceUpper, Ratio = NA_real_, PercentChange = NA_real_, PValue, AdjustedPValue = NA_real_, Scale = "metric_difference")],
  perturbation_endpoints[, .(Section = "5.5", SentenceKey = "Perturbation_qc_endpoint", Metric = as.character(Metric), prep, Context = NA_character_, Estimate = Difference, SE = DifferenceSE, Lower = DifferenceLower, Upper = DifferenceUpper, Ratio = NA_real_, PercentChange = NA_real_, PValue, AdjustedPValue = NA_real_, Scale = "metric_difference")]
), fill = TRUE)
write_csv(key_numbers, "20_RESULTS_KEY_NUMBERS.csv")

direction_check <- copy(key_numbers)
direction_check[, Direction := fifelse(Lower > 0, "Increase", fifelse(Upper < 0, "Decrease", "Uncertain"))]
direction_check[, Supported := Direction != "Uncertain"]
write_csv(direction_check, "21_RESULTS_DIRECTION_CHECK.csv")


# =============================================================================
# FINAL PUBLICATION FIGURES
# Figure numbering and file names below match the Paper 2 manuscript.
# =============================================================================

rule_panel <- function(panel_title, mode) {
  source_labels <- switch(
    mode,
    threshold = c("s[i](t)", "s[m](t)", "s[k](t)"),
    active = c("s[i](t)", "s[r(j)](t)==1", "s[k](t)"),
    inactive = c("s[i](t)", "s[r(j)](t)==0", "s[k](t)")
  )

  output_label <- if (mode == "active") "b[j]" else "s[j](t+1)"

  source_nodes <- data.table(
    x = c(0, 0, 0),
    y = c(1.20, 0.60, 0.00),
    label = source_labels,
    Selected = c(FALSE, mode != "threshold", FALSE)
  )

  target_node <- data.table(x = 2.20, y = 0.60, label = "Node j")
  output_node <- data.table(x = 4.30, y = 0.60, label = output_label)

  edge_data <- data.table(
    x = source_nodes$x + 0.34,
    y = source_nodes$y,
    xend = target_node$x - 0.43,
    yend = target_node$y,
    Selected = source_nodes$Selected
  )

  edge_data[, `:=`(
    Colour = fifelse(
      mode == "active" & Selected,
      "#D55E00",
      fifelse(
        mode == "active" & !Selected,
        "#B8B8B8",
        fifelse(mode == "inactive" & Selected, "#A0A0A0", "#3C3C3C")
      )
    ),
    Width = fifelse(Selected, 1.15, fifelse(mode == "active", 0.48, 0.68)),
    LineType = fifelse(mode == "inactive" & Selected, "dashed", "solid")
  )]

  condition_labels <- switch(
    mode,
    threshold = "z[j]==0",
    active = "z[j]==1*','~~s[r(j)](t)==1",
    inactive = "z[j]==1*','~~s[r(j)](t)==0"
  )

  equation_labels <- switch(
    mode,
    threshold = "I[j](t)==sum(a[ij]*s[i](t),i)*','~~s[j](t+1)==bold(1)*group('[',I[j](t)>theta[j],']')",
    active = "s[j](t+1)==b[j]*','~~b[j]==1~'(activation),'~~b[j]==0~'(inhibition)'",
    inactive = "I[j](t)==sum(a[ij]*s[i](t),i)*','~~s[j](t+1)==bold(1)*group('[',I[j](t)>theta[j]^fb,']')"
  )

  condition_data <- data.table(x = 2.20, y = 1.52, label = condition_labels)
  equation_data <- data.table(x = 2.20, y = -0.42, label = equation_labels)

  ggplot() +
    geom_segment(
      data = edge_data,
      aes(
        x = x, y = y, xend = xend, yend = yend,
        colour = Colour, linewidth = Width, linetype = LineType
      ),
      arrow = grid::arrow(type = "closed", length = grid::unit(1.8, "mm")),
      show.legend = FALSE
    ) +
    scale_colour_identity() +
    scale_linewidth_identity() +
    scale_linetype_identity() +
    geom_segment(
      aes(
        x = target_node$x + 0.43,
        y = target_node$y,
        xend = output_node$x - 0.43,
        yend = output_node$y
      ),
      linewidth = 0.78,
      colour = if (mode == "active") "#D55E00" else "#222222",
      arrow = grid::arrow(type = "closed", length = grid::unit(1.8, "mm"))
    ) +
    geom_label(
      data = source_nodes,
      aes(x, y, label = label),
      parse = TRUE,
      size = 4.00,
      linewidth = 0.30,
      label.padding = grid::unit(1.45, "mm"),
      fill = "white"
    ) +
    geom_label(
      data = target_node,
      aes(x, y, label = label),
      size = 4.00,
      linewidth = 0.34,
      label.padding = grid::unit(1.65, "mm"),
      fill = "#E8F1F8"
    ) +
    geom_label(
      data = output_node,
      aes(x, y, label = label),
      parse = TRUE,
      size = 4.00,
      linewidth = 0.34,
      label.padding = grid::unit(1.65, "mm"),
      fill = if (mode == "active") "#FCE8D5" else "#F3F3F3"
    ) +
    geom_text(
      data = condition_data,
      aes(x, y, label = label),
      parse = TRUE,
      size = 4.40,
      fontface = "bold"
    ) +
    geom_text(
      data = equation_data,
      aes(x, y, label = label),
      parse = TRUE,
      size = 4.15
    ) +
    {
      if (mode == "inactive") {
        annotate(
          "text",
          x = 1.08,
          y = 0.60,
          label = "×",
          size = 6.0,
          colour = "#777777",
          fontface = "bold"
        )
      } else {
        NULL
      }
    } +
    coord_cartesian(
      xlim = c(-0.72, 5.02),
      ylim = c(-0.74, 1.88),
      clip = "off"
    ) +
    labs(title = panel_title) +
    theme_void(base_size = 10.5, base_family = "sans") +
    theme(
      plot.title = element_text(
        size = 12.8,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 3)
      ),
      plot.margin = margin(5, 6, 5, 6)
    )
}

# -----------------------------------------------------------------------------
# Figure 1: update-rule scenarios
# -----------------------------------------------------------------------------

scenario_1 <- rule_panel("Signed-threshold rule", "threshold")
scenario_2 <- rule_panel("Canalising regulator active", "active")
scenario_3 <- rule_panel("Canalising regulator inactive", "inactive")

figure_1 <- (
  scenario_1 /
    (scenario_2 + scenario_3)
) +
  plot_layout(heights = c(0.95, 1.05)) +
  plot_annotation(tag_levels = "A")

save_figure(
  figure_1,
  "FIGURE_01_UPDATE_RULE_SCENARIOS",
  height = 6.5,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Figure 2: observed period, fixed points, and adjusted rule effects
# -----------------------------------------------------------------------------

observed[, prep_factor := factor(
  sprintf("%.2f", prep),
  levels = names(prep_palette)
)]

figure_2a <- ggplot(
  observed,
  aes(
    qc,
    GeometricMeanPeriod,
    colour = prep_factor,
    fill = prep_factor,
    linetype = prep_factor,
    shape = prep_factor,
    group = prep_factor
  )
) +
  geom_ribbon(
    aes(
      ymin = LowerGeometricMeanPeriod,
      ymax = UpperGeometricMeanPeriod
    ),
    alpha = 0.11,
    colour = NA,
    linetype = 0
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.1) +
  scale_colour_manual(
    values = prep_palette,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = prep_palette,
    guide = "none",
    drop = TRUE
  ) +
  scale_linetype_manual(
    values = prep_linetypes,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_shape_manual(
    values = c(
      "0.10" = 16,
      "0.24" = 17,
      "0.30" = 15,
      "0.41" = 18,
      "0.50" = 8
    ),
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_x_continuous(
    breaks = qc_levels,
    labels = format_qc_axis(qc_levels)
  ) +
  labs(
    x = expression(q[c]),
    y = "Observed geometric mean period"
  ) +
  paper2_theme()

figure_2b <- ggplot(
  observed,
  aes(
    qc,
    FixedPointRate,
    colour = prep_factor,
    fill = prep_factor,
    linetype = prep_factor,
    shape = prep_factor,
    group = prep_factor
  )
) +
  geom_ribbon(
    aes(
      ymin = pmax(0, LowerFixedPointRate),
      ymax = pmin(1, UpperFixedPointRate)
    ),
    alpha = 0.11,
    colour = NA,
    linetype = 0
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.1) +
  scale_colour_manual(
    values = prep_palette,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = prep_palette,
    guide = "none",
    drop = TRUE
  ) +
  scale_linetype_manual(
    values = prep_linetypes,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_shape_manual(
    values = c(
      "0.10" = 16,
      "0.24" = 17,
      "0.30" = 15,
      "0.41" = 18,
      "0.50" = 8
    ),
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_x_continuous(
    breaks = qc_levels,
    labels = format_qc_axis(qc_levels)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, NA)
  ) +
  labs(
    x = expression(q[c]),
    y = "Fixed-point rate"
  ) +
  paper2_theme()

qc_prep_emm[, prep_factor := factor(
  sprintf("%.2f", prep),
  levels = names(prep_palette)
)]

figure_2c <- ggplot(
  qc_prep_emm,
  aes(
    qc,
    EstimateY,
    colour = prep_factor,
    fill = prep_factor,
    linetype = prep_factor,
    shape = prep_factor,
    group = prep_factor
  )
) +
  geom_ribbon(
    aes(ymin = LowerY, ymax = UpperY),
    alpha = 0.10,
    colour = NA,
    linetype = 0
  ) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.0) +
  scale_colour_manual(
    values = prep_palette,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = prep_palette,
    guide = "none",
    drop = TRUE
  ) +
  scale_linetype_manual(
    values = prep_linetypes,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_shape_manual(
    values = c(
      "0.10" = 16,
      "0.24" = 17,
      "0.30" = 15,
      "0.41" = 18,
      "0.50" = 8
    ),
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_x_continuous(
    breaks = qc_levels,
    labels = format_qc_axis(qc_levels)
  ) +
  labs(
    x = expression(q[c]),
    y = "Adjusted mean log period"
  ) +
  paper2_theme()

figure_2 <- wrap_plots(
  figure_2a,
  figure_2b,
  figure_2c,
  ncol = 3,
  guides = "collect"
) +
  plot_annotation(tag_levels = "A")

save_figure(
  figure_2,
  "FIGURE_02_OBSERVED_AND_ADJUSTED_RULE_EFFECTS",
  height = 4.25,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Figure 3: rewiring effect
# Manuscript uses q_c = 0, 0.50, 1 and p_rep = 0.10, 0.24, 0.50.
# -----------------------------------------------------------------------------

selected_qc <- qc_levels[qc_levels %in% c(0, 0.5, 1)]
selected_prep <- prep_levels[prep_levels %in% c(0.10, 0.24, 0.50)]

figure_3_data <- p_qc_prep_emm[
  qc %in% selected_qc &
    prep %in% selected_prep
]

figure_3_data[, `:=`(
  p_factor = factor(
    sprintf("%.2f", p),
    levels = sprintf("%.2f", p_levels)
  ),
  qc_factor = factor(
    sprintf("%.2f", qc),
    levels = names(qc_palette)
  ),
  prep_strip = factor(
    sprintf("p[rep] == %.2f", prep),
    levels = sprintf("p[rep] == %.2f", selected_prep)
  )
)]

figure_3 <- ggplot(
  figure_3_data,
  aes(
    p_factor,
    EstimateY,
    colour = qc_factor,
    shape = qc_factor,
    group = qc_factor
  )
) +
  geom_errorbar(
    aes(ymin = LowerY, ymax = UpperY),
    width = 0.10,
    linewidth = 0.48,
    position = position_dodge(width = 0.12)
  ) +
  geom_line(
    linewidth = 0.82,
    position = position_dodge(width = 0.12)
  ) +
  geom_point(
    size = 2.0,
    position = position_dodge(width = 0.12)
  ) +
  facet_wrap(
    ~prep_strip,
    ncol = 3,
    labeller = label_parsed
  ) +
  scale_colour_manual(
    values = qc_palette,
    name = expression(q[c]),
    drop = TRUE
  ) +
  scale_shape_manual(
    values = qc_shapes,
    name = expression(q[c]),
    drop = TRUE
  ) +
  labs(
    x = expression("Rewiring probability "*p),
    y = "Adjusted mean log period"
  ) +
  paper2_theme(9.5)

save_figure(
  figure_3,
  "FIGURE_03_REWIRING_EFFECT",
  height = 4.8,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Figure 4: network size and kappa
# -----------------------------------------------------------------------------

N_kappa_emm[, kappa_factor := factor(
  kappa,
  levels = kappa_levels
)]

kappa_linetypes <- setNames(
  rep(series_linetypes, length.out = length(kappa_levels)),
  as.character(kappa_levels)
)

kappa_shapes <- setNames(
  rep(series_shapes, length.out = length(kappa_levels)),
  as.character(kappa_levels)
)

figure_4b_data <- N_kappa_qc_emm[N == max(N_levels)]
figure_4c_data <- N_kappa_prep_emm[N == max(N_levels)]

figure_4b_data[, qc_factor := factor(
  sprintf("%.2f", qc),
  levels = names(qc_palette)
)]

figure_4c_data[, prep_factor := factor(
  sprintf("%.2f", prep),
  levels = names(prep_palette)
)]

figure_4a <- ggplot(
  N_kappa_emm,
  aes(
    N,
    EstimateY,
    colour = kappa_factor,
    fill = kappa_factor,
    linetype = kappa_factor,
    shape = kappa_factor,
    group = kappa_factor
  )
) +
  geom_ribbon(
    aes(ymin = LowerY, ymax = UpperY),
    alpha = 0.07,
    colour = NA,
    linetype = 0
  ) +
  geom_line(linewidth = 0.78) +
  geom_point(size = 1.8) +
  scale_colour_viridis_d(
    option = "D",
    end = 0.88,
    name = expression(kappa)
  ) +
  scale_fill_viridis_d(
    option = "D",
    end = 0.88,
    guide = "none"
  ) +
  scale_linetype_manual(
    values = kappa_linetypes,
    name = expression(kappa)
  ) +
  scale_shape_manual(
    values = kappa_shapes,
    name = expression(kappa)
  ) +
  scale_x_continuous(breaks = N_levels) +
  labs(
    x = expression(N),
    y = "Adjusted mean log period"
  ) +
  paper2_theme()

figure_4b <- ggplot(
  figure_4b_data,
  aes(
    kappa,
    EstimateY,
    colour = qc_factor,
    fill = qc_factor,
    shape = qc_factor,
    group = qc_factor
  )
) +
  geom_ribbon(
    aes(ymin = LowerY, ymax = UpperY),
    alpha = 0.09,
    colour = NA
  ) +
  geom_line(linewidth = 0.78) +
  geom_point(size = 1.9) +
  scale_colour_manual(
    values = qc_palette,
    name = expression(q[c]),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = qc_palette,
    guide = "none",
    drop = TRUE
  ) +
  scale_shape_manual(
    values = qc_shapes,
    name = expression(q[c]),
    drop = TRUE
  ) +
  scale_x_continuous(breaks = kappa_levels) +
  labs(
    x = expression(kappa),
    y = "Adjusted mean log period",
    subtitle = paste0("N = ", max(N_levels))
  ) +
  paper2_theme() +
  theme(
    plot.subtitle = element_text(
      size = 9.5,
      hjust = 0.5
    )
  )

figure_4c <- ggplot(
  figure_4c_data,
  aes(
    kappa,
    EstimateY,
    colour = prep_factor,
    fill = prep_factor,
    linetype = prep_factor,
    shape = prep_factor,
    group = prep_factor
  )
) +
  geom_ribbon(
    aes(ymin = LowerY, ymax = UpperY),
    alpha = 0.09,
    colour = NA,
    linetype = 0
  ) +
  geom_line(linewidth = 0.78) +
  geom_point(size = 1.9) +
  scale_colour_manual(
    values = prep_palette,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = prep_palette,
    guide = "none",
    drop = TRUE
  ) +
  scale_linetype_manual(
    values = prep_linetypes,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_shape_manual(
    values = c(
      "0.10" = 16,
      "0.24" = 17,
      "0.30" = 15,
      "0.41" = 18,
      "0.50" = 8
    ),
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_x_continuous(breaks = kappa_levels) +
  labs(
    x = expression(kappa),
    y = "Adjusted mean log period",
    subtitle = paste0("N = ", max(N_levels))
  ) +
  paper2_theme() +
  theme(
    plot.subtitle = element_text(
      size = 9.5,
      hjust = 0.5
    )
  )

figure_4 <- (
  figure_4a /
    (figure_4b + figure_4c)
) +
  plot_annotation(tag_levels = "A")

save_figure(
  figure_4,
  "FIGURE_04_NETWORK_SIZE_AND_KAPPA",
  height = 7.2,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Shared line-panel helper for Figures 5, 6 and S1
# -----------------------------------------------------------------------------

line_summary_panel <- function(
    DT,
    metric_name,
    y_label,
    panel_tag = NULL,
    legend_position = "right") {

  local <- copy(DT[Metric == metric_name])

  local[, prep_factor := factor(
    sprintf("%.2f", prep),
    levels = names(prep_palette)
  )]

  plot_object <- ggplot(
    local,
    aes(
      qc,
      Estimate,
      colour = prep_factor,
      fill = prep_factor,
      linetype = prep_factor,
      shape = prep_factor,
      group = prep_factor
    )
  ) +
    geom_ribbon(
      aes(ymin = Lower, ymax = Upper),
      alpha = 0.11,
      colour = NA,
      linetype = 0
    ) +
    geom_line(linewidth = 0.82) +
    geom_point(size = 1.9) +
    scale_colour_manual(
      values = prep_palette,
      name = expression(p[rep]),
      drop = TRUE
    ) +
    scale_fill_manual(
      values = prep_palette,
      guide = "none",
      drop = TRUE
    ) +
    scale_linetype_manual(
      values = prep_linetypes,
      name = expression(p[rep]),
      drop = TRUE
    ) +
    scale_shape_manual(
      values = c(
        "0.10" = 16,
        "0.24" = 17,
        "0.30" = 15,
        "0.41" = 18,
        "0.50" = 8
      ),
      name = expression(p[rep]),
      drop = TRUE
    ) +
    scale_x_continuous(
      breaks = sort(unique(local$qc)),
      labels = format_qc_axis(sort(unique(local$qc)))
    ) +
    labs(
      x = expression(q[c]),
      y = y_label,
      tag = panel_tag
    ) +
    paper2_theme(9.5) +
    theme(
      legend.position = legend_position
    )

  if (!is.null(panel_tag)) {
    plot_object <- plot_object +
      theme(
        # Keep panel letters inside the upper-left of each plotting panel.
        # This prevents C and D from crowding long vertical y-axis labels.
        plot.tag.position = c(0.055, 0.975),
        plot.tag = element_text(
          face = "bold",
          size = 10.5,
          hjust = 0,
          vjust = 1
        ),
        plot.margin = margin(
          t = 7,
          r = 5,
          b = 5,
          l = 5
        )
      )
  }

  plot_object
}

# -----------------------------------------------------------------------------
# Figure 5: basin organisation
# -----------------------------------------------------------------------------

figure_5 <- wrap_plots(
  line_summary_panel(
    basin_cell_estimates,
    "DistinctAttractors",
    "Number of attractors reached"
  ),
  line_summary_panel(
    basin_cell_estimates,
    "LargestBasinFraction",
    "Largest observed basin fraction"
  ),
  line_summary_panel(
    basin_cell_estimates,
    "BasinEntropy",
    "Observed basin entropy"
  ),
  ncol = 3,
  guides = "collect"
) +
  plot_annotation(tag_levels = "A")

save_figure(
  figure_5,
  "FIGURE_05_BASIN_ORGANISATION",
  height = 4.2,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Figure 6: perturbation recovery
# -----------------------------------------------------------------------------

figure_6 <- wrap_plots(
  line_summary_panel(
    perturbation_cell_estimates,
    "ReturnProbability",
    "Return probability"
  ),
  line_summary_panel(
    perturbation_cell_estimates,
    "RestrictedMeanReturnTime",
    "Restricted mean return time (updates)"
  ),
  line_summary_panel(
    perturbation_cell_estimates,
    "HammingFractionT25",
    "Hamming fraction at 25 updates"
  ),
  ncol = 3,
  guides = "collect"
) +
  plot_annotation(tag_levels = "A")

save_figure(
  figure_6,
  "FIGURE_06_PERTURBATION_RECOVERY",
  height = 4.2,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Figure S1: six supplementary perturbation/recovery measures
# Corrected panel-tag placement: C and D no longer overlap y-axis labels.
# -----------------------------------------------------------------------------

s1_metrics <- c(
  "ReturnProbability",
  "RestrictedMeanReturnTime",
  "SuccessfulMeanReturnTime",
  "SuccessfulMedianReturnTime",
  "HammingFractionT1",
  "HammingFractionT10"
)

s1_labels <- c(
  ReturnProbability =
    "Return probability",
  RestrictedMeanReturnTime =
    "Restricted mean return time (updates)",
  SuccessfulMeanReturnTime =
    "Mean return time among returned trials (updates)",
  SuccessfulMedianReturnTime =
    "Median return time among returned trials (updates)",
  HammingFractionT1 =
    "Hamming fraction at 1 update",
  HammingFractionT10 =
    "Hamming fraction at 10 updates"
)

figure_s1_panels <- Map(
  function(metric_name, panel_tag) {
    line_summary_panel(
      perturbation_cell_estimates,
      metric_name,
      s1_labels[[metric_name]],
      panel_tag = panel_tag,
      legend_position = "bottom"
    )
  },
  s1_metrics,
  LETTERS[1:6]
)

figure_s1 <- wrap_plots(
  plotlist = figure_s1_panels,
  ncol = 3,
  guides = "collect"
)

save_figure(
  figure_s1,
  "FIGURE_S1_PERTURBATION_MEASURES",
  height = 7.0,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Figure S2: representative rewiring probabilities
# -----------------------------------------------------------------------------

representative_p <- p_levels[
  vapply(
    p_levels,
    function(x) any(abs(x - c(0.01, 0.20, 0.60)) < 1e-10),
    logical(1)
  )
]

figure_s2_data <- copy(
  p_qc_prep_emm[
    p %in% representative_p
  ]
)

figure_s2_data[, `:=`(
  p_strip = factor(
    sprintf("p == %.2f", p),
    levels = sprintf("p == %.2f", c(0.01, 0.20, 0.60))
  ),
  prep_factor = factor(
    sprintf("%.2f", prep),
    levels = names(prep_palette)
  )
)]

figure_s2 <- ggplot(
  figure_s2_data,
  aes(
    qc,
    EstimateY,
    colour = prep_factor,
    fill = prep_factor,
    linetype = prep_factor,
    shape = prep_factor,
    group = prep_factor
  )
) +
  geom_ribbon(
    aes(ymin = LowerY, ymax = UpperY),
    alpha = 0.09,
    colour = NA,
    linetype = 0
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.9) +
  facet_wrap(
    ~p_strip,
    ncol = 3,
    labeller = label_parsed
  ) +
  scale_colour_manual(
    values = prep_palette,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_fill_manual(
    values = prep_palette,
    guide = "none",
    drop = TRUE
  ) +
  scale_linetype_manual(
    values = prep_linetypes,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_shape_manual(
    values = c(
      "0.10" = 16,
      "0.24" = 17,
      "0.30" = 15,
      "0.41" = 18,
      "0.50" = 8
    ),
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_x_continuous(
    breaks = qc_levels,
    labels = format_qc_axis(qc_levels)
  ) +
  labs(
    x = expression(q[c]),
    y = "Adjusted mean log period"
  ) +
  paper2_theme(9.5)

save_figure(
  figure_s2,
  "FIGURE_S2_RULE_EFFECTS_AT_REPRESENTATIVE_P",
  height = 4.8,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Figures S3A and S3B
# -----------------------------------------------------------------------------

pair_plot <- function(
    DT,
    x_variable,
    colour_variable,
    x_label,
    colour_label,
    panel_title) {

  local <- copy(DT)
  colour_values <- sort(unique(safe_numeric(local[[colour_variable]])))

  colour_levels <- if (colour_variable == "kappa") {
    sprintf("%d", as.integer(colour_values))
  } else {
    sprintf("%.2f", colour_values)
  }

  local[, colour_factor := factor(
    if (colour_variable == "kappa") {
      sprintf("%d", as.integer(safe_numeric(get(colour_variable))))
    } else {
      sprintf("%.2f", safe_numeric(get(colour_variable)))
    },
    levels = colour_levels
  )]

  local_linetypes <- setNames(
    rep(series_linetypes, length.out = length(colour_levels)),
    colour_levels
  )

  local_shapes <- setNames(
    rep(series_shapes, length.out = length(colour_levels)),
    colour_levels
  )

  ggplot(
    local,
    aes(
      x = .data[[x_variable]],
      y = EstimateY,
      colour = colour_factor,
      fill = colour_factor,
      linetype = colour_factor,
      shape = colour_factor,
      group = colour_factor
    )
  ) +
    geom_ribbon(
      aes(ymin = LowerY, ymax = UpperY),
      alpha = 0.07,
      colour = NA,
      linetype = 0
    ) +
    geom_line(linewidth = 0.75) +
    geom_point(size = 1.7) +
    scale_colour_viridis_d(
      option = "D",
      end = 0.88,
      name = colour_label
    ) +
    scale_fill_viridis_d(
      option = "D",
      end = 0.88,
      guide = "none"
    ) +
    scale_linetype_manual(
      values = local_linetypes,
      name = colour_label
    ) +
    scale_shape_manual(
      values = local_shapes,
      name = colour_label
    ) +
    labs(
      x = x_label,
      y = "Adjusted mean log period",
      subtitle = panel_title
    ) +
    paper2_theme(9.5) +
    theme(
      plot.subtitle = element_text(
        size = 9.5,
        face = "bold",
        hjust = 0.5
      )
    )
}

figure_s3a <- wrap_plots(
  pair_plot(
    N_kappa_emm,
    "N",
    "kappa",
    expression(N),
    expression(kappa),
    expression(N~"x"~kappa~"on observed support")
  ),
  pair_plot(
    N_p_emm,
    "N",
    "p",
    expression(N),
    expression(p),
    bquote(N~"x"~p~"at"~kappa == .(min(kappa_levels)))
  ),
  pair_plot(
    N_qc_emm,
    "N",
    "qc",
    expression(N),
    expression(q[c]),
    bquote(N~"x"~q[c]~"at"~kappa == .(min(kappa_levels)))
  ),
  pair_plot(
    N_prep_emm,
    "N",
    "prep",
    expression(N),
    expression(p[rep]),
    bquote(N~"x"~p[rep]~"at"~kappa == .(min(kappa_levels)))
  ),
  ncol = 2
) +
  plot_annotation(tag_levels = "A")

save_figure(
  figure_s3a,
  "FIGURE_S3A_NETWORK_SIZE_PAIRWISE_VIEWS",
  height = 7.0,
  width = 7.8
)

figure_s3b <- wrap_plots(
  pair_plot(
    kappa_p_emm,
    "kappa",
    "p",
    expression(kappa),
    expression(p),
    bquote(kappa~"x"~p~"at"~N == .(max(N_levels)))
  ),
  pair_plot(
    kappa_qc_emm,
    "kappa",
    "qc",
    expression(kappa),
    expression(q[c]),
    bquote(kappa~"x"~q[c]~"at"~N == .(max(N_levels)))
  ),
  pair_plot(
    kappa_prep_emm,
    "kappa",
    "prep",
    expression(kappa),
    expression(p[rep]),
    bquote(kappa~"x"~p[rep]~"at"~N == .(max(N_levels)))
  ),
  ncol = 2
) +
  plot_annotation(tag_levels = "A")

save_figure(
  figure_s3b,
  "FIGURE_S3B_KAPPA_PAIRWISE_VIEWS",
  height = 6.7,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Figure S4: kappa x q_c x p_rep interaction at N = 100
# -----------------------------------------------------------------------------

figure_s4_data <- kappa_qc_prep_endpoints[
  N == max(N_levels)
]

figure_s4_data[, prep_factor := factor(
  sprintf("%.2f", prep),
  levels = names(prep_palette)
)]

figure_s4 <- ggplot(
  figure_s4_data,
  aes(
    kappa,
    DeltaY,
    colour = prep_factor,
    linetype = prep_factor,
    shape = prep_factor,
    group = prep_factor
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = 2,
    linewidth = 0.45,
    colour = "#555555"
  ) +
  geom_errorbar(
    aes(ymin = LowerDeltaY, ymax = UpperDeltaY),
    width = 0.10,
    linewidth = 0.55,
    position = position_dodge(width = 0.08)
  ) +
  geom_line(
    linewidth = 0.82,
    position = position_dodge(width = 0.08)
  ) +
  geom_point(
    size = 2.1,
    position = position_dodge(width = 0.08)
  ) +
  scale_colour_manual(
    values = prep_palette,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_linetype_manual(
    values = prep_linetypes,
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_shape_manual(
    values = c(
      "0.10" = 16,
      "0.24" = 17,
      "0.30" = 15,
      "0.41" = 18,
      "0.50" = 8
    ),
    name = expression(p[rep]),
    drop = TRUE
  ) +
  scale_x_continuous(breaks = kappa_levels) +
  labs(
    title = paste0(
      "Connectivity-dependent canalisation effect at N = ",
      max(N_levels)
    ),
    x = expression(kappa),
    y = expression(
      Delta*" mean log period: "*q[c]*" = 1 minus "*q[c]*" = 0"
    )
  ) +
  paper2_theme() +
  theme(
    plot.title = element_text(
      size = 10,
      hjust = 0.5,
      face = "plain"
    )
  )

save_figure(
  figure_s4,
  "FIGURE_S4_KAPPA_QC_PREP_INTERACTION",
  height = 4.4,
  width = 7.8
)

# -----------------------------------------------------------------------------
# Figure source-data files
# -----------------------------------------------------------------------------

figure_1_source <- data.table(
  Panel = c("A", "B", "C"),
  Scenario = c(
    "Signed-threshold rule",
    "Canalising regulator active",
    "Canalising regulator inactive"
  ),
  Condition = c(
    "z_j = 0",
    "z_j = 1 and s_r(j)(t) = 1",
    "z_j = 1 and s_r(j)(t) = 0"
  ),
  Update = c(
    "Signed threshold using all inputs",
    "s_j(t+1) = b_j",
    "Fallback signed threshold when the dominant regulator is inactive"
  ),
  Interpretation = c(
    "All inputs contribute",
    "The selected edge sign fixes the output",
    "The dominant regulator contributes zero and the fallback threshold is used"
  )
)

write_csv(
  figure_1_source,
  "FIGURE_01_SOURCE_DEFINITIONS.csv"
)

write_csv(
  rbindlist(
    list(
      copy(observed)[, Panel := "A_B_observed"],
      copy(qc_prep_emm)[, Panel := "C_adjusted"]
    ),
    fill = TRUE
  ),
  "FIGURE_02_SOURCE_DATA.csv"
)

write_csv(
  figure_3_data,
  "FIGURE_03_SOURCE_DATA.csv"
)

write_csv(
  rbindlist(
    list(
      copy(N_kappa_emm)[, Panel := "A_N_by_kappa"],
      copy(figure_4b_data)[, Panel := "B_kappa_by_qc_at_Nmax"],
      copy(figure_4c_data)[, Panel := "C_kappa_by_prep_at_Nmax"]
    ),
    fill = TRUE
  ),
  "FIGURE_04_SOURCE_DATA.csv"
)

write_csv(
  basin_cell_estimates,
  "FIGURE_05_SOURCE_DATA.csv"
)

write_csv(
  perturbation_cell_estimates[
    Metric %in% c(
      "ReturnProbability",
      "RestrictedMeanReturnTime",
      "HammingFractionT25"
    )
  ],
  "FIGURE_06_SOURCE_DATA.csv"
)

write_csv(
  perturbation_cell_estimates[
    Metric %in% s1_metrics
  ],
  "FIGURE_S1_SOURCE_DATA.csv"
)

write_csv(
  figure_s2_data,
  "FIGURE_S2_SOURCE_DATA.csv"
)

write_csv(
  rbindlist(
    list(
      copy(N_kappa_emm)[, Panel := "A_N_by_kappa"],
      copy(N_p_emm)[, Panel := "B_N_by_p"],
      copy(N_qc_emm)[, Panel := "C_N_by_qc"],
      copy(N_prep_emm)[, Panel := "D_N_by_prep"]
    ),
    fill = TRUE
  ),
  "FIGURE_S3A_SOURCE_DATA.csv"
)

write_csv(
  rbindlist(
    list(
      copy(kappa_p_emm)[, Panel := "A_kappa_by_p"],
      copy(kappa_qc_emm)[, Panel := "B_kappa_by_qc"],
      copy(kappa_prep_emm)[, Panel := "C_kappa_by_prep"]
    ),
    fill = TRUE
  ),
  "FIGURE_S3B_SOURCE_DATA.csv"
)

write_csv(
  figure_s4_data,
  "FIGURE_S4_SOURCE_DATA.csv"
)

write_latex(model_comparison, "TABLE_MODEL_COMPARISON.tex")
write_latex(model_performance, "TABLE_PRIMARY_MODEL_PERFORMANCE.tex")
write_latex(observed, "TABLE_OBSERVED_PERIOD_FIXED_POINT.tex")
write_latex(qc_contrasts, "TABLE_ADJUSTED_QC_CONTRASTS.tex")
write_latex(prep_contrasts, "TABLE_ADJUSTED_PREP_CONTRASTS.tex")
write_latex(p_endpoints, "TABLE_REWIRING_CONTRASTS.tex")
write_latex(kappa_qc_prep_endpoints, "TABLE_KAPPA_QC_PREP_CONTRASTS.tex")
write_latex(kappa_prep_contrasts, "TABLE_KAPPA_PREP_CONTRASTS.tex")
write_latex(basin_endpoints, "TABLE_BASIN_ENDPOINT_CONTRASTS.tex")
write_latex(perturbation_endpoints, "TABLE_PERTURBATION_ENDPOINT_CONTRASTS.tex")

figure_manifest <- data.table(
  Figure = c(
    "Figure 1",
    "Figure 2",
    "Figure 3",
    "Figure 4",
    "Figure 5",
    "Figure 6",
    "Figure S1",
    "Figure S2",
    "Figure S3A",
    "Figure S3B",
    "Figure S4"
  ),
  FileStem = c(
    "FIGURE_01_UPDATE_RULE_SCENARIOS",
    "FIGURE_02_OBSERVED_AND_ADJUSTED_RULE_EFFECTS",
    "FIGURE_03_REWIRING_EFFECT",
    "FIGURE_04_NETWORK_SIZE_AND_KAPPA",
    "FIGURE_05_BASIN_ORGANISATION",
    "FIGURE_06_PERTURBATION_RECOVERY",
    "FIGURE_S1_PERTURBATION_MEASURES",
    "FIGURE_S2_RULE_EFFECTS_AT_REPRESENTATIVE_P",
    "FIGURE_S3A_NETWORK_SIZE_PAIRWISE_VIEWS",
    "FIGURE_S3B_KAPPA_PAIRWISE_VIEWS",
    "FIGURE_S4_KAPPA_QC_PREP_INTERACTION"
  ),
  WidthInches = c(
    7.8, 7.8, 7.8, 7.8, 7.8, 7.8,
    7.8, 7.8, 7.8, 7.8, 7.8
  ),
  HeightInches = c(
    6.5, 4.25, 4.8, 7.2, 4.2, 4.2,
    7.0, 4.8, 7.0, 6.7, 4.4
  ),
  Formats = "PDF;SVG;PNG600dpi"
)
write_csv(figure_manifest, "FIGURE_MANIFEST.csv")

figure_captions <- data.table(
  Figure = figure_manifest$Figure,
  Title = c(
    "Update-rule scenarios",
    "Observed and adjusted effects of canalisation and inhibitory-edge probability",
    "Rewiring modifies the rule effects",
    "Network size and neighbourhood size modify period dynamics",
    "Basin summaries estimated from sampled initial states",
    "Recovery following a one-bit perturbation",
    "Perturbation and recovery summaries",
    "Rule effects at representative rewiring probabilities",
    "Network-size pairwise views from the primary model",
    "Neighbourhood-size pairwise views from the primary model",
    "Connectivity-dependent canalisation effect"
  ),
  Caption = c(
    "Examples of the update-rule scenarios. A node with z_j = 0 uses the signed-threshold rule. For a node with z_j = 1, an active dominant input fixes the target state according to the sign of the dominant edge; when the dominant input is inactive, the node uses the fallback signed-threshold rule.",
    "Observed and adjusted effects of canalisation and inhibitory-edge probability. Panel A shows the observed geometric mean period, panel B shows the observed fixed-point rate, and panel C shows the adjusted mean log period from the primary model, across q_c and p_rep. Each independently generated graph-rule replicate contributes equally to panels A and B. Panel C predictions are averaged equally over controlled-parameter combinations present in the simulations.",
    "Rewiring modifies the rule effects. Adjusted mean log period is shown across the six controlled rewiring probabilities for q_c = 0, 0.50, and 1.00. Panels show the representative inhibitory-edge probabilities p_rep = 0.10, 0.24, and 0.50. Predictions are averaged equally over the simulated N-kappa combinations.",
    "Network size and neighbourhood size modify period dynamics. Panel A shows adjusted mean log period across N and the actual Watts-Strogatz neighbourhood parameter kappa, restricted to simulated N-kappa combinations. Panels B and C show the interactions of kappa with q_c and p_rep, respectively, at N = 100, where every kappa value was simulated.",
    "Basin summaries estimated from sampled initial states in the fixed-system design. Panels show the number of attractors reached, the largest observed basin fraction, and observed basin entropy.",
    "Recovery following a one-bit perturbation in the fixed-system design. Panels show return probability within 2000 updates, restricted mean return time, and the Hamming fraction after 25 updates.",
    "Perturbation and recovery summaries. Panels show return probability, restricted mean return time, mean and median return time among successful returns, and Hamming fractions after 1 and 10 updates.",
    "Adjusted canalisation and inhibitory-edge-probability effects at representative low, intermediate, and high rewiring probabilities, p = 0.01, 0.20, and 0.60.",
    "Network-size pairwise views from the primary model. The N-kappa panel is restricted to observed support. The N by p, N by q_c, and N by p_rep panels are evaluated at kappa = 1, the common neighbourhood-size support across all N.",
    "Neighbourhood-size pairwise views from the primary model. The kappa by p, kappa by q_c, and kappa by p_rep panels are evaluated at N = 100, where every kappa value was simulated.",
    "Connectivity-dependent canalisation effect from the primary model. The figure shows the adjusted q_c = 1 versus q_c = 0 contrast across kappa and p_rep at N = 100. Negative values indicate shorter periods under full canalisation."
  )
)
write_csv(figure_captions, "FIGURE_TITLES_AND_CAPTIONS.csv")

writeLines(capture.output(sessionInfo()), file.path(AUDIT_DIR, "SESSION_INFO.txt"))
writeLines(capture.output({
  cat("Output directory:", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE), "\n")
  cat("Period input:", normalizePath(PERIOD_FILE, winslash = "/", mustWork = FALSE), "\n")
  cat("Basin input:", normalizePath(BASIN_FILE, winslash = "/", mustWork = FALSE), "\n")
  cat("Perturbation input:", normalizePath(PERTURBATION_FILE, winslash = "/", mustWork = FALSE), "\n")
  cat("Hamming input:", normalizePath(HAMMING_FILE, winslash = "/", mustWork = FALSE), "\n")
  cat("Analysis version:", ANALYSIS_VERSION, "\n")
  cat("Selected primary model:", primary_model_name, "\n")
  cat("Structural parameter used in the model: kappa = max(floor(d_input / 2), 1)\n")
  cat("Primary formula:", paste(deparse(formula(primary_model)), collapse = " "), "\n")
  cat("Graph rows:", nrow(period_graph), "\n")
  cat("Design cells:", nrow(design_cells), "\n")
  cat("Basin systems:", nrow(basin_system), "\n")
  cat("Perturbation systems:", nrow(perturbation_system), "\n")
  cat("Perturbation trials:", nrow(perturbation_trial), "\n")
  cat("CSV tables:", length(list.files(TABLE_DIR, pattern = "\\.csv$")), "\n")
  cat("LaTeX tables:", length(list.files(LATEX_DIR, pattern = "\\.tex$")), "\n")
  cat("Figure files:", length(list.files(FIGURE_DIR)), "\n")
}), file.path(AUDIT_DIR, "RUN_SUMMARY.txt"))

cat("Completed.\n")
cat("Output directory:", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE), "\n")
