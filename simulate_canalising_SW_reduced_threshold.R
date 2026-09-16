

suppressPackageStartupMessages({
  library(Rcpp)
  library(igraph)
  library(data.table)
})

# -------------------------------------------------------------------------
# Load the C++ detector
# -------------------------------------------------------------------------

load_thesis_detector <- function(
    cpp_file = "detect_attractor_canalising_reduced_threshold.cpp",
    rebuild = FALSE) {
  if (!file.exists(cpp_file)) {
    stop("Cannot find C++ detector: ", normalizePath(cpp_file, mustWork = FALSE))
  }
  Rcpp::sourceCpp(cpp_file, rebuild = rebuild)
  invisible(TRUE)
}

if (!exists("detect_attractor_rcpp", mode = "function")) {
  default_cpp <- file.path(getwd(), "detect_attractor_canalising_reduced_threshold.cpp")
  if (file.exists(default_cpp)) {
    load_thesis_detector(default_cpp, rebuild = FALSE)
  }
}

# -------------------------------------------------------------------------
# General helpers
# -------------------------------------------------------------------------

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

safe_fwrite <- function(x, file) {
  ensure_dir(dirname(file))
  data.table::fwrite(x, file)
  invisible(file)
}

state_to_string <- function(x) paste(as.integer(x), collapse = "")

seed_from_index <- function(base_seed, index, stream_offset = 0) {

  modulus <- .Machine$integer.max - 1
  x <- (as.double(base_seed) + as.double(stream_offset) + as.double(index)) %% modulus
  as.integer(x + 1)
}

initial_state_bernoulli_half <- function(n) {
  as.integer(rbinom(n, size = 1L, prob = 0.5))
}

kappa_from_d <- function(d) {
  max(floor(as.integer(d) / 2), 1L)
}

all_d_values <- function(n) {
  seq_len(floor(n / 10))
}

representative_d_values_fixed <- function(n) {
  dmax <- floor(n / 10)
  unique(as.integer(c(1L, ceiling(dmax / 2), dmax)))
}

# -------------------------------------------------------------------------
#  Watts--Strogatz construction
# -------------------------------------------------------------------------

ring_lattice_edges <- function(n, kappa) {
  if (n < 3L) stop("n must be at least 3.")
  if (kappa < 1L || kappa > floor((n - 1L) / 2L)) {
    stop("kappa must satisfy 1 <= kappa <= floor((n-1)/2).")
  }

  out <- matrix(NA_integer_, nrow = n * kappa, ncol = 2L)
  z <- 0L
  for (u in seq_len(n)) {
    for (r in seq_len(kappa)) {
      z <- z + 1L
      v <- ((u - 1L + r) %% n) + 1L
      out[z, ] <- c(u, v)
    }
  }
  out
}

generate_ws_undirected_thesis <- function(n, d, p_rewire) {
  kappa <- kappa_from_d(d)
  original_edges <- ring_lattice_edges(n, kappa)

  U <- matrix(FALSE, nrow = n, ncol = n)
  for (r in seq_len(nrow(original_edges))) {
    u <- original_edges[r, 1L]
    v <- original_edges[r, 2L]
    U[u, v] <- TRUE
    U[v, u] <- TRUE
  }

  # Consider each original lattice edge once. The first endpoint u is retained.
  for (r in seq_len(nrow(original_edges))) {
    u <- original_edges[r, 1L]
    v <- original_edges[r, 2L]

    if (runif(1) < p_rewire) {
      current_neighbours <- which(U[u, ])
      candidates <- setdiff(seq_len(n), c(u, current_neighbours))

      if (length(candidates) > 0L) {
        w <- candidates[sample.int(length(candidates), 1L)]

        U[u, v] <- FALSE
        U[v, u] <- FALSE
        U[u, w] <- TRUE
        U[w, u] <- TRUE
      }
    }
  }

  if (any(diag(U))) stop("Self-loop created unexpectedly.")
  if (!all(U == t(U))) stop("Undirected adjacency lost symmetry.")

  expected_edges <- n * kappa
  actual_edges <- sum(U[upper.tri(U)])
  if (actual_edges != expected_edges) {
    stop("WS construction did not preserve N*kappa edges.")
  }

  g_und <- igraph::graph_from_adjacency_matrix(
    U * 1L, mode = "undirected", diag = FALSE
  )

  list(g_undirected = g_und, U = U, kappa = kappa)
}

orient_edges_once <- function(g_undirected) {
  ed <- igraph::as_edgelist(g_undirected, names = FALSE)

  if (nrow(ed) == 0L) {
    return(igraph::make_empty_graph(
      n = igraph::vcount(g_undirected), directed = TRUE
    ))
  }

  reverse_edge <- runif(nrow(ed)) < 0.5
  src <- ifelse(reverse_edge, ed[, 2L], ed[, 1L])
  dst <- ifelse(reverse_edge, ed[, 1L], ed[, 2L])

  igraph::graph_from_edgelist(
    cbind(as.integer(src), as.integer(dst)),
    directed = TRUE
  )
}

assign_edge_signs <- function(g_directed, p_rep) {
  n <- igraph::vcount(g_directed)
  A <- matrix(0, nrow = n, ncol = n)
  ed <- igraph::as_edgelist(g_directed, names = FALSE)

  if (nrow(ed) > 0L) {
    signs <- ifelse(runif(nrow(ed)) < p_rep, -1, 1)
    A[cbind(ed[, 1L], ed[, 2L])] <- signs
  }

  A
}

calculate_graph_metrics <- function(g_undirected, g_directed) {
  msp <- suppressWarnings(
    igraph::mean_distance(g_directed, directed = TRUE, unconnected = TRUE)
  )
  if (!is.finite(msp)) msp <- 0

  C <- suppressWarnings(
    igraph::transitivity(g_undirected, type = "global", isolates = "zero")
  )
  if (!is.finite(C)) C <- 0

  list(
    N = igraph::vcount(g_directed),
    m = igraph::ecount(g_directed),
    AverageTotalDegree = mean(igraph::degree(g_directed, mode = "all")),
    AverageInDegree = mean(igraph::degree(g_directed, mode = "in")),
    ClusteringCoefficient = C,
    MeanShortestPath = msp
  )
}

# -------------------------------------------------------------------------
# Depth-1 dominant-regulator assignment
# -------------------------------------------------------------------------

assign_depth1_rule <- function(A, q_c) {
  n <- ncol(A)

  z <- integer(n)
  dominant <- integer(n)
  alpha <- integer(n)
  beta <- integer(n)

  for (j in seq_len(n)) {
    incoming <- which(A[, j] != 0)

    if (length(incoming) > 0L && runif(1) < q_c) {
      # Correct for both one-element and multi-element incoming sets.
      dom <- incoming[sample.int(length(incoming), 1L)]

      if (A[dom, j] == 0) {
        stop("Internal error: dominant regulator is not an incoming neighbour.")
      }

      z[j] <- 1L
      dominant[j] <- as.integer(dom)
      alpha[j] <- 1L
      beta[j] <- if (A[dom, j] > 0) 1L else 0L
    }
  }

  list(
    canalizing_active = z,
    dominant_regulator = dominant,
    canalizing_input = alpha,
    canalized_output = beta
  )
}

validate_rule <- function(A, rule) {
  active <- which(rule$canalizing_active == 1L)
  if (length(active) == 0L) return(invisible(TRUE))

  dom <- rule$dominant_regulator[active]
  if (any(dom < 1L | dom > nrow(A))) {
    stop("Invalid dominant-regulator index.")
  }

  vals <- A[cbind(dom, active)]
  if (any(vals == 0)) {
    stop("At least one assigned dominant regulator is not an incoming neighbour.")
  }

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Threshold construction
# -------------------------------------------------------------------------

draw_threshold_vector <- function(A, rule, sigma = 1e-4) {
  validate_rule(A, rule)

  n <- ncol(A)
  eps <- rnorm(n, mean = 0, sd = sigma)
  sums <- colSums(A)
  theta <- 0.5 * sums + eps

  active <- which(rule$canalizing_active == 1L)
  if (length(active) > 0L) {
    dom <- rule$dominant_regulator[active]
    theta[active] <-
      0.5 * (sums[active] - A[cbind(dom, active)]) + eps[active]
  }

  as.numeric(theta)
}

# -------------------------------------------------------------------------
# System constructors
# -------------------------------------------------------------------------

build_graph_rule_replicate <- function(n, d, p_rewire, p_rep, q_c) {
  ws <- generate_ws_undirected_thesis(n, d, p_rewire)
  g_dir <- orient_edges_once(ws$g_undirected)
  A <- assign_edge_signs(g_dir, p_rep)
  rule <- assign_depth1_rule(A, q_c)
  validate_rule(A, rule)

  list(
    n = n,
    d = d,
    kappa = ws$kappa,
    p = p_rewire,
    p_rep = p_rep,
    q_c = q_c,
    g_undirected = ws$g_undirected,
    g_directed = g_dir,
    A = A,
    rule = rule,
    metrics = calculate_graph_metrics(ws$g_undirected, g_dir)
  )
}

build_fixed_system <- function(
    n, d, p_rewire, p_rep, q_c, sigma = 1e-4) {

  obj <- build_graph_rule_replicate(n, d, p_rewire, p_rep, q_c)

  # One threshold vector is retained for the entire fixed Boolean system.
  obj$thresholds <- draw_threshold_vector(obj$A, obj$rule, sigma)
  obj$sigma <- sigma
  obj
}

# -------------------------------------------------------------------------
# Detector wrappers
# -------------------------------------------------------------------------

detect_from_system <- function(
    system, initial_state, thresholds,
    max_iters = 1000000L, return_states = FALSE) {

  if (!exists("detect_attractor_rcpp", mode = "function")) {
    stop("C++ detector is not loaded. Run load_thesis_detector().")
  }

  detect_attractor_rcpp(
    A = system$A,
    initial_state = as.integer(initial_state),
    thresholds = as.numeric(thresholds),
    canalizing_active = as.integer(system$rule$canalizing_active),
    dominant_regulator = as.integer(system$rule$dominant_regulator),
    canalizing_input = as.integer(system$rule$canalizing_input),
    canalized_output = as.integer(system$rule$canalized_output),
    max_iters = as.integer(max_iters),
    return_states = isTRUE(return_states)
  )
}

one_step_from_system <- function(system, state, thresholds) {
  one_step_update_rcpp(
    A = system$A,
    state = as.integer(state),
    thresholds = as.numeric(thresholds),
    canalizing_active = as.integer(system$rule$canalizing_active),
    dominant_regulator = as.integer(system$rule$dominant_regulator),
    canalizing_input = as.integer(system$rule$canalizing_input),
    canalized_output = as.integer(system$rule$canalized_output)
  )
}

# -------------------------------------------------------------------------
# Period design
# -------------------------------------------------------------------------

summarise_period_trajectories <- function(periods, attempted) {
  good <- is.finite(periods)

  if (!any(good)) {
    return(data.table(
      AttemptedTrajectories = attempted,
      DetectedTrajectories = 0L,
      DetectionRate = 0,
      MeanPeriod = NA_real_,
      MeanLogPeriod = NA_real_,
      GeometricMeanPeriod = NA_real_,
      FixedPointRate = NA_real_
    ))
  }

  L <- periods[good]
  Y <- mean(log(L))

  data.table(
    AttemptedTrajectories = attempted,
    DetectedTrajectories = length(L),
    DetectionRate = length(L) / attempted,
    MeanPeriod = mean(L),
    MeanLogPeriod = Y,
    GeometricMeanPeriod = exp(Y),
    FixedPointRate = mean(L == 1)
  )
}

period_replicate_row <- function(
    system,
    n_trajectories = 25L,
    sigma = 1e-4,
    max_iters = 1000000L,
    base_seed = 15092026L,
    graph_global_id,
    graph_replicate = NA_integer_) {

  periods <- rep(NA_real_, n_trajectories)

  for (r in seq_len(n_trajectories)) {
    trajectory_global_id <-
      (as.double(graph_global_id) - 1) * n_trajectories + r

    set.seed(seed_from_index(
      base_seed,
      trajectory_global_id,
      stream_offset = 50000000
    ))

    S0 <- initial_state_bernoulli_half(system$n)

    # Thesis requirement: redraw thresholds for every period trajectory.
    thresholds <- draw_threshold_vector(
      system$A, system$rule, sigma = sigma
    )

    out <- detect_from_system(
      system,
      initial_state = S0,
      thresholds = thresholds,
      max_iters = max_iters,
      return_states = FALSE
    )

    if (isTRUE(out$found)) {
      periods[r] <- as.numeric(out$attractor_period)
    }
  }

  m <- system$metrics

  cbind(
    data.table(
      N = system$n,
      d = system$d,
      kappa = system$kappa,
      p = system$p,
      p_rep = system$p_rep,
      q_c = system$q_c,
      GraphReplicate = graph_replicate,
      NumberOfEdges = m$m,
      AverageTotalDegree = m$AverageTotalDegree,
      AverageInDegree = m$AverageInDegree,
      ClusteringCoefficient = m$ClusteringCoefficient,
      MeanShortestPath = m$MeanShortestPath,
      CanalisingAssigned = sum(system$rule$canalizing_active == 1L)
    ),
    summarise_period_trajectories(periods, n_trajectories)
  )
}

paper2_period_grid <- function(
    n_values = seq(10, 100, by = 10),
    p_values = c(0.01, 0.05, 0.10, 0.20, 0.40, 0.60),
    p_rep_values = c(0.10, 0.24, 0.30, 0.41, 0.50),
    q_c_values = c(0, 0.25, 0.50, 0.75, 1)) {

  rows <- list()
  z <- 0L

  for (n in n_values) {
    for (d in all_d_values(n)) {
      for (p in p_values) {
        for (p_rep in p_rep_values) {
          for (q_c in q_c_values) {
            z <- z + 1L
            rows[[z]] <- data.table(
              N = n, d = d, p = p, p_rep = p_rep, q_c = q_c
            )
          }
        }
      }
    }
  }

  rbindlist(rows)
}

run_paper2_period_design <- function(
    num_graph_replicates = 100L,
    num_initial_states_per_graph = 25L,
    sigma = 1e-4,
    max_iters = 1000000L,
    output_dir = "paper2_period_reconstructed",
    base_seed = 15092026L,
    overwrite = FALSE) {

  ensure_dir(output_dir)
  grid <- paper2_period_grid()

  if (nrow(grid) != 8250L) {
    stop("Period grid should contain 8,250 original cells; found ", nrow(grid))
  }

  expected_rows <- nrow(grid) * num_graph_replicates
  message("Expected graph-rule replicates: ", format(expected_rows, big.mark = ","))

  for (cell in seq_len(nrow(grid))) {
    g <- grid[cell]

    out_file <- file.path(
      output_dir,
      sprintf(
        "period_N%03d_d%02d_p%.2f_prep%.2f_qc%.2f.csv",
        g$N, g$d, g$p, g$p_rep, g$q_c
      )
    )

    if (file.exists(out_file) && !overwrite) next

    cell_rows <- vector("list", num_graph_replicates)

    for (rep_id in seq_len(num_graph_replicates)) {
      graph_global_id <-
        (as.double(cell) - 1) * num_graph_replicates + rep_id

      set.seed(seed_from_index(
        base_seed,
        graph_global_id,
        stream_offset = 0
      ))

      system <- build_graph_rule_replicate(
        n = g$N,
        d = g$d,
        p_rewire = g$p,
        p_rep = g$p_rep,
        q_c = g$q_c
      )

      cell_rows[[rep_id]] <- period_replicate_row(
        system,
        n_trajectories = num_initial_states_per_graph,
        sigma = sigma,
        max_iters = max_iters,
        base_seed = base_seed,
        graph_global_id = graph_global_id,
        graph_replicate = rep_id
      )
    }

    safe_fwrite(rbindlist(cell_rows), out_file)

    if (cell %% 50L == 0L || cell == nrow(grid)) {
      message("Completed period cells: ", cell, " / ", nrow(grid))
    }
  }

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Fixed-system basin summaries
# -------------------------------------------------------------------------

basin_summary_for_system <- function(
    system,
    n_initial_states = 200L,
    max_iters = 1000000L,
    base_seed = 15092026L,
    system_index,
    system_id = NA_character_,
    return_run_table = FALSE) {

  ids <- rep(NA_character_, n_initial_states)
  periods <- rep(NA_real_, n_initial_states)

  for (r in seq_len(n_initial_states)) {
    basin_run_global_id <-
      (as.double(system_index) - 1) * n_initial_states + r

    set.seed(seed_from_index(
      base_seed,
      basin_run_global_id,
      stream_offset = 150000000
    ))

    S0 <- initial_state_bernoulli_half(system$n)

    out <- detect_from_system(
      system,
      initial_state = S0,
      thresholds = system$thresholds,
      max_iters = max_iters,
      return_states = FALSE
    )

    if (isTRUE(out$found)) {
      ids[r] <- as.character(out$attractor_id)
      periods[r] <- as.numeric(out$attractor_period)
    }
  }

  found <- !is.na(ids)
  n_found <- sum(found)

  if (n_found == 0L) {
    summary <- data.table(
      BasinStarts = n_initial_states,
      BasinDetected = 0L,
      BasinDetectionRate = 0,
      MeanLogPeriod = NA_real_,
      GeometricMeanPeriod = NA_real_,
      FixedPointRate = NA_real_,
      K_obs = 0L,
      LargestObservedBasinFraction = NA_real_,
      ObservedBasinEntropy = NA_real_
    )
  } else {
    counts <- table(ids[found])
    pi_hat <- as.numeric(counts) / sum(counts)
    Y <- mean(log(periods[found]))

    summary <- data.table(
      BasinStarts = n_initial_states,
      BasinDetected = n_found,
      BasinDetectionRate = n_found / n_initial_states,
      MeanLogPeriod = Y,
      GeometricMeanPeriod = exp(Y),
      FixedPointRate = mean(periods[found] == 1),
      K_obs = length(counts),
      LargestObservedBasinFraction = max(pi_hat),
      ObservedBasinEntropy = -sum(pi_hat * log(pi_hat))
    )
  }

  if (!return_run_table) return(summary)

  list(
    summary = summary,
    runs = data.table(
      SystemID = system_id,
      BasinRun = seq_len(n_initial_states),
      Found = found,
      AttractorID = ids,
      Period = periods
    )
  )
}

# -------------------------------------------------------------------------
# One-bit perturbation trials
# -------------------------------------------------------------------------

hamming_fraction <- function(a, b) {
  mean(as.integer(a) != as.integer(b))
}

make_state_hash_environment <- function(strings) {
  e <- new.env(hash = TRUE, parent = emptyenv())
  for (s in strings) assign(s, TRUE, envir = e)
  e
}

run_one_bit_trial_from_detected_attractor <- function(
    system,
    attractor_detection,
    return_cap = 2000L) {

  if (!isTRUE(attractor_detection$found)) {
    stop("A detected attractor is required.")
  }

  attr_strings <- as.character(attractor_detection$attractor_state_strings)
  if (length(attr_strings) == 0L) {
    stop("Attractor states were not returned by the detector.")
  }

  # Thesis convention: first stored state of the detected attractor cycle.
  base_state <- as.integer(attractor_detection$first_attractor_state)

  # Flip exactly one uniformly selected node.
  flip_node <- sample.int(system$n, 1L)
  pert_state <- base_state
  pert_state[flip_node] <- 1L - pert_state[flip_node]

  unpert_state <- base_state
  attr_set <- make_state_hash_environment(attr_strings)

  return_time <- NA_integer_
  h1 <- h10 <- h25 <- NA_real_

  for (t in seq_len(return_cap)) {
    pert_state <- one_step_from_system(
      system, pert_state, system$thresholds
    )

    if (t <= 25L) {
      unpert_state <- one_step_from_system(
        system, unpert_state, system$thresholds
      )

      h <- hamming_fraction(unpert_state, pert_state)
      if (t == 1L) h1 <- h
      if (t == 10L) h10 <- h
      if (t == 25L) h25 <- h
    }

    if (is.na(return_time)) {
      key <- state_to_string(pert_state)
      if (exists(key, envir = attr_set, inherits = FALSE)) {
        return_time <- t
      }
    }

    if (t >= 25L && !is.na(return_time)) break
  }

  returned <- !is.na(return_time)

  data.table(
    FlipNode = flip_node,
    Returned = returned,
    ReturnTime = if (returned) return_time else NA_integer_,
    RestrictedReturnTime = if (returned) return_time else return_cap,
    Hamming1 = h1,
    Hamming10 = h10,
    Hamming25 = h25
  )
}

perturbation_summary_for_system <- function(
    system,
    n_trials = 30L,
    max_iters = 1000000L,
    return_cap = 2000L,
    max_base_attempts = 600L,
    base_seed = 15092026L,
    system_index,
    system_id = NA_character_) {

  trial_rows <- vector("list", n_trials)
  valid <- 0L
  attempts <- 0L

  while (valid < n_trials && attempts < max_base_attempts) {
    attempts <- attempts + 1L

    perturb_attempt_global_id <-
      (as.double(system_index) - 1) * max_base_attempts + attempts

    set.seed(seed_from_index(
      base_seed,
      perturb_attempt_global_id,
      stream_offset = 200000000
    ))

    S0 <- initial_state_bernoulli_half(system$n)

    out <- detect_from_system(
      system,
      initial_state = S0,
      thresholds = system$thresholds,
      max_iters = max_iters,
      return_states = TRUE
    )

    if (!isTRUE(out$found)) next

    valid <- valid + 1L

    trial <- run_one_bit_trial_from_detected_attractor(
      system,
      attractor_detection = out,
      return_cap = return_cap
    )

    trial[, `:=`(
      SystemID = system_id,
      Trial = valid,
      BaseSearchAttempt = attempts,
      BaseAttractorID = as.character(out$attractor_id),
      BaseAttractorPeriod = as.integer(out$attractor_period)
    )]

    setcolorder(
      trial,
      c(
        "SystemID", "Trial", "BaseSearchAttempt",
        "BaseAttractorID", "BaseAttractorPeriod",
        "FlipNode", "Returned", "ReturnTime",
        "RestrictedReturnTime", "Hamming1",
        "Hamming10", "Hamming25"
      )
    )

    trial_rows[[valid]] <- trial
  }

  if (valid < n_trials) {
    stop(
      "Only ", valid, " valid perturbation trials were obtained after ",
      attempts, " base-attractor attempts for system ", system_id, "."
    )
  }

  trials <- rbindlist(trial_rows)
  successful_times <- trials[Returned == TRUE, ReturnTime]

  summary <- data.table(
    PerturbationTrials = nrow(trials),
    BaseSearchAttempts = attempts,
    ReturnProbability = mean(trials$Returned),
    RestrictedMeanReturnTime = mean(trials$RestrictedReturnTime),
    SuccessfulMeanReturnTime =
      if (length(successful_times)) mean(successful_times) else NA_real_,
    SuccessfulMedianReturnTime =
      if (length(successful_times)) median(successful_times) else NA_real_,
    MeanHamming1 = mean(trials$Hamming1),
    MeanHamming10 = mean(trials$Hamming10),
    MeanHamming25 = mean(trials$Hamming25)
  )

  list(summary = summary, trials = trials)
}

# -------------------------------------------------------------------------
# Fixed-system design: exactly 7,200 systems
# -------------------------------------------------------------------------

paper2_fixed_grid <- function(
    n_values = c(20L, 60L, 100L),
    p_values = c(0.01, 0.10, 0.60),
    p_rep_values = c(0.30, 0.50),
    q_c_values = c(0, 0.25, 0.50, 0.75, 1)) {

  rows <- list()
  z <- 0L

  for (n in n_values) {
    for (d in representative_d_values_fixed(n)) {
      for (p in p_values) {
        for (p_rep in p_rep_values) {
          for (q_c in q_c_values) {
            z <- z + 1L
            rows[[z]] <- data.table(
              N = n,
              d = d,
              kappa = kappa_from_d(d),
              p = p,
              p_rep = p_rep,
              q_c = q_c
            )
          }
        }
      }
    }
  }

  rbindlist(rows)
}

run_paper2_fixed_system_design <- function(
    systems_per_parameter_combination = 30L,
    basin_initial_states = 200L,
    perturbation_trials = 30L,
    sigma = 1e-4,
    max_iters = 1000000L,
    return_cap = 2000L,
    max_base_attempts = 600L,
    output_dir = "paper2_fixed_system_reconstructed",
    base_seed = 15092026L,
    overwrite = FALSE,
    save_basin_runs = FALSE,
    save_perturbation_trials = TRUE) {

  ensure_dir(output_dir)

  grid <- paper2_fixed_grid()
  if (nrow(grid) != 240L) {
    stop("Fixed-system original grid should have 240 combinations; found ",
         nrow(grid))
  }

  expected_systems <- nrow(grid) * systems_per_parameter_combination
  if (expected_systems != 7200L) {
    stop("This thesis design should generate exactly 7,200 fixed systems.")
  }

  message("Expected fixed Boolean systems: 7,200")

  for (cell in seq_len(nrow(grid))) {
    g <- grid[cell]

    stem <- sprintf(
      "fixed_N%03d_d%02d_p%.2f_prep%.2f_qc%.2f",
      g$N, g$d, g$p, g$p_rep, g$q_c
    )

    system_file <- file.path(output_dir, paste0(stem, "_systems.csv"))
    trial_file <- file.path(output_dir, paste0(stem, "_perturbation_trials.csv"))
    basin_file <- file.path(output_dir, paste0(stem, "_basin_runs.csv"))

    if (file.exists(system_file) && !overwrite) next

    system_rows <- vector("list", systems_per_parameter_combination)
    trial_rows <- if (save_perturbation_trials) vector("list", systems_per_parameter_combination) else NULL
    basin_rows <- if (save_basin_runs) vector("list", systems_per_parameter_combination) else NULL

    for (sys_rep in seq_len(systems_per_parameter_combination)) {
      system_global_id <-
        (as.double(cell) - 1) * systems_per_parameter_combination + sys_rep

      system_id <- sprintf(
        "N%d_d%d_p%.2f_prep%.2f_qc%.2f_sys%02d",
        g$N, g$d, g$p, g$p_rep, g$q_c, sys_rep
      )

      set.seed(seed_from_index(
        base_seed,
        system_global_id,
        stream_offset = 100000000
      ))

      system <- build_fixed_system(
        n = g$N,
        d = g$d,
        p_rewire = g$p,
        p_rep = g$p_rep,
        q_c = g$q_c,
        sigma = sigma
      )

      basin <- basin_summary_for_system(
        system,
        n_initial_states = basin_initial_states,
        max_iters = max_iters,
        base_seed = base_seed,
        system_index = system_global_id,
        system_id = system_id,
        return_run_table = save_basin_runs
      )

      basin_summary <- if (save_basin_runs) basin$summary else basin

      pert <- perturbation_summary_for_system(
        system,
        n_trials = perturbation_trials,
        max_iters = max_iters,
        return_cap = return_cap,
        max_base_attempts = max_base_attempts,
        base_seed = base_seed,
        system_index = system_global_id,
        system_id = system_id
      )

      m <- system$metrics

      system_rows[[sys_rep]] <- cbind(
        data.table(
          SystemID = system_id,
          N = g$N,
          d = g$d,
          kappa = g$kappa,
          p = g$p,
          p_rep = g$p_rep,
          q_c = g$q_c,
          SystemReplicate = sys_rep,
          NumberOfEdges = m$m,
          AverageTotalDegree = m$AverageTotalDegree,
          AverageInDegree = m$AverageInDegree,
          ClusteringCoefficient = m$ClusteringCoefficient,
          MeanShortestPath = m$MeanShortestPath,
          CanalisingAssigned = sum(system$rule$canalizing_active == 1L)
        ),
        basin_summary,
        pert$summary
      )

      if (save_perturbation_trials) {
        trial_rows[[sys_rep]] <- pert$trials
      }

      if (save_basin_runs) {
        basin_rows[[sys_rep]] <- basin$runs
      }
    }

    safe_fwrite(rbindlist(system_rows, fill = TRUE), system_file)

    if (save_perturbation_trials) {
      safe_fwrite(rbindlist(trial_rows, fill = TRUE), trial_file)
    }

    if (save_basin_runs) {
      safe_fwrite(rbindlist(basin_rows, fill = TRUE), basin_file)
    }

    if (cell %% 20L == 0L || cell == nrow(grid)) {
      message("Completed fixed-system cells: ", cell, " / ", nrow(grid))
    }
  }

  invisible(TRUE)
}

# -------------------------------------------------------------------------
# Design checks
# -------------------------------------------------------------------------

verify_thesis_design_counts <- function() {
  pg <- paper2_period_grid()
  fg <- paper2_fixed_grid()

  out <- data.table(
    Quantity = c(
      "Period original cells",
      "Period graph-rule replicates",
      "Period attempted trajectories",
      "Fixed-system original combinations",
      "Fixed Boolean systems",
      "Fixed-system basin runs",
      "Perturbation trials"
    ),
    Expected = c(
      nrow(pg),
      nrow(pg) * 100L,
      nrow(pg) * 100L * 25L,
      nrow(fg),
      nrow(fg) * 30L,
      nrow(fg) * 30L * 200L,
      nrow(fg) * 30L * 30L
    ),
    ThesisTarget = c(
      8250L,
      825000L,
      20625000L,
      240L,
      7200L,
      1440000L,
      216000L
    )
  )

  out[, Matches := Expected == ThesisTarget]
  out
}

run_small_smoke_test <- function(base_seed = 12345L) {
  set.seed(base_seed)

  gr <- build_graph_rule_replicate(
    n = 10, d = 1, p_rewire = 0.10,
    p_rep = 0.30, q_c = 0.50
  )

  pr <- period_replicate_row(
    gr,
    n_trajectories = 3,
    sigma = 1e-4,
    max_iters = 10000,
    base_seed = base_seed,
    graph_global_id = 1,
    graph_replicate = 1
  )

  fs <- build_fixed_system(
    n = 20, d = 1, p_rewire = 0.10,
    p_rep = 0.30, q_c = 0.50,
    sigma = 1e-4
  )

  bs <- basin_summary_for_system(
    fs,
    n_initial_states = 5,
    max_iters = 10000,
    base_seed = base_seed,
    system_index = 1,
    system_id = "smoke"
  )

  pt <- perturbation_summary_for_system(
    fs,
    n_trials = 2,
    max_iters = 10000,
    return_cap = 100,
    max_base_attempts = 20,
    base_seed = base_seed,
    system_index = 1,
    system_id = "smoke"
  )

  list(
    design_counts = verify_thesis_design_counts(),
    period_test = pr,
    basin_test = bs,
    perturbation_test = pt$summary,
    perturbation_trials = pt$trials
  )
}
