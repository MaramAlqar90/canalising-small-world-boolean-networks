
suppressPackageStartupMessages({
  library(Rcpp)
  library(igraph)
  library(data.table)
})

source("simulate_canalising_SW_reduced_threshold.R")

run_period_design_safe <- function(workers = 1L, ...) {
  if (workers != 1L) {
    warning(
      "Historical parallel RNG/provenance is not verified. ",
      "Running the thesis-consistent reconstruction sequentially."
    )
  }
  run_paper2_period_design(...)
}

run_fixed_system_design_safe <- function(workers = 1L, ...) {
  if (workers != 1L) {
    warning(
      "Historical parallel RNG/provenance is not verified. ",
      "Running the thesis-consistent reconstruction sequentially."
    )
  }
  run_paper2_fixed_system_design(...)
}

execute_sw_graph_replicates_parallel <- function(...) {
  run_period_design_safe(...)
}

execute_sw_basin_and_perturbation_subset <- function(...) {
  run_fixed_system_design_safe(...)
}
