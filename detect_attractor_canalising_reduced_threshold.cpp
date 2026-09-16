
#include <Rcpp.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <algorithm>

using namespace Rcpp;

static std::string state_vec_to_string_cpp(const IntegerVector& x) {
  std::string s;
  s.reserve(x.size());
  for (int i = 0; i < x.size(); ++i) {
    s.push_back(x[i] == 0 ? '0' : '1');
  }
  return s;
}

// [[Rcpp::export]]
IntegerVector one_step_update_rcpp(
    NumericMatrix A,
    IntegerVector state,
    NumericVector thresholds,
    IntegerVector canalizing_active,
    IntegerVector dominant_regulator,
    IntegerVector canalizing_input,
    IntegerVector canalized_output) {

  const int n = A.ncol();

  if (A.nrow() != n) stop("A must be square.");
  if (state.size() != n) stop("state length mismatch.");
  if (thresholds.size() != n) stop("thresholds length mismatch.");
  if (canalizing_active.size() != n) stop("canalizing_active length mismatch.");
  if (dominant_regulator.size() != n) stop("dominant_regulator length mismatch.");
  if (canalizing_input.size() != n) stop("canalizing_input length mismatch.");
  if (canalized_output.size() != n) stop("canalized_output length mismatch.");

  IntegerVector next_state(n);

  for (int j = 0; j < n; ++j) {
    bool forced = false;

    if (canalizing_active[j] == 1) {
      // R stores the regulator index in 1-based form.
      const int dom1 = dominant_regulator[j];
      if (dom1 < 1 || dom1 > n) {
        stop("Assigned dominant regulator is outside 1..N.");
      }
      const int dom = dom1 - 1;

      // Enforce the thesis rule v_{r(j)} in N^-(v_j).
      if (A(dom, j) == 0.0) {
        stop("Assigned dominant regulator is not an incoming neighbour.");
      }

      if (state[dom] == canalizing_input[j]) {
        next_state[j] = canalized_output[j];
        forced = true;
      }
    }

    if (!forced) {
      double input = 0.0;
      for (int i = 0; i < n; ++i) {
        input += A(i, j) * state[i];
      }
      next_state[j] = (input > thresholds[j]) ? 1 : 0;
    }
  }

  return next_state;
}


// [[Rcpp::export]]
List detect_attractor_rcpp(
    NumericMatrix A,
    IntegerVector initial_state,
    NumericVector thresholds,
    IntegerVector canalizing_active,
    IntegerVector dominant_regulator,
    IntegerVector canalizing_input,
    IntegerVector canalized_output,
    int max_iters = 1000000,
    bool return_states = false) {

  const int n = A.ncol();

  if (A.nrow() != n) stop("A must be square.");
  if (initial_state.size() != n) stop("initial_state length mismatch.");
  if (thresholds.size() != n) stop("thresholds length mismatch.");
  if (max_iters < 1) stop("max_iters must be positive.");

  std::unordered_map<std::string, int> first_seen;
  first_seen.reserve(4096);

  IntegerVector current = clone(initial_state);
  first_seen[state_vec_to_string_cpp(current)] = 0;

  for (int t = 1; t <= max_iters; ++t) {
    IntegerVector next_state = one_step_update_rcpp(
      A, current, thresholds,
      canalizing_active, dominant_regulator,
      canalizing_input, canalized_output
    );

    const std::string key = state_vec_to_string_cpp(next_state);
    auto it = first_seen.find(key);

    if (it != first_seen.end()) {
      const int t1 = it->second;
      const int period = t - t1;

      // The repeated state is the state first stored at the beginning
      // of the detected cycle. Reconstruct one full cycle from it.
      IntegerVector cycle_state = clone(next_state);
      std::vector<std::string> cycle_strings;
      cycle_strings.reserve(period);

      IntegerMatrix cycle_states;
      if (return_states) {
        cycle_states = IntegerMatrix(period, n);
      } else {
        cycle_states = IntegerMatrix(0, n);
      }

      std::string min_state;
      IntegerVector first_attractor_state = clone(cycle_state);

      for (int r = 0; r < period; ++r) {
        const std::string cs = state_vec_to_string_cpp(cycle_state);
        cycle_strings.push_back(cs);

        if (r == 0 || cs < min_state) min_state = cs;

        if (return_states) {
          for (int i = 0; i < n; ++i) {
            cycle_states(r, i) = cycle_state[i];
          }
        }

        cycle_state = one_step_update_rcpp(
          A, cycle_state, thresholds,
          canalizing_active, dominant_regulator,
          canalizing_input, canalized_output
        );
      }

      // In a deterministic map, two distinct attractor cycles cannot share a state.
      // The lexicographically smallest state plus the period is therefore a
      // phase-independent identifier for the detected attractor.
      const std::string canonical_id =
        std::string("A_") + min_state + std::string("_L") + std::to_string(period);

      CharacterVector cycle_string_out;
      if (return_states) {
        cycle_string_out = CharacterVector(period);
        for (int r = 0; r < period; ++r) {
          cycle_string_out[r] = cycle_strings[r];
        }
      } else {
        cycle_string_out = CharacterVector(0);
      }

      return List::create(
        Named("found") = true,
        Named("iteration_1") = t1,
        Named("iteration_2") = t,
        Named("transient_length") = t1,
        Named("attractor_period") = period,
        Named("attractor_id") = canonical_id,
        Named("first_attractor_state") = first_attractor_state,
        Named("attractor_states") = cycle_states,
        Named("attractor_state_strings") = cycle_string_out
      );
    }

    first_seen[key] = t;
    current = clone(next_state);
  }

  return List::create(
    Named("found") = false,
    Named("iteration_1") = R_NilValue,
    Named("iteration_2") = R_NilValue,
    Named("transient_length") = R_NilValue,
    Named("attractor_period") = R_NilValue,
    Named("attractor_id") = NA_STRING,
    Named("first_attractor_state") = IntegerVector(0),
    Named("attractor_states") = IntegerMatrix(0, n),
    Named("attractor_state_strings") = CharacterVector(0)
  );
}
