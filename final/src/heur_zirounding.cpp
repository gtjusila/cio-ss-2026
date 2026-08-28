#include "heur_zirounding.hpp"

#include <scip/pub_lp.h>
#include <scip/pub_var.h>
#include <scip/scip_branch.h>
#include <scip/scip_lp.h>
#include <scip/scip_numerics.h>
#include <scip/scip_solvingstats.h>
#include <scip/type_heur.h>
#include <scip/type_result.h>
#include <scip/type_retcode.h>
#include <scip/type_var.h>

#include <cmath>
#include <cstddef>
#include <memory>
#include <ranges>
#include <span>
#include <vector>

#include "utils.hpp"

static constexpr SCIP_Real min_shift =
    1e-4;  // Minimum change to the incumbent for a move to be considered
static constexpr SCIP_Real min_slack =
    1e-5;  // Minimum slack that a constraint should have before move can be
           // generated from it
static constexpr SCIP_Real max_iter = 1'000'000;

/* Do not change start */
struct SolDeleter {
  SCIP* scip;
  SolDeleter(SCIP* scip) : scip(scip) {};
  void operator()(SCIP_Sol* p) { SCIPfreeSol(scip, &p); };
};
using SolPtr = std::unique_ptr<SCIP_SOL, SolDeleter>;
bool add_solution_helper(SCIP* scip, SCIP_HEUR* heur,
                         std::span<SCIP_COL* const> lp_cols,
                         std::span<SCIP_Real> incumbent) {
  SolDeleter deleter(scip);
  SolPtr sol(nullptr, deleter);
  CALL_CHECK(SCIPcreateSol(scip, std::out_ptr(sol), heur));
  for (auto [col, val] : std::ranges::zip_view(lp_cols, incumbent)) {
    CALL_CHECK(SCIPsetSolVal(scip, sol.get(), SCIPcolGetVar(col), val));
  }
  unsigned int stored;
  CALL_CHECK(
      SCIPtrySol(scip, sol.get(), FALSE, TRUE, TRUE, TRUE, TRUE, &stored));
  return stored;
}
/* Do not change end*/

bool is_any_fractional_remaining(SCIP* scip, std::span<SCIP_COL* const> lp_cols,
                                 std::span<SCIP_Real> incumbent) {
  // HINT use SCIPcolIsIntegral to check if col is binary or integer
  // Use SCIPisFeasIntegral to check if value is within
  // integer tolerances
  for (auto [col, val] : std::ranges::zip_view(lp_cols, incumbent)) {
    if (SCIPcolIsIntegral(col) && !SCIPisFeasIntegral(scip, val)) return true;
  }
  return false;
}

SCIP_Real compute_incumbent_activity(SCIP_Row* row,
                                     std::span<SCIP_Real> incumbent) {
  // Compute incumbent activity for the current row
  // HINT create a coeff span using SCIProwGetVals, SCIProwGetNNonz and a cols
  // span using SCIProwGetCols, and SCIProwGetNNonz. Finally use SCIPcolGetLPPos
  // row activity = row constant + sum over nonzero coeff*value. The constant
  // term can be obtianed via SCIProwGetConstant
  SCIP_Real activity = SCIProwGetConstant(row);
  std::span<const SCIP_Real> vals(SCIProwGetVals(row), SCIProwGetNNonz(row));
  std::span<SCIP_COL* const> cols(SCIProwGetCols(row), SCIProwGetNNonz(row));
  for (auto [val, col] : std::ranges::zip_view(vals, cols)) {
    auto col_idx = SCIPcolGetLPPos(col);
    activity += val * incumbent[col_idx];
  }
  return activity;
}

SCIP_DECL_HEUREXEC(ZIRoundHeur::scip_exec) {
  /* Dont touch start from here */
  *result = SCIP_DIDNOTRUN;
  assert(SCIPhasCurrentNodeLP(scip));
  assert(SCIPgetNLPs(scip) != nlpsolve);

  // Do not call heuristic if node was already detected to be infeasible
  if (nodeinfeasible) return SCIP_OKAY;

  // Only call heuristic if an optimal LP solution is at hand
  if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL) return SCIP_OKAY;

  // Only call heuristic if LP objective value is smaller than the cutoff bound
  if (SCIPisGE(scip, SCIPgetLPObjval(scip), SCIPgetCutoffbound(scip)))
    return SCIP_OKAY;

  // If LP Solution satisfy integrality constraint then skip the heuristic
  if (SCIPgetNLPBranchCands(scip) == 0) return SCIP_OKAY;

  // Make sure we are at a new LP solution before heuristic is called
  nlpsolve = SCIPgetNLPs(scip);

  *result = SCIP_DIDNOTFIND;
  /* Don't touch end here */

  /* Start editing from here onwards*/
  // Get LP informations, construct spans from SCIPgetLPRows, SCIPgetNLPRows,
  // SCIPgetLPCols, SCIPgetNLPCols
  std::span<SCIP_ROW* const> lp_rows(SCIPgetLPRows(scip), SCIPgetNLPRows(scip));
  std::span<SCIP_COL* const> lp_cols(SCIPgetLPCols(scip), SCIPgetNLPCols(scip));

  // Read the incumbent LP solution from SCIP see SCIPcolGetPrimsol
  std::vector<SCIP_Real> incumbent(std::ssize(lp_cols));
  for (int i = 0; i < std::ssize(lp_cols); ++i) {
    incumbent[i] = SCIPcolGetPrimsol(lp_cols[i]);
  }

  // Loop while there is any fractional variable remaining in the incumbent
  // solution and iter have not reached max iter
  int iter_count = 0;
  while (is_any_fractional_remaining(scip, lp_cols, incumbent)) {
    // check if var is fractional see is_fractional_remaining hint
    bool found_shift = false;

    for (auto [idx, col] : std::views::enumerate(lp_cols)) {
      // compute how far up / down this variable can be rounded
      // start with var_ub and var_lb via SCIPcolGetUb/ SCIPcolGetLb
      auto col_incumbent = incumbent[idx];
      if (!SCIPcolIsIntegral(col) || SCIPisFeasIntegral(scip, col_incumbent))
        continue;
      auto var_ub = SCIPcolGetUb(col);
      auto var_lb = SCIPcolGetLb(col);
      // for each row containing this variable, recompute the row's current
      // activity and see how it limits this variable's value
      // Hint use SCIPcolGetVals and SCIPcolGetNNonz to get the coefficient of
      // the columns on rows that are non zero and use SCIPcolGetRows to get the
      // rows for which the variable coefficients is nonzero
      std::span<const SCIP_Real> vals(SCIPcolGetVals(col),
                                      SCIPcolGetNNonz(col));
      std::span<SCIP_Row* const> rows(SCIPcolGetRows(col),
                                      SCIPcolGetNNonz(col));

      // Loop over the rows update var_ub and var_lb
      // Hint 1: during the for loop over the rows check if row is in LP using
      // if (!SCIProwIsInLP(row)) continue;
      // Hint 2: Use SCIProwGetLhs and SCIProwGetRhs to get left and right and
      // side remember to check if they are + or - infity using
      // SCIPisInfinity, note it applies if l is - SCIPinfinity(scip.get())
      // then SCIPisInifinity(-l) is true
      // Hint 3: SCIP provide with MIN and MAX function
      // Determine update rule for var_ub and var_lb
      for (auto [val, row] : std::ranges::zip_view(vals, rows)) {
        if (!SCIProwIsInLP(row)) continue;
        if (val == 0.0) continue;

        auto activity = compute_incumbent_activity(row, incumbent);

        auto l = SCIProwGetLhs(row);
        auto r = SCIProwGetRhs(row);
        auto base = col_incumbent - activity / val;

        if (val > 0.0) {
          if (!SCIPisInfinity(scip, r)) var_ub = MIN(var_ub, base + r / val);
          if (!SCIPisInfinity(scip, -l)) var_lb = MAX(var_lb, base + l / val);
        } else {
          if (!SCIPisInfinity(scip, -l)) var_ub = MIN(var_ub, base + l / val);
          if (!SCIPisInfinity(scip, r)) var_lb = MAX(var_lb, base + r / val);
        }
      }

      // check if the next integer up or down is reachable without breaking
      // anything use SCIPfeasCeil/Floor to get next integer up or down here,
      // since we don't need to ceil if the value is solution is already near
      // enough to integral To do LE and GE with tolerances use SCIPisLE and
      // SCIPisGE To do LT and GT with tolerances use SCIPisLT and SCIPisGT If
      // var is up and down roundable pick objective improving direction If var
      // is only up or only down roundable pick roundable direction If var is
      // neither roundable determine whether staying, going to the upper bound,
      // going to the lower bound give the best fractionality reduction
      // Tie can be broken by selecting the objective improving direction

      auto up_val = SCIPfeasCeil(scip, incumbent[idx]);
      auto down_val = SCIPfeasFloor(scip, incumbent[idx]);
      if (SCIPisLE(scip, up_val, var_ub)) {
        // round up
        incumbent[idx] = up_val;
        found_shift = true;
      } else if (SCIPisGE(scip, down_val, var_lb)) {
        // check if var can be rounded down
        incumbent[idx] = down_val;
        found_shift = true;
      } else {
        // neither roundable: pick best fractionality of staying/up/down
        auto frac_dist = [&](SCIP_Real v) {
          auto frac = SCIPfeasFrac(scip, v);
          return MIN(frac, 1.0 - frac);
        };

        auto best = incumbent[idx];
        auto best_frac = frac_dist(best);

        auto up_frac = frac_dist(var_ub);
        if (SCIPisLT(scip, up_frac, best_frac - min_shift)) {
          best = var_ub;
          best_frac = up_frac;
        }

        auto down_frac = frac_dist(var_lb);
        if (SCIPisLT(scip, down_frac, best_frac - min_shift)) {
          best = var_lb;
          best_frac = down_frac;
        }

        if (SCIPisGT(scip, ABS(best - incumbent[idx]), min_shift)) {
          incumbent[idx] = best;
          found_shift = true;
        }
      }
    }

    if (!found_shift) break;
    if ((++iter_count) >= max_iter) break;
  }

  // Add Sol if any is found
  // Use
  // bool stored = add_solution_helper(scip, heur, lp_cols, incumbent);
  // if(stored) *result = SCIP_FOUNDSOL;
  if (!is_any_fractional_remaining(scip, lp_cols, incumbent)) {
    bool stored = add_solution_helper(scip, heur, lp_cols, incumbent);
    if (stored) {
      *result = SCIP_FOUNDSOL;
    }
  }

  return SCIP_OKAY;
}
