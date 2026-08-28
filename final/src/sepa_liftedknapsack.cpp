#include "sepa_liftedknapsack.hpp"

#include <scip/def.h>
#include <scip/pub_lp.h>
#include <scip/scip_numerics.h>
#include <scip/type_result.h>
#include <scip/type_var.h>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <ranges>
#include <span>
#include <utility>
#include <vector>

#include "utils.hpp"
/* Do not change start */
struct RowDeleter {
  SCIP* scip;
  RowDeleter(SCIP* scip) : scip(scip) {};
  void operator()(SCIP_ROW* p) { SCIPreleaseRow(scip, &p); }
};
using RowPtr = std::unique_ptr<SCIP_ROW, RowDeleter>;

void LiftedKnapsackSepa::add_cut(SCIP* scip, SCIP_SEPA* sepa,
                                 std::span<int> indices,
                                 std::span<SCIP_Real> coefficients,
                                 SCIP_Real lhs, SCIP_Real rhs, bool forcecut,
                                 SCIP_Result* result) {
  RowDeleter deleter(scip);
  RowPtr new_row(nullptr, deleter);
  CALL_CHECK(SCIPcreateEmptyRowSepa(
      scip, std::out_ptr(new_row), sepa,
      std::format("lifted_knapsack_{}", cut_count++).c_str(), lhs, rhs, false,
      false, false));
  CALL_CHECK(SCIPcacheRowExtensions(scip, new_row.get()));
  std::span<SCIP_COL*> cols(SCIPgetLPCols(scip), SCIPgetNLPCols(scip));
  for (auto [ind, coef] : std::views::zip(indices, coefficients)) {
    if (SCIPisZero(scip, coef)) continue;
    CALL_CHECK(
        SCIPaddVarToRow(scip, new_row.get(), SCIPcolGetVar(cols[ind]), coef));
  }
  CALL_CHECK(SCIPflushRowExtensions(scip, new_row.get()));

  if (SCIPisCutEfficacious(scip, nullptr, new_row.get())) {
    SCIP_Bool cutoff;
    CALL_CHECK(SCIPaddRow(scip, new_row.get(), forcecut, &cutoff));
    // CALL_CHECK(SCIPaddPoolCut(scip, new_row.get()));
    *result = SCIP_SEPARATED;
  }
}
/*Do not change end*/

SCIP_DECL_SEPAEXECLP(LiftedKnapsackSepa::scip_execlp) {
  /* Do not touch start*/
  *result = SCIP_DIDNOTRUN;

  // Check if LP is solved or fail because of some other condition (e.g.
  // numerics)
  if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL) return SCIP_OKAY;
  // Check if LP posess a basis, this might not be true if LP comes from an
  // interior point method
  if (!SCIPisLPSolBasic(scip)) return SCIP_OKAY;

  /*We are ready to proceed*/
  *result = SCIP_DIDNOTFIND;

  // Get LProws pointer
  std::span<SCIP_Row* const> lp_rows(SCIPgetLPRows(scip), SCIPgetNLPRows(scip));

  // We store pointers to rows corresponding to knapsack constraints
  if (!knapsack_constraints) {
    std::vector<SCIP_Row*> temp_knapsack_rows;
    for (auto row : lp_rows) {
      // We check if the row is an original row (i.e. it is part of the original
      // problem)
      if (SCIProwIsModifiable(row) ||
          SCIProwGetOrigintype(row) != SCIP_ROWORIGINTYPE_CONS ||
          SCIProwIsLocal(row))
        continue;

      // Next we check if the row is a knapsack constraint
      // A knapsack constraint is a linear constraint of the form:
      // sum_{i in I} a_i * x_i <= b
      // where a_i > 0 for all i in I, a_i in Z, b >= 0, x_i are binary
      // Although in theory rows may be the other way around i.e. >= with
      // negative b and negative coefficient we ignore the case for now
      auto constant = SCIProwGetConstant(row);
      auto lhs = SCIProwGetLhs(row);
      auto rhs = SCIProwGetRhs(row);
      if (!SCIPisInfinity(scip, -lhs) && !SCIPisInfinity(scip, rhs))
        continue;  // Row is two-sided
      if (!SCIPisInfinity(scip, -lhs)) {
        continue;  // Row is a <= constraint
      }
      if (!SCIPisGT(scip, rhs - constant, 0)) {
        continue;  // RHsi is not positive
      }
      // b must also be integral
      if (!SCIPisIntegral(scip, rhs - constant)) continue;
      // Check if all nonzero coefficient corresponds to binary variable
      std::span<SCIP_Col*> nonzero_col(SCIProwGetCols(row),
                                       SCIProwGetNNonz(row));
      bool all_binary = true;
      for (auto col : nonzero_col) {
        auto var = SCIPcolGetVar(col);
        if (!(SCIPvarGetType(var) == SCIP_VARTYPE_BINARY)) all_binary = false;
      }

      if (!all_binary) continue;

      // Nonzero coefficient must all be integer and positive
      std::span<SCIP_Real> vals(SCIProwGetVals(row), SCIProwGetNNonz(row));
      if (std::ranges::any_of(vals, [scip](SCIP_Real x) {
            return !SCIPisIntegral(scip, x) || !SCIPisGT(scip, x, 0);
          }))
        continue;
      // Row is knapsack!!!
      temp_knapsack_rows.push_back(row);
    }
    // Store knapsack constraints
    knapsack_constraints = std::move(temp_knapsack_rows);
  }

  /* Still don't touch this part */

  // Iterate through all knapsack constraints

  for (auto row : knapsack_constraints.value()) {
    // Necessary check for technical reasons.
    if (SCIProwGetNLPNonz(row) != SCIProwGetNNonz(row)) continue;

    // Now the code actually start
    std::span<SCIP_COL* const> row_cols(SCIProwGetCols(row),
                                        SCIProwGetNNonz(row));
    std::span<SCIP_Real> row_coeffs(SCIProwGetVals(row), SCIProwGetNNonz(row));
    auto n_k = SCIProwGetNNonz(row);

    /*
     * We have a knapsack row so the row is
     * a[0]*x[i_1] + a[1]*x[i_2] + ... + a[k]*x[i_k] <= b
     * where a >= 0, {i_1, ..., i_k} c {1,..., n}
     * to simplify, we can pretend the row is of the form
     * a[0]*x[0] + ... + a[k]*x[k] <= b
     * and only do the mapping back at the end
     * the vector a is already available as row_coefs
     */
    // Get the vector x
    std::vector<SCIP_Real> lp_solution;
    lp_solution.reserve(n_k);
    for (auto col : row_cols) lp_solution.push_back(SCIPcolGetPrimsol(col));
    // Get the mapping m -> i_m
    std::vector<int> indices;
    indices.reserve(n_k);
    for (auto col : row_cols) indices.push_back(SCIPcolGetLPPos(col));
    auto rhs = SCIProwGetRhs(row) - SCIProwGetConstant(row);
    // Items taken should contain your cover
    std::vector<int> items_taken;

    // Start here :)
    // Check that the row is not trivially satisfiable(by setting all xs 0 to 1
    // otherwise continue)
    auto total_weight = 0.0;
    for (auto coef : row_coeffs) total_weight += coef;
    if (SCIPisFeasLE(scip, total_weight, rhs)) continue;

    // for our purposes now the row is ax_0 + ax_1 + ... + ax_n <= b, where
    // b= rhs and all x_i are binary

    auto accumulated_weight = 0.0;
    auto accumulated_cost = 0.0;

    // We do multiple past for clarity
    // First past get all LP solutions where x_j = 1

    // For each of the remaining vars
    // We construct a vector of pairs (cost per unit of weight, indices)
    // Hint: std::vector<std::pair<SCIP_Real,int>> list;
    // list.push_back(std::make_pair(quantity, index));
    std::vector<std::pair<SCIP_Real, int>> cost_per_unit_of_weight_indices;
    for (int i = 0; i < n_k; ++i) {
      if (SCIPisFeasEQ(scip, lp_solution[i], 1)) continue;
      auto cost = 1 - lp_solution[i];
      auto weight = row_coeffs[i];
      cost_per_unit_of_weight_indices.push_back(
          std::make_pair(cost / weight, i));
    }

    // Sort using std::ranges::sort by default ranges is sorted in ascending
    // order it is sorted by first checking the first argument

    std::ranges::sort(cost_per_unit_of_weight_indices);

    // start grabbing while accumulated_weight is not greater than b+1
    // Use auto [ratio, index] = vars[i]; to extract data from pair

    // ...
    for (auto [ratio, index] : cost_per_unit_of_weight_indices) {
      if (SCIPisFeasGT(scip, accumulated_weight, rhs + 1)) break;
      accumulated_weight += row_coeffs[index];
      accumulated_cost += 1 - lp_solution[index];
      items_taken.push_back(index);
    }

    // We might overshoot so see if item could have been removed
    /* less efficient
    std::vector<int> items_taken_copy;
    for (int i = std::ssize(items_taken) - 1; i >= 0; --i) {
      auto index = items_taken[i];
      if (SCIPisFeasGT(scip, accumulated_weight - row_coeffs[index], rhs + 1)) {
        accumulated_weight -= row_coeffs[index];
        accumulated_cost -= 1 - incumbent_solution[index];
      } else {
        items_taken_copy.push_back(index);
      }
    }
    items_taken = std::move(items_taken_copy);
    */

    std::erase_if(items_taken, [&](int i) {
      if (SCIPisGE(scip, accumulated_weight - row_coeffs[i], rhs + 1)) {
        accumulated_weight -= row_coeffs[i];
        accumulated_cost -= 1 - lp_solution[i];
        return true;
      }
      return false;
    });

    // No cut can be separated (at least as we heuristically can see it from
    // this constraint)

    if (SCIPisFeasGE(scip, accumulated_cost, 1)) continue;

    // End of code you need to fill

    // The following code push the knapsack cover cut to SCIP
    std::vector<SCIP_Real> coeffs(std::ssize(items_taken), 1);
    std::vector<int> original_indices;
    original_indices.reserve(std::ssize(items_taken));
    for (auto item : items_taken) original_indices.push_back(indices[item]);

    add_cut(scip, sepa, original_indices, coeffs, -SCIPinfinity(scip),
            std::size(items_taken) - 1, true, result);
  }
  return SCIP_OKAY;
}
