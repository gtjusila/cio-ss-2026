#include "branch_driebeekpenalties.hpp"

#include <scip/pub_lp.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <ranges>
#include <span>
#include <vector>

#include "utils.hpp"
enum SourceType { Col, Row };
enum BaseStat { Basic, AtUpBound, AtLowBound, AtZero };

static BaseStat convert_base_status(SCIP_BASESTAT status,
                                    bool reverse_bounds = false) {
  switch (status) {
    case SCIP_BASESTAT_BASIC:
      return Basic;
    case SCIP_BASESTAT_UPPER:
      return reverse_bounds ? AtLowBound : AtUpBound;
    case SCIP_BASESTAT_LOWER:
      return reverse_bounds ? AtUpBound : AtLowBound;
    case SCIP_BASESTAT_ZERO:
      return AtZero;
  }
  return AtZero;
}

static std::vector<SCIP_Real> get_tableau_row(
    SCIP* scip, std::span<SCIP_COL* const> lp_cols,
    std::span<SCIP_ROW* const> lp_rows, std::span<const int> basic_indices,
    int basic_idx) {
  auto num_cols = static_cast<int>(lp_cols.size());
  auto num_rows = static_cast<int>(lp_rows.size());
  assert(basic_idx >= 0 && basic_idx < num_cols + num_rows);

  auto basis_position = std::ranges::find(basic_indices, basic_idx);
  assert(basis_position != basic_indices.end());
  auto basis_row = static_cast<int>(basis_position - basic_indices.begin());

  std::vector<SCIP_Real> tableau_row(num_cols + num_rows);
  CALL_CHECK(SCIPgetLPBInvRow(scip, basis_row, tableau_row.data() + num_cols,
                              nullptr, nullptr));
  CALL_CHECK(SCIPgetLPBInvARow(scip, basis_row, tableau_row.data() + num_cols,
                               tableau_row.data(), nullptr, nullptr));
  return tableau_row;
}

static void branch_on_column(SCIP* scip, SCIP_RESULT* result,
                             std::span<SCIP_VAR* const> lp_variables,
                             int column_index, SCIP_Real down_child_lower_bound,
                             SCIP_Real up_child_lower_bound,
                             SCIP_Real down_child_priority,
                             SCIP_Real up_child_priority) {
  if (column_index < 0) {
    *result = SCIP_DIDNOTFIND;
    return;
  }

  assert(column_index < std::ssize(lp_variables));
  *result = SCIP_BRANCHED;
  SCIP_NODE* down_node;
  SCIP_NODE* up_node;
  CALL_CHECK(SCIPbranchVar(scip, lp_variables[column_index], &down_node,
                           nullptr, &up_node));

  if (down_node != nullptr) {
    CALL_CHECK(
        SCIPupdateNodeLowerbound(scip, down_node, down_child_lower_bound));
    CALL_CHECK(SCIPchgChildPrio(scip, down_node, down_child_priority));
  }
  if (up_node != nullptr) {
    CALL_CHECK(SCIPupdateNodeLowerbound(scip, up_node, up_child_lower_bound));
    CALL_CHECK(SCIPchgChildPrio(scip, up_node, up_child_priority));
  }
}

SCIP_DECL_BRANCHEXECLP(DriebeekPenalties::scip_execlp) {
  /* Do not edit, SETUP */
  *result = SCIP_DIDNOTRUN;
  if (SCIPgetLPSolstat(scip) != SCIP_LPSOLSTAT_OPTIMAL) return SCIP_OKAY;
  if (!SCIPisLPSolBasic(scip)) return SCIP_OKAY;

  std::vector<int> scip_basis_indices(SCIPgetNLPRows(scip));
  {
    CALL_CHECK(SCIPgetLPBasisInd(scip, scip_basis_indices.data()));
  }

  std::span<SCIP_COL* const> lp_cols(SCIPgetLPCols(scip), SCIPgetNLPCols(scip));
  std::span<SCIP_ROW* const> lp_rows(SCIPgetLPRows(scip), SCIPgetNLPRows(scip));

  auto num_rows = std::ssize(lp_rows);
  auto num_cols = std::ssize(lp_cols);
  // Store basic variables as indices in the extended [columns | row slacks]
  // space while preserving their simplex basis-row order.
  auto basic_indices =
      scip_basis_indices | std::views::transform([num_cols](auto idx) {
        return static_cast<int>(idx >= 0 ? idx : num_cols + (-idx - 1));
      }) |
      std::ranges::to<std::vector<int>>();
  auto lp_variables =
      lp_cols |
      std::views::transform([](auto col) { return SCIPcolGetVar(col); }) |
      std::ranges::to<std::vector<SCIP_VAR*>>();
  auto is_integer_constrained =
      lp_cols | std::views::transform([](auto col) {
        return static_cast<bool>(SCIPcolIsIntegral(col));
      }) |
      std::ranges::to<std::vector<bool>>();
  is_integer_constrained.reserve(num_cols + num_rows);
  for (auto row : lp_rows)
    is_integer_constrained.push_back(SCIProwIsIntegral(row) &&
                                     !SCIProwIsModifiable(row));
  auto lp_solution =
      lp_cols |
      std::views::transform([](auto col) { return SCIPcolGetPrimsol(col); }) |
      std::ranges::to<std::vector<SCIP_Real>>();
  lp_solution.reserve(num_cols + num_rows);
  for (auto row : lp_rows)
    lp_solution.push_back(SCIProwGetConstant(row) -
                          SCIPgetRowLPActivity(scip, row));

  auto lower_bounds =
      lp_cols |
      std::views::transform([](auto col) { return SCIPcolGetLb(col); }) |
      std::ranges::to<std::vector<SCIP_Real>>();
  lower_bounds.reserve(num_cols + num_rows);
  for (auto row : lp_rows)
    lower_bounds.push_back(SCIProwGetConstant(row) - SCIProwGetRhs(row));

  auto upper_bounds =
      lp_cols |
      std::views::transform([](auto col) { return SCIPcolGetUb(col); }) |
      std::ranges::to<std::vector<SCIP_Real>>();
  upper_bounds.reserve(num_cols + num_rows);
  for (auto row : lp_rows)
    upper_bounds.push_back(SCIProwGetConstant(row) - SCIProwGetLhs(row));

  std::vector<SourceType> var_source(num_cols, Col);
  var_source.resize(num_cols + num_rows, Row);
  std::vector<BaseStat> base_stats;
  base_stats.reserve(num_cols + num_rows);
  for (auto col : lp_cols)
    base_stats.push_back(convert_base_status(SCIPcolGetBasisStatus(col)));
  for (auto row : lp_rows)
    base_stats.push_back(convert_base_status(SCIProwGetBasisStatus(row), true));

  std::vector<SCIP_Real> reduced_costs(num_cols + num_rows, 0.0);
  for (auto i : std::views::iota(0, static_cast<int>(num_cols))) {
    if (base_stats[i] == AtUpBound)
      reduced_costs[i] = MIN(SCIPgetColRedcost(scip, lp_cols[i]), 0.0);
    else if (base_stats[i] == AtLowBound)
      reduced_costs[i] = MAX(SCIPgetColRedcost(scip, lp_cols[i]), 0.0);
  }
  for (auto r : std::views::iota(0, static_cast<int>(num_rows))) {
    auto i = num_cols + r;
    if (base_stats[i] == AtUpBound)
      reduced_costs[i] = MIN(-SCIProwGetDualsol(lp_rows[r]), 0.0);
    else if (base_stats[i] == AtLowBound)
      reduced_costs[i] = MAX(-SCIProwGetDualsol(lp_rows[r]), 0.0);
  }

  auto lp_optimal_value = SCIPgetLPObjval(scip);
  auto cutoff = SCIPgetCutoffbound(scip);
  auto have_cutoff = !SCIPisInfinity(scip, cutoff);

  auto best_col = -1;
  auto best_penalty = -SCIPinfinity(scip);
  auto best_least_down = SCIPinfinity(scip);
  auto best_least_up = SCIPinfinity(scip);
  /* End of do not edit, SETUP*/

  /*
   *
   * Lets recap what have been provided
   * Our problem
   *
   * min c^T x
   * s.t. LHS <= Ax + constant <= RHS
   * lb <= x <= ub
   *
   * have been converted into
   *
   * min c^Tx
   * s.t. Ax + slack = 0
   * lb <= x <= ub
   * -RHS + constant <= slack <= -LHS + constant
   *
   * So our LP now has num_cols + num_rows variables.
   * Provided is
   * lp_solution,
   * lower_bounds (adjusted),
   * upper_bounds (adjusted),
   * reduced_costs,
   * var_source (row or column),
   * base_stats, and
   * is_integer_constrained
   *
   * lp_variables contains only num_cols variables
   */

  for (auto idx : basic_indices) {
    // Check preconditions
    // skip if var_source is row since we cannot branch on those
    // skip if variable is not integer constrained
    // skip if variable is basic and integer constrained but is already integral
    if (var_source[idx] == Row) continue;

    // Row slacks have no corresponding entry in lp_variables.
    auto var = lp_variables[idx];
    auto sol = lp_solution[idx];
    if (!is_integer_constrained[idx]) continue;
    if (SCIPisFeasIntegral(scip, sol)) continue;

    // Get tableau_coeffcient for the row correspond to the basic variable
    auto tableau_coeff =
        get_tableau_row(scip, lp_cols, lp_rows, basic_indices, idx);

    /*
     * General bound-aware Driebeek penalty:
     *
     *   down_gap = x_i* - floor(x_i*)
     *   up_gap   = ceil(x_i*) - x_i*
     *
     * For each nonbasic variable j:
     *
     *   down_delta[j] =  down_gap / tableau_coeff[j]
     *   up_delta[j]   = -up_gap   / tableau_coeff[j]
     *
     * A movement delta is infeasible if x_j is fixed, if delta > 0 while x_j
     * is at its upper bound, or if delta < 0 while x_j is at its lower bound.
     *
     *   movement_penalty(j, delta) =
     *     infinity,                 if the movement is infeasible
     *     0,                        if x_j has AtZero basis status
     *     reduced_costs[j] * delta, otherwise
     *
     * In this minimization problem, movement_penalty is the minimum amount by
     * which the objective value must worsen (increase) when x_j moves by delta.
     *
     *   P_i_down = min over nonbasic j of
     *              movement_penalty(j, down_delta[j])
     *
     *   P_i_up   = min over nonbasic j of
     *              movement_penalty(j, up_delta[j])
     *
     * Thus, P_i_down and P_i_up estimate the minimum objective deterioration
     * required to reach the down and up branches, respectively.
     *
     *   Task 1. compute P_i_down and P_i_up and store them in least_up_penalty
     * and least_down_penalty
     *
     */
    auto least_up_penalty = SCIPinfinity(scip);
    auto least_down_penalty = SCIPinfinity(scip);

    // TODO: Compute least_up_penalty and least_down_penalty.

    /*
    * Our driebeek penalty give lowerbounds on the LP objective if we branch up
    * and down. SCIP also give cutoff value. If the LP relaxation of the node
    * have a value higher than the cutoff value we can prune the node
    */
    auto up_pruneable =
        SCIPisInfinity(scip, least_up_penalty) ||
        (have_cutoff &&
         SCIPisGE(scip, lp_optimal_value + least_up_penalty, cutoff));
    auto down_pruneable =
        SCIPisInfinity(scip, least_down_penalty) ||
        (have_cutoff &&
         SCIPisGE(scip, lp_optimal_value + least_down_penalty, cutoff));

    if (up_pruneable && down_pruneable) {
      *result = SCIP_CUTOFF;
      return SCIP_OKAY;
    } else if (up_pruneable) {
      CALL_CHECK(SCIPchgVarUb(scip, var, SCIPfeasFloor(scip, sol)));
      *result = SCIP_REDUCEDDOM;
      return SCIP_OKAY;
    } else if (down_pruneable) {
      CALL_CHECK(SCIPchgVarLb(scip, var, SCIPfeasCeil(scip, sol)));
      *result = SCIP_REDUCEDDOM;
      return SCIP_OKAY;
    } else {
      auto penalty = MAX(least_up_penalty, least_down_penalty);
      if (penalty > best_penalty) {
        best_penalty = penalty;
        best_col = idx;
        best_least_up = least_up_penalty;
        best_least_down = least_down_penalty;
      }
    }
  }

  // Call branch on column function
  // branch_on_column(
  //     SCIP* scip, SCIP_RESULT* result,
  //     std::span<SCIP_VAR* const> lp_variables,
  //     int column_index, column to branch on
  //     SCIP_Real down_child_lower_bound, // A lower bound on the LP optimal
  //     value of the down child
  //     SCIP_Real up_child_lower_bound, // A lower bound on the LP optimal
  //     value of the down child
  //     SCIP_Real down_child_priority, // Value to
  //     determine if down child or up child will be visited first
  //     SCIP_Real up_child_priority // same here
  //     );
  if (best_col >= 0) {
    branch_on_column(scip, result, lp_variables, best_col,
                     lp_optimal_value + best_least_down,
                     lp_optimal_value + best_least_up, -best_least_down,
                     -best_least_up);
  } else {
    *result = SCIP_DIDNOTFIND;
  }

  return SCIP_OKAY;
}
