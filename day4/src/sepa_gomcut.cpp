#include "sepa_gomcut.hpp"

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

void GomCutSepa::add_cut(SCIP* scip, SCIP_SEPA* sepa,
                                 std::span<int> indices,
                                 std::span<SCIP_Real> coefficients,
                                 SCIP_Real lhs, SCIP_Real rhs, bool forcecut,
                                 SCIP_Result* result) {
  RowDeleter deleter(scip);
  RowPtr new_row(nullptr, deleter);
  CALL_CHECK(SCIPcreateEmptyRowSepa(
      scip, std::out_ptr(new_row), sepa,
      std::format("gom_cut_{}", cut_count++).c_str(), lhs, rhs, false,
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

SCIP_DECL_SEPAEXECLP(GomCutSepa::scip_execlp) {
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

  // Get LProws, LPcols pointer
  std::span<SCIP_Row* > lp_rows(SCIPgetLPRows(scip), SCIPgetNLPRows(scip));
  std::span<SCIP_Col* > lp_cols(SCIPgetLPCols(scip), SCIPgetNLPCols(scip));

  // not sure whether vector is acceptable for SCIPgetLPBasisInd
  std::vector<int> basis_indices(std::ssize(lp_rows));

  // I need later for BInvA
  int n_cols = SCIPgetNLPCols(scip);
  std::vector<SCIP_Real> row_coefs(n_cols);
  std::vector<int> row_inds(n_cols);

  // instead of &, using data return the pointer
  SCIPgetLPBasisInd(scip, basis_indices.data());

  for (auto [idx, col] : std::views::enumerate(lp_rows)) {


    int basis_idx = basis_indices[idx];
    if (basis_idx < 0) { continue; }

    auto basic_col = lp_cols[basis_idx];
    auto basic_var = SCIPcolGetVar(basic_col);

    if (SCIPvarGetType(basic_var) == SCIP_VARTYPE_CONTINUOUS) { continue; }
    SCIP_Real basic_value = SCIPcolGetPrimsol(basic_col);

    if (SCIPisFeasIntegral(scip, basic_value)) {continue;}
    int n_inds = 0;

    /* arg Parameters
    scip	SCIP data structure
    r	row number
    binvrow	row in B^-1 from prior call to SCIPgetLPBInvRow(), or NULL
    coefs	array to store the coefficients of the row
    inds	array to store the non-zero indices, or NULL
    ninds	pointer to store the number of non-zero indices, or NULL (-1: if we do not store sparsity informations)
    */
    SCIPgetLPBInvARow(scip, idx, nullptr, row_coefs.data(), row_inds.data(), &n_inds);

    auto frac_basic_value = SCIPfrac(scip, basic_value);
    if (SCIPisFeasEQ(scip, frac_basic_value, 0.0)) { continue; }

    std::vector<int> cut_indices;
    std::vector<SCIP_Real> cut_coeffs;

    for (int k = 0; k < n_inds; ++k) {
        // get the data we prev requested
        int col_idx = row_inds[k];
        if (col_idx == basis_idx) { continue; }

        SCIP_Real coef = row_coefs[k];
        SCIP_Real frac_coef = SCIPfrac(scip, coef);

        if (SCIPisFeasEQ(scip, frac_coef, 0.0)) { continue; }

        // auto non_basic_col = lp_cols[col_idx];
        // auto non_basic_var = SCIPcolGetVar(non_basic_col);

        cut_indices.push_back(col_idx);
        cut_coeffs.push_back(frac_coef);
    }
  // The following code push the knapsack cover cut to SCIP
    if (!cut_indices.empty()) {
      add_cut(scip, sepa, cut_indices, cut_coeffs,
            frac_basic_value, SCIPinfinity(scip), false, result);
    }

  }
  return SCIP_OKAY;
}
