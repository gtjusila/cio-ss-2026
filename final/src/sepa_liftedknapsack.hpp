#pragma once
#include <objscip/objsepa.h>
#include <scip/type_retcode.h>

#include <optional>
#include <span>
#include <vector>

class LiftedKnapsackSepa : public scip::ObjSepa {
 public:
  LiftedKnapsackSepa(SCIP* scip)
      : scip::ObjSepa(scip, "liftedknapsack", "Lifted Knapsack Separator",
                      99999, 1, 1.0, 0, false) {};
  SCIP_DECL_SEPAEXECLP(scip_execlp) override;

 private:
  void add_cut(SCIP* scip, SCIP_SEPA* sepa, std::span<int> indices,
               std::span<SCIP_Real> coefficients, SCIP_Real lhs, SCIP_Real rhs,
               bool forcecut, SCIP_Result* result);
  int cut_count = 0;
  std::optional<std::vector<SCIP_Row*>> knapsack_constraints = std::nullopt;
};
