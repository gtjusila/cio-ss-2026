#pragma once
#include <objscip/objsepa.h>
#include <scip/type_retcode.h>

#include <optional>
#include <span>
#include <vector>

class GomCutSepa : public scip::ObjSepa {
 public:
  GomCutSepa(SCIP* scip)
      : scip::ObjSepa(scip, "gomcut", "Gomory Cut Separator",
                      99999, 1, 1.0, 0, false) {};
  SCIP_DECL_SEPAEXECLP(scip_execlp) override;

 private:
  void add_cut(SCIP* scip, SCIP_SEPA* sepa, std::span<int> indices,
               std::span<SCIP_Real> coefficients, SCIP_Real lhs, SCIP_Real rhs,
               bool forcecut, SCIP_Result* result);
  int cut_count = 0;
};
