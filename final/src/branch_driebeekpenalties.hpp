#pragma once

#include <objscip/objbranchrule.h>
#include <scip/type_branch.h>
#include <scip/type_retcode.h>

class DriebeekPenalties : public scip::ObjBranchrule {
 public:
  DriebeekPenalties(SCIP* scip)
      : scip::ObjBranchrule(scip, "driebeekpenalties", "Driebeek Penalties", 0,
                            -1, 1.0) {};
  ~DriebeekPenalties() = default;
  SCIP_DECL_BRANCHEXECLP(scip_execlp) override;
};