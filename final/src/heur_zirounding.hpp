#pragma once

#include <objscip/objheur.h>
#include <scip/type_heur.h>
#include <scip/type_retcode.h>

class ZIRoundHeur : public scip::ObjHeur {
 public:
  ZIRoundHeur(SCIP* scip)
      : scip::ObjHeur(scip, "ziround", "ZI Rounding Heuristic", 'z', 99999, 1,
                      0, -1, SCIP_HEURTIMING_DURINGLPLOOP, false) {};
  ~ZIRoundHeur() = default;
  SCIP_DECL_HEUREXEC(scip_exec) override;

 private:
  int nlpsolve = -1;
};
