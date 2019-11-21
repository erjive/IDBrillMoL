
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void IDBrillMoL_RegisterVars(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_INT ierr = 0, group, rhs, var;

  /* phi and rhs_phi */
  group   = CCTK_GroupIndex("IDBrillMoL::scalarevolvemol");
  rhs   = CCTK_GroupIndex("IDBrillMoL::rhsscalar");
  ierr += MoLRegisterEvolvedGroup(group, rhs);

  /* pi and rhs_pi */
/*  var   = CCTK_VarIndex("IDBrillMoL::pi");
  rhs   = CCTK_VarIndex("IDBrill::rhs_pi");
  ierr += MoLRegisterEvolved(var, rhs);*/
  /*return ;*/
  if (ierr) CCTK_ERROR("Problems registering with MoL");

}
