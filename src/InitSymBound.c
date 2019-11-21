/* Set the symmetries for the conformal factor */


#include "cctk.h"
#include "cctk_Arguments.h"

#include "Symmetry.h"

static const char *rcsid = "$Header$";

CCTK_FILEVERSION(CactusWave_IDBrillMoL_InitSymBound_c)

void IDBrillMoL_InitSymBound(CCTK_ARGUMENTS);

void IDBrillMoL_InitSymBound(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;

  int sym[3];

  sym[0] = 1;
  sym[1] = 1;
  sym[2] = 1;

  SetCartSymVN(cctkGH, sym,"IDBrillMoL::phi");
  SetCartSymVN(cctkGH, sym,"IDBrillMoL::pi");

  SetCartSymVN(cctkGH,sym, "IDBrillMoL::rhs_phi" );
  SetCartSymVN(cctkGH,sym, "IDBrillMoL::rhs_pi" );

  return;
}
