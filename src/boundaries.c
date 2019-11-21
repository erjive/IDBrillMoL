#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

/*
#ifndef DEBUG_MOL
#define DEBUG_MOL
#endif
*/

static const char *rcsid = "$Header$";

void IDBrillMoL_Boundaries(CCTK_ARGUMENTS);

void IDBrillMoL_Boundaries(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
 /* DECLARE_CCTK_FUNCTIONS;*/

  CCTK_INT ierr;

  CCTK_INT one;

  /*shorthands*/
  one = 1;

  /*According to ProcaEvolve Thorn, the radiative Bcs are applied
    with the NewRad infraestructure. To enforce the symmetry Bcs,
    we should register all the Bcs as 'none'.*/

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,
       "IDBrillMoL::scalarevolvemol", "none");
  if (ierr < 0)
  {
       CCTK_WARN(0, "Failed to register BC for IDBrillMoL::phi,pi!");
  }
  /*ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,
       "IDBrillMoL::pi", "none");
  if (ierr < 0)
  {
       CCTK_WARN(0, "Failed to register BC for IDBrillMoL::pi!");
  }
*/
  return;
}
