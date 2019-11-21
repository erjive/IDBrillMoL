#include <stddef.h>

#include "stdio.h"
#include "stdlib.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "math.h"
#include "carpet.hh"

void tolerance(CCTK_ARGUMENTS)
  /* Calculate the diference of Pi=dphi/dt between two consecutive iterations
     If we reach the tolerance, stop the wave-like evolution for the conformal factor */
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int vindex;             /* grid variable index */
  int target_proc;        /* processor to hold the result */
  int reduction_handle;   /* handle for reduction operator */
  char *reduction_name;   /* reduction operator to use */

  CCTK_REAL result;       // reduction value of iteration i
  static CCTK_REAL result_p; //reduction value of iteration i-1
  double diff_it;         // difference of mean values at consecutive iterations
  int ierr = -1;

  //int i,j,k;
  //int index;
  //int istart, jstart, kstart, iend, jend, kend;
  CCTK_REAL lres , gres;     //Local and global residual


  if (cctk_iteration>2){

    DECLARE_CCTK_ARGUMENTS;
    /*    istart = 1;
          jstart = 1;
          kstart = 1;

          iend = cctk_lsh[0];
          jend = cctk_lsh[1];
          kend = cctk_lsh[2];
          */

    lres = 0.0;
    gres = 0.0;

    BEGIN_GLOBAL_MODE(cctkGH) {
      ENTER_LEVEL_MODE(cctkGH, 0) {
        BEGIN_LOCAL_MAP_LOOP(cctkGH,CCTK_GF) {
          BEGIN_LOCAL_COMPONENT_LOOP(cctkGH, CCTK_GF) {
            DECLARE_CCTK_ARGUMENTS;
            for(int k=1;k<cctk_lsh[2];k++) {
              for(int j=1;j<cctk_lsh[1];j++) {
                for(int i=1;i<cctk_lsh[0];i++) {
                  int index= CCTK_GFINDEX3D(cctkGH,i,j,k);
                  lres += abs(phi[index]-phi_p[index]);
                }
              }
            }
            gres = lres/((cctk_lsh[0]-1)*(cctk_lsh[1]-1)*(cctk_lsh[2]-1));
          } END_LOCAL_COMPONENT_LOOP;
        } END_LOCAL_MAP_LOOP;
      } LEAVE_LEVEL_MODE;
    } END_GLOBAL_MODE;
    //gres = lres/((iend-1)*(jend-1)*(kend-1));


    //diff_it = fabs(result-result_p);

    if ( gres <= tol_pi )
    {

      printf("Tolerance reached after %d iterations, tolerance=%e \n ",cctk_iteration,gres);
    /*BEGIN_GLOBAL_MODE(cctkGH) {
      ENTER_LEVEL_MODE(cctkGH, 0) {
      ierr = CCTK_OutputVarAs(cctkGH,"IDBrillMoL::phi","brillpsi");
      } LEAVE_LEVEL_MODE;
    } END_GLOBAL_MODE;
    */
      CCTK_TerminateNext (cctkGH);
    }



    /*Save the mean value of Pi of the current iteration */
    //result_p = result;

  }
}
