 /*@@
   @file      Startup.c
   @date      23 Sep 2019
   @author    E. Jimenez
   @desc
   Register banner - a straight copy of WaveToy.
   @enddesc
 @@*/


#include "cctk.h"

static const char *rcsid = "$Header$";

//CCTK_FILEVERSION(CactusExamples_WaveMoL_Startup_c)

int IDBrillMoL_Startup(void);

int IDBrillMoL_Startup(void)
{

   const char *banner = "IDBrillMoL: Quick and dirty elliptic solver for the conformal factor of Brill waves";

   CCTK_RegisterBanner(banner);

   return 0;
}
