#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "math.h"

void setupbrilldata2D(CCTK_ARGUMENTS)

/* Calculate the Brill source of the wave-like equation
for the conformal factor of Brill waves.*/

{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
/*  DECLARE_CCTK_FUNCTIONS;*/

  int i,j,k;
  int index;
  int istart, jstart, kstart, iend, jend, kend;
  CCTK_REAL x1,y1,z1,rho1;
  istart = 0;
  jstart = 0;
  kstart = 0;

  iend = cctk_lsh[0];
  jend = cctk_lsh[1];
  kend = cctk_lsh[2];

  /* Calculate the right hand sides. */

#ifdef DEBUG_MOL
  printf("About to loop.\n");
#endif

  for (k=kstart; k<kend; k++)
  {
    for (j=jstart; j<jend; j++)
    {
      for (i=istart; i<iend; i++)
      {
        index = CCTK_GFINDEX3D(cctkGH,i,j,k);


               x1 = x[index];
               y1 = y[index];
               z1 = z[index];

               rho1 = sqrt(x1*x1 + y1*y1);

               if (CCTK_EQUALS(q_function,"exp"))
               {
                  brillsource[index] = 0.0;
               }
               else if (CCTK_EQUALS(q_function,"eppley"))
               {
                  brillsource[index] = 0.0;
               }
               else if (CCTK_EQUALS(q_function,"gundlach"))
               {
                  brillsource[index] = (0.5*gundlacha/(pow(gundlachsigmar,4)*pow(gundlachsrho,2)) *
                                       expf(-(pow(rho1,2) + pow(z1,2) - pow(gundlachr0,2) )/pow(gundlachsigmar,2)) *
                                       (pow(gundlachsigmar,4) - 6.*pow(gundlachsigmar,2)*pow(rho1,2) +
                                       2.*pow(rho1,4) + 2.*pow(rho1,2)*pow(z1,2)));
               }
               else
               {
                  CCTK_WARN(0,"Brill wave data type not recognised");
               }

              /*printf("%f %f %f %f %s\n",x1,y1,z1,brillsource[index],qfunction);*/

      }
    }
  }

#ifdef DEBUG_MOL
  printf("Done with loop\n");
#endif

  return;

}
