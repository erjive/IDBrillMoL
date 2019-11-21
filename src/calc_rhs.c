#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "math.h"
#include "Boundary.h"
void calc_rhs(CCTK_ARGUMENTS)

/* Calculate the rhs of the wave-like equation
for the conformal factor of Brill waves.*/

{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
/*  DECLARE_CCTK_FUNCTIONS;*/

  int i,j,k;
  int index;
  int istart, jstart, kstart, iend, jend, kend;
  CCTK_REAL dx,dy,dz,dx2,dy2,dz2;
  CCTK_REAL dx2i,dy2i,dz2i;
  CCTK_REAL itwelfth,fthird,fivehalf;

  /* Set up shorthands */
  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);

  dx2 = dx*dx;
  dy2 = dy*dy;
  dz2 = dz*dz;

  dx2i = 1.0/dx2;
  dy2i = 1.0/dy2;
  dz2i = 1.0/dz2;


  istart = 1;
  jstart = 1;
  kstart = 1;

  iend = cctk_lsh[0];
  jend = cctk_lsh[1];
  kend = cctk_lsh[2];

  itwelfth = 1.0/12.0;
  fthird   = 4.0/3.0;
  fivehalf = 5.0/2.0;


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

              /*Add the laplacian to the r.h.s of the wave equation.*/
              /*This is only second order accurate. */

               rhs_pi[index] = brillsource[index]*phi[index]+
                                   /* dx2i*(phi[CCTK_GFINDEX3D(cctkGH, i-1, j  , k  )]-2.0*phi[index]+phi[CCTK_GFINDEX3D(cctkGH, i+1, j  , k)])+
                                    dy2i*(phi[CCTK_GFINDEX3D(cctkGH, i  , j-1, k  )]-2.0*phi[index]+phi[CCTK_GFINDEX3D(cctkGH, i  , j+1, k)])+
                                    dz2i*(phi[CCTK_GFINDEX3D(cctkGH, i  , j  , k-1)]-2.0*phi[index]+phi[CCTK_GFINDEX3D(cctkGH, i  , j  , k+1)])-*/
                               dx2i*(-itwelfth*phi[CCTK_GFINDEX3D(cctkGH, i-2, j  , k  )] +fthird*phi[CCTK_GFINDEX3D(cctkGH, i-1, j  , k  )]
                                     -itwelfth*phi[CCTK_GFINDEX3D(cctkGH, i+2, j  , k  )] +fthird*phi[CCTK_GFINDEX3D(cctkGH, i+1, j  , k  )]
                                     -fivehalf*phi[index])+
                               dy2i*(-itwelfth*phi[CCTK_GFINDEX3D(cctkGH, i, j-2 , k  )] +fthird*phi[CCTK_GFINDEX3D(cctkGH, i, j-1  , k  )]
                                     -itwelfth*phi[CCTK_GFINDEX3D(cctkGH, i, j+2 , k  )] +fthird*phi[CCTK_GFINDEX3D(cctkGH, i, j+1  , k  )]
                                     -fivehalf*phi[index])+
                               dz2i*(-itwelfth*phi[CCTK_GFINDEX3D(cctkGH, i, j , k-2  )] +fthird*phi[CCTK_GFINDEX3D(cctkGH, i, j  , k-1  )]
                                     -itwelfth*phi[CCTK_GFINDEX3D(cctkGH, i, j , k+2  )] +fthird*phi[CCTK_GFINDEX3D(cctkGH, i, j  , k+1  )]
                                     -fivehalf*phi[index])-
                               k_damping*pi[index] ;
               /* We consider the wave equation as a first order system,
               where dphi/dt = pi  and d pi/dt = brillsource*phi+nabla phi */

               rhs_phi[index] = pi[index];
/*               rhs_pi [index] = brillsource[index];*/
/*              printf("%f %f %f %f\n",x[index],y[index],z[index],phi[CCTK_GFINDEX3D(cctkGH, i-1, j  , k  )]);*/

#ifdef DEBUG_MOL
        printf("ijk %i %i %i, brillsource %f\n",i,j,k,brillsource[index]);
#endif
      }
    }
  }

#ifdef DEBUG_MOL
  printf("Done with loop\n");
#endif

  return;

}


void calc_rhs_bdry(CCTK_ARGUMENTS)

/* Apply the radiative boundary conditions with the Newrad thorn interface.*/
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
/*  DECLARE_CCTK_FUNCTIONS;*/
  int sw[3];
  int ierr=-1;

  /*Stencil width */

  sw[0]=1;
  sw[1]=1;
  sw[2]=1;

/*  ierr = NewRad_Apply(cctkGH, pi, rhs_pi, 1.0, 1.0, n_brill);*/
    ierr = NewRad_Apply(cctkGH, pi, rhs_pi, 0.0, 1.0, n_brill);
/*    ierr = BndRadiativeVN(cctkGH,sw,0.0,1.0,"IDBrillMoL::pi","IDBrillMoL::pi");*/
/*    ierr = BndRadiativeDirVN(cctkGH,1,1,0.0,1.0,"IDBrillMoL::pi","IDBrillMoL::pi");*/
/*    ierr = BndRadiativeDirVN(cctkGH,1,1,0.0,1.0,"IDBrillMoL::pi","IDBrillMoL::pi");*/

/*    ierr = BndRadiativeDirVN(cctkGH,sw,3,0.0,1.0,"IDBrillMoL::pi","IDBrillMoL::pi");*/
  if (ierr<0)
  {
      CCTK_WARN(0,"Boundary conditions not applied - Radiative Bcs not applied!");
  }
}
