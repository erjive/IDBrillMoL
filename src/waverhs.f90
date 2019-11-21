#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

    implicit none

    subroutine setupbrilldata2D(CCTK_ARGUMENTS)

c     Set up axisymmetric Brill data for elliptic solve. The elliptic
c     equation we need to solve is:
c
c     __                 2      2
c     \/  psi  +  psi ( d q  + d   q ) / 4  =  0
c       f                z      rho
c
c     where:
c
c     __
c     \/  =  Flat space Laplacian.
c       f

      implicit none

      DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_PARAMETERS

      integer i,j,k

      CCTK_REAL x1,y1,z1,rho1
      CCTK_REAL brillq,eps
      CCTK_REAL zp,zm,rhop,rhom
      CCTK_REAL zero,one

      external brillq

c     Numbers.

      zero = 0.0D0
      one  = 1.0D0

c  /* Set up shorthands */
      dx = CCTK_DELTA_SPACE(0);
      dy = CCTK_DELTA_SPACE(1);
      dz = CCTK_DELTA_SPACE(2);
      dt = CCTK_DELTA_TIME;

      dx2 = dx*dx;
      dy2 = dy*dy;
      dz2 = dz*dz;
      dt2 = dt*dt;

      dx2i = one/dx2;
      dy2i = one/dy2;
      dz2i = one/dz2;

      do k=1,cctk_lsh(3)-1
         do j=1,cctk_lsh(2)-1
            do i=1,cctk_lsh(1)-1

               x1 = x(i,j,k)
               y1 = y(i,j,k)
               z1 = z(i,j,k)

               rho1 = sqrt(x1*x1 + y1*y1)

               brillpsi(i,j,k) = one

               if (CCTK_EQUALS(q_function,"exp")) then

                  brillsource(i,j,k) = 0.0

               else if (CCTK_EQUALS(q_function,"eppley")) then
                  brillsource(i,j,k) = 0.0
               else if (CCTK_EQUALS(q_function,"gundlach")) then

                  brillsource(i,j,k) = (0.5*gundlacha/(gundlachsigmar**4*gundlachsrho**2) * &
                                       Exp((rho1**2 + z1**2 -gundlachr0**2 )/gundlachsigmar**2) * &
                                       (gundlachsigmar**4 - 6.*gundlachsigmar**2*rho1**2 + &
                                       2.*rho1**4 + 2.*rho1**2*z1**2))*phi(i,j,k)

               else
                  call CCTK_WARN(0,"Brill wave data type not recognised")
               end if

c              Add the laplacian to the r.h.s of the wave equation.

               brillsource(i,j,k) = brillsource(i,j,k)+ &
                                    dx2i*(phi(i-1,j  ,k)-2.0*phi(i,j,k)+phi(i+1,j,k))+&
                                    dy2i*(phi(i  ,j-1,k)-2.0*phi(i,j,k)+phi(i  ,j+1,k))+&
                                    dz2i*(phi(i,j,k-1)-2.0*phi(i,j,k)+phi(i,j,k+1))

            end do
         end do
      end do

