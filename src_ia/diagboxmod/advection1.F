c =====================================================================
c
c     SUBROUTINE ADVECTION1
c 
c     PURPOSE: this routine calculates the contribution of the
c                  advection to the mixed layer budget. Exact behavior
c                  depends on sol_scheme. In the 'std' solution scheme,
c                  everything is calculated; in any other solution scheme
c                  only the horizontal gradients are determined
c                  
c
c     VARIABLES:
c 
c     IN : k              : time loop counter
c          ts             : time step
c          sol_scheme     : solution scheme
c          u              : horizontal velocity
c          dc_dh          : horizontal C gradient
c          ddc13_dh       : horizontal dC13 gradient
c          c              : obs c concentration
c          c13            : obs c13 concentraton
c          h              : obs mixed layer depth
c     
c
c     OUT: dc13_dh        : horizontal gradient
c          D_C            : modified array holding the rates of change for C
c          D_C13          : modified array holding the rates of change for C13
c          fluxes         : modified array holding the fluxes
c 
c     REVISION:
c
c     date     author  remarks
c   
c      2.08.96   ng    first implementation
c      9.09.96   ng    use accurate formula for dc13_dh
c     17.10.99   ng    included sol_scheme, adapting for interannual var.
c                          runs
cc
c =====================================================================
c
      subroutine advection1(k,ts,sol_scheme,
     $     u,dc_dh,ddc13_dh,
     $     c,c13,h,
     $     dc13_dh,
     $     D_C,D_c13,fluxes)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
# include "input.h"
# include "adv_params.h"
# include "c13_params.h"
      integer k
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision rates_to_fluxes, c13_to_delta! functions
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision dc13_o,dc13_dh_old
      double precision sec_day, day_year
      parameter (sec_day = 24.d0 * 3600.d0)
      parameter (day_year = 365.)

c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     calculate dc13_dh from ddc13_dh
c      same formula as for dc13dz, use more accurate formula
c --------------------------------------------------------------------
c
      dC13_o = c13_to_delta(c13,c)
c
c      dc13_dh_old = r_std/(1+r_std) * 
c     $     (c * ddC13_dh*1.d-3 + dc_dh*(dc13_o*1.d-3 + 1.d0))
      dc13_dh = r_std/(1+r_std*dc13_o*1.d-3+r_std) * 
     $     ((dc13_o*1.d-3 + 1.d0)*dc_dh +
     $      c / (1+r_std*dc13_o*1.d-3+r_std) * ddc13_dh * 1.d-3)
c
c      write(*,*) 'k : ',k
c      write(*,*) 'dc13_dh,dc13_dh_old : ',dc13_dh,dc13_dh_old
c
c --------------------------------------------------------------------
c     calculate rate of change and fluxes due to advection
c       multiply with ts and sec_day because u is given in m s-1
c       17.10.99 ng added day_year since ts is in years
c     do this only for old standard solution scheme
c --------------------------------------------------------------------
c
      if (sol_scheme .eq. 'std') then
         D_C(4,k) = - u * dc_dh * ts * sec_day * day_year
         fluxes(4,k) = rates_to_fluxes(D_C(4,k),ts,h)
         D_c13(4,k) = -u * dc13_dh * ts * sec_day * day_year
      else
         D_C(4,k) = 0.0d0
         fluxes(4,k) =  0.0d0
         D_c13(4,k) =  0.0d0
      endif
          
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

