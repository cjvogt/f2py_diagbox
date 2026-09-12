c =====================================================================
c
c     SUBROUTINE ADVECTION2
c 
c     PURPOSE: this routine calculates the contribution of the
c                  advection to the mixed layer budget. Only executed
c                  if sol_scheme is 'adv', it also computes the horizontal
c                  velocity u to close the budget
c                  
c
c     VARIABLES:
c 
c     IN : k              : time loop counter
c          ts             : time step
c          sol_scheme     : solution scheme
c          dc_dh          : horizontal C gradient
c          ddc13_dh       : horizontal dC13 gradient
c          c              : obs c concentration
c          c13            : obs c13 concentraton
c          h              : obs mixed layer depth
c     
c
c     OUT: u              : horizontal velocity
c          D_C            : modified array holding the rates of change for C
c          D_C13          : modified array holding the rates of change for C13
c          fluxes         : modified array holding the fluxes
c 
c     REVISION:
c
c     date     author  remarks
c   
c     17.10.99   ng    first implementation, based on advection.F
c     05.05.02   ng    set u to 0 in case dc_dh = 0

c =====================================================================
c
      subroutine advection2(k,ts,sol_scheme,
     $     dc_dh,dc13_dh,
     $     c,c13,h,
     $     u_var,D_C,D_c13,fluxes)
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
# include "solution.h"
      integer k
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision rates_to_fluxes ! functions
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer i
      double precision sec_day, day_year
      double precision sum_c
      parameter (sec_day = 24.d0 * 3600.d0)
      parameter (day_year = 365.)

c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     calculate horizontal velocity
c --------------------------------------------------------------------
c 
      sum_c = 0.d0
      do i = 1,3
         sum_c = sum_c + D_C(i,k)
      enddo
      sum_c = sum_c + D_C(5,k)
      if (dc_dh .ne. 0.d0) then
         u_var(k) = - ((D_C(7,k) - sum_c)/dc_dh) 
     $        / (ts * sec_day * day_year) ! gives m s-1
      else
          u_var(k) = 0.0d0
       endif
         
c     
c --------------------------------------------------------------------
c     calculate rate of change and fluxes due to advection
c       multiply with ts and sec_day because u is given in m s-1
c       17.10.99 ng added day_year since ts is in years
c     do this only for old standard solution scheme
c --------------------------------------------------------------------
c
      D_C(4,k) = - u_var(k) * dc_dh * ts * sec_day * day_year
      fluxes(4,k) = rates_to_fluxes(D_C(4,k),ts,h)
      D_c13(4,k) = - u_var(k) * dc13_dh * ts * sec_day * day_year
c      
c --------------------------------------------------------------------
c     for debugging purposes
c --------------------------------------------------------------------
c
#ifdef debugging
c
      if (mod(k,90) .eq. 0) then
         write(*,*)
         write(*,*) '--ADVECTION2: '
         write(*,*) '  k : ',k
         write(*,*) '  dc_dh                     : ',dc_dh
         write(*,*) '  u_var(k)                  : ',u_var(k)
         write(*,*) '  D_C(7,k)                  : ',D_C(7,k)
         write(*,*) '  sum_c                     : ',sum_c
         write(*,*)
         write(*,*) '  dc_dt_adv                 : ',D_C(4,k)
         write(*,*) '  D_C13_adv                 : ',D_C13(4,k)
      endif
#endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

