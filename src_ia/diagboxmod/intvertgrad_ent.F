c =====================================================================
c
c     SUBROUTINE INTVERTGRAD_ENT
c 
c     PURPOSE: this routine calculates the contribution of entrainment
c                  on the mixed layer budget. This routine is the
c                  version where the integrated vertical gradient 
c                  scheme developed by ng is used. 
c                  
c
c     VARIABLES:
c 
c     IN : k              : time loop counter
c          time           : actual time
c          ts             : time step
c          const_vertgrad : logical indicating use of constant vert grad
c          dCdz_const     : value of constant vertical C gradient
c          dCdz_var       : value of variable vertical C gradient
c          D_dcdz         : correction term for dc/dz
c          ddC13_dC_bml   : ratio of dd13C/dz to dC/dz below mixed layer
c          ent_fact       : entrainment multiplication factor
c
c          h_new          : mixed layer depth at timestep t+ts
c          h_old          : mixed layer depth at tiemstep t
c
c          C              : obs C in the mixed layer
c          C13            : obs C13 in the mixed layer
c          h              : obs mixed layer depth
c
c     OUT: D_C            : modified array holding the rates of change for C
c          D_C13          : modified array holding the rates of change for C13
c          fluxes         : modified array holding the fluxes
c 
c     REVISION:
c
c     date     author  remarks
c   
c      2.08.96   ng    first implementation
c      9.09.96   ng    use accurate formula for dc13_dz
c     17.10.99   ng    adapted for interannual variability runs
c
c =====================================================================
c
      subroutine intvertgrad_ent(k,time,ts,
     $     const_vertgrad,dCdz_const,dCdz_var,D_dcdz,
     $     ddC13_dC_bml,ent_fact,
     $     h_new,h_old,
     $     c,c13,h,
     $     D_C,D_c13,fluxes)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
# include "input.h"
# include "diffent_params.h"
# include "c13_params.h"
      integer k
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision rates_to_fluxes,c13_to_delta ! functions
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision dc13_o,ddc13dz,dcdz,dc13dz
      double precision integral_dcdz,integral_dc13dz
      double precision c_zero,c13_zero,c_bml,c13_bml
      logical first
      data first /.true./
c
c ---------------------------------------------------------------------
c     variables to be saved between calls
c ---------------------------------------------------------------------
c
      save first,integral_dcdz
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     determine vertical C gradient
c --------------------------------------------------------------------
c
      if (const_vertgrad) then
         dCdz = dCdz_const
      else
         dCdz = dCdz_var + D_dcdz
      endif
c
c --------------------------------------------------------------------
c     calulate the vertical C13 gradient from the vertical
c       dC13 gradient (approximate formula)
c       units [umol kg-1 m-1]
c     changed to accurate formula
c --------------------------------------------------------------------
c    
      dC13_o = c13_to_delta(c13,c)
      ddC13dz = dCdz * ddC13_dC_bml ! unit [per mil m-1]
c      dc13dz = r_std/(1+r_std) * 
c     $     (c * ddC13dz*1.d-3 + dcdz*(dc13_o*1.d-3 + 1.d0))

      dc13dz = r_std/(1+r_std*dc13_o*1.d-3+r_std) * 
     $     ((dc13_o*1.d-3 + 1.d0)*dcdz +
     $      c / (1+r_std*dc13_o*1.d-3+r_std) * ddc13dz * 1.d-3)
c
c --------------------------------------------------------------------
c     determine integral 
c --------------------------------------------------------------------
c
#ifndef backw_sim
      if (h_new .ge. h_old) then
#else
      if (h_new .le. h_old) then
#endif
         if (first) then
            integral_dcdz = 0.d0
            integral_dc13dz = 0.d0
            c_zero = c
            c13_zero = c13
c            c_zero = 2026            ! value at minimum
c            c13_zero = 22.5467       ! value at minimum
            first = .false.
         else
            integral_dcdz = integral_dcdz + 
     $           (h_new-h_old)*dcdz
            integral_dc13dz = integral_dc13dz + 
     $           (h_new-h_old)*dc13dz
         endif
      else
         if (.not. first) then
            integral_dcdz = 0.d0
            integral_dc13dz = 0.d0
            c_zero = 0.d0
            c13_zero = 0.d0
            c_bml = 0.d0
            c13_bml = 0.d0
          endif       
          first = .true.
      endif
c
c --------------------------------------------------------------------
c     calculate contribution only during depeening of the 
c       mixed layer
c --------------------------------------------------------------------
c  
#ifndef backw_sim
      if (h_new .ge. h_old) then
#else
      if (h_new .le. h_old) then
#endif
c
c --------------------------------------------------------------------
c     determine the different contributions to entrainment
c --------------------------------------------------------------------
c
         c_bml = c_zero + integral_dcdz
         c13_bml = c13_zero + integral_dc13dz
c
         write(30,*) k,c_bml,c13_bml
c
         if (c_bml .le. c) then
         endif
c
         D_c(3,k) = ent_fact * 1.d0/h * (h_new - h_old) * (c_bml - c)
         fluxes(3,k) = rates_to_fluxes(D_C(3,k),ts,h)
         D_c13(3,k) = ent_fact * 1.d0/h * (h_new - h_old) * 
     $        (c13_bml - c13)
c
c --------------------------------------------------------------------
c     shoaling period
c --------------------------------------------------------------------
c
      else
         D_c(3,k) = 0.d0
         fluxes(3,k) = 0.d0
         D_c13(3,k) = 0.d0
      endif
c
c --------------------------------------------------------------------
c     for debugging purposes
c --------------------------------------------------------------------
c
c
      if (mod(k,30) .eq. 0) then
         write(*,*)
         write(*,*) 'k : ',k
         write(*,*) 'first        : ',first
         write(*,*) 'h_new, h_old : ',h_new, h_old
         write(*,*) 'dcdz         : ',dcdz
         write(*,*) 'dc13dz       : ',dc13dz
         write(*,*) 'c_zero       : ',c_zero
         write(*,*) 'c13_zero     : ',c13_zero
         write(*,*) 'int_dcdz     : ',integral_dcdz
         write(*,*) 'int_dc13dz   : ',integral_dc13dz
         write(*,*) 'c_bml        : ',c_bml
         write(*,*) 'c            : ',c
         write(*,*) 'c13_bml      : ',c13_bml
         write(*,*) 'c13          : ',c13
         write(*,*) 'rate         : ',D_C(3,k)
         write(*,*) 'rate (c13)   : ',D_C13(3,k)
         write(*,*) 'flux         : ',fluxes(3,k)
         write(*,*)
      endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

