c =====================================================================
c
c     SUBROUTINE VERT_DIFF
c 
c     PURPOSE: this routine calculates the contribution of vertical
c                 diffusion on the mixed layer budget
c
c     VARIABLES:
c 
c     IN : k              : time loop counter
c          ts             : time step
c          ent_scheme     : entrainment scheme
c          h_th           : depth at which thermocline values for 
c                                 mc entrainment are taken
c          c_th           : C concentration at h_th [umol kg-1]
c          dC13_th        : dC13 concentration at h_th [per mil]
c          const_kz       : logical indicating use of constant Kz
c          kz_const       : value of constant Kz
c          kz_var         : value of seasonally varying kz
c          const_vertgrad : logical indicating use of constant vert grad
c          dCdz_const     : value of constant vertical C gradient
c          dCdz_var       : value of variable vertical C gradient
c          D_dcdz         : correction for vertical gradient
c          ddC13_dC_bml   : ratio of dd13C/dz to dC/dz below mixed layer
c          diff_fact      : vertical diffusion mulitplication factor
c
c          C              : obs C in the mixed layer
c          C13            : obs C13 in the mixed layer
c          h              : obs mixed layer depth
c          dCdz           : vertical C gradient below mixed layer
c
c     OUT: D_C            : modified array holding the rates of change for C
c          D_C13          : modified array holding the rates of change for C13
c          fluxes         : modified array holding the fluxes
c 
c     REVISION:
c
c     date     author  remarks
c   
c      1.08.96   ng    first implementation
c      9.09.96   ng    use accurate formula for dc13_dz
c     17.10.99   ng    adapted for interannual variability runs, time is
c                         now in years rather than days
c
c =====================================================================
c
      subroutine vert_diff(k,ts,
     $     ent_scheme,
     $     h_th,c_th,dc13_th,
     $     const_kz,kz_const,kz_var,
     $     const_vertgrad,dCdz_const,dCdz_var,D_dcdz,
     $     ddC13_dC_bml,diff_fact,
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
      double precision fluxes_to_rates,c13_to_delta ! functions
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision rho
      parameter (rho = 1026.2d0)
      double precision kz,dcdz,dC13dz,ddC13dz,dc13_o
      double precision flux_diff,flux_c13_diff
      double precision dc13_diff,rstar
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     evaluate value of kz [m2 s-1]
c --------------------------------------------------------------------
c        
      if (const_kz) then
         kz = kz_const
      else
         kz = kz_var
      endif
c
c --------------------------------------------------------------------
c     evaluate vertical C and dC13 gradient
c --------------------------------------------------------------------
c    
      dC13_o = c13_to_delta(c13,c)
c
      if (const_vertgrad) then
         dCdz = dCdz_const
         ddC13dz = dCdz * ddC13_dC_bml
      else
         if (ent_scheme .eq. 'mc') then
            dCdz = (C_th-C)/(h_th-h)
            ddC13dz = (dC13_th-dC13_o)/(h_th-h)
         else
            dCdz = dCdz_var + D_dCdz
            ddC13dz = dCdz * ddC13_dC_bml
         endif
      endif
c
c --------------------------------------------------------------------
c     calulate the vertical C13 gradient from the vertical
c       dC13 gradient (approximate formula)
c       units [umol kg-1 m-1]
c     changed to accurate formula
c --------------------------------------------------------------------
c   
c      dc13dz = r_std/(1+r_std) * 
c     $     (c * ddC13dz*1.d-3 + dcdz*(dc13_o*1.d-3 + 1.d0))
c
      dc13dz = r_std/(1+r_std*dc13_o*1.d-3+r_std) * 
     $     ((dc13_o*1.d-3 + 1.d0)*dcdz +
     $      c / (1+r_std*dc13_o*1.d-3+r_std) * ddc13dz * 1.d-3)
c
c      write(*,*) 'ddC13/dz, dc13/dz : ',ddc13dz,dc13dz
c
c --------------------------------------------------------------------
c     calculate total C flux using first fick law
c       convert to mol m-2 s-1 and include directly also
c       multiplication factor diff_fact
c --------------------------------------------------------------------
c
      flux_diff = Kz * dCdz   ! units [m s-1 umol kg-1]
      fluxes(2,k) = flux_diff * diff_fact * 1.d-6 * rho ! units [mol m-2 s-1]

      D_C(2,k) = fluxes_to_rates(fluxes(2,k),ts, h)
c
c --------------------------------------------------------------------
c     calculate total C13 flux using first fick law
c       convert to mol m-2 s-1 and include directly also
c       multiplication factor diff_fact
c --------------------------------------------------------------------
c        
      flux_c13_diff = Kz * dc13dz   ! units [m s-1 umol kg-1]
      flux_c13_diff = flux_c13_diff * diff_fact * 1.d-6 * rho 
                                ! units [mol m-2 s-1]
c
      D_c13(2,k) = fluxes_to_rates(flux_c13_diff,ts, h) ! units [umol/kg]
c
c --------------------------------------------------------------------
c     for debugging purposes
c --------------------------------------------------------------------
c
#ifdef debugging
      if (mod(k,90) .eq. 0) then
         write(*,*) 
         write(*,*) '--VERT_DIFFUSION :'
         write(*,*) '  k : ',k
         write(*,*) '  flux_c13_diff : ',flux_c13_diff
         write(*,*) '  D_c13_diff      : ',D_c13(2,k)
c     
         if (D_c(2,k) .ne. 0.d0) then
            rstar = D_c13(2,k) / D_c(2,k)
            dc13_diff = 1.d3 * (1.d0/r_std*(rstar/(1.d0-rstar)) - 1.d0)
         else
            dc13_diff = 0.d0
         endif
         write(*,*) '  dc13_diff : ',dc13_diff
      endif
#endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

