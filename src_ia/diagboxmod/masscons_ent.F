c =====================================================================
c
c     SUBROUTINE MASSCONS_ENT
c 
c     PURPOSE: this routine calculates the contribution of entrainment
c                  on the mixed layer budget. This routine is the
c                  version where the mass conservation scheme of 
c                  T. Stocker is used. Basic idea is mass conservation 
c                  in the combined system mixed layer and submixed 
c                  layer down to h_th. 
c
c     VARIABLES:
c 
c     IN : k              : time loop counter
c          ts             : time step
c          h_th           : depth at which thermocline values for 
c                                 mc entrainment are taken
c          c_th           : C concentration at h_th [umol kg-1]
c          dC13_th        : dC13 concentration at h_th [per mil]
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
c      1.08.96   ng    first implementation
c
c =====================================================================
c
      subroutine masscons_ent(k,ts,
     $     h_th,c_th,dc13_th,ent_fact,
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
      double precision rates_to_fluxes,delta_to_c13 ! functions
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision c13_th
      double precision dc13_ent,rstar
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
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
         D_C(3,k) = ent_fact * (c_th - c)*(h_new-h_old)/
     $        (h_old + h_th)
         fluxes(3,k) = rates_to_fluxes(D_C(3,k),ts,h)
c
c --------------------------------------------------------------------
c     calculate the corresponding change for C13
c --------------------------------------------------------------------
c
         c13_th = delta_to_c13(dc13_th,c_th)
c
         D_c13(3,k) = ent_fact * (c13_th - c13)*(h_new - h_old)/
     $        (h_old + h_th)
c
c --------------------------------------------------------------------
c     shoaling period
c --------------------------------------------------------------------
c
      else
         D_C(3,k) = 0.d0
         fluxes(3,k) = 0.d0
         D_c13(3,k) = 0.d0
      endif
c
c --------------------------------------------------------------------
c     for debugging purposes
c --------------------------------------------------------------------
c
#ifdef debugging
      if (mod(k,90) .eq. 0) then
         write(*,*) 'k : ',k
         write(*,*) 'D_c13_ent      : ',D_c13(3,k)
c     
         if (D_c(3,k) .ne. 0.d0) then
            rstar = D_c13(3,k) / D_c(3,k)
            dc13_ent = 1.d3 * (1.d0/r_std*(rstar/(1.d0-rstar)) - 1.d0)
         else
            dc13_ent = 0.d0
         endif
         write(*,*) 'dc13_ent : ',dc13_ent
      endif
#endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

