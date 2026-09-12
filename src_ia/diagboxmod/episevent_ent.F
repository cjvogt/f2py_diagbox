c =====================================================================
c
c     SUBROUTINE EPISEVENT_ENT
c 
c     PURPOSE: this routine calculates the contribution of entrainment
c                  on the mixed layer budget. This routine is the
c                  version where the episodic event scheme of C.D.
c                  Keeling is used. It is the version used for the
c                  station 's' paper.
c                  
c
c     VARIABLES:
c 
c     IN : k              : time loop counter
c          time           : actual time
c          ts             : time step
c          dc13_coeff      : character specifying which dC13 coefficients
c                              are used       
c                            <o1> : original data with 1st order fit   
c                            <o2> : original data with 2nd order fit   
c                            <a2> : artficial data with 2nd order fit   
c                            <n1> : adjusted data with 1st order fit
c                            <n2> : adjusted data with 2nd order fit   
c                            <m2> : adj. manip. data with 2nd o. fit  
c          h_coeff         : character specifying which h coeff are taken
c          dt_ent_re      : entrainment recurrence interval
c          const_vertgrad : logical indicating use of constant vert grad
c          dCdz_const     : value of constant vertical C gradient
c          dCdz_var       : value of variable vertical C gradient
c          D_dcdz         : correction for dc/dz
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
c      1.08.96   ng    first implementation
c      9.09.96   ng    use accurate formula for dc13_dz
c     18.10.99   ng    adatped for interannual variability runs
c
c =====================================================================
c
      subroutine episevent_ent(k,time,ts,
     $     dc13_coeff,h_coeff,dt_ent_re,
     $     const_vertgrad,dCdz_const,dCdz_var,D_dcdz,
     $     ddC13_dC_bml,ent_fact,
     $     h_new,h_old,
     $     c,c13,h,
     $     D_C,D_c13,fluxes,error)
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
# include "harm_coeff.h"
# include "diffent_params.h"
# include "c13_params.h"
      integer k
      logical error
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision rates_to_fluxes,delta_to_c13 ! functions
      double precision c13_to_delta
      double precision harmonic,spline
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision dc13_o,ddc13dz,dcdz,dc13dz
      double precision t_pre,t_pst,h_pre,h_pst
      double precision c_pre,c_pst,c13_pre,c13_pst,dc13_tmp
      double precision dh_ent,c_lb,c13_lb
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
c     determine the mixed layer depth at time t_pst and t_pre
c --------------------------------------------------------------------
c
      t_pre = time - 0.5d0*dt_ent_re + 0.5d0*ts
      t_pst = time + 0.5d0*dt_ent_re + 0.5d0*ts
c
      h_pre = harmonic(t_pre,c_h,7) + spline(t_pre,'h',error)
      h_pst = harmonic(t_pst,c_h,7) + spline(t_pre,'h',error)
         
      if ((h_new .ge. h_old) .and. (h_pst .gt. h_pre)) then
c
c --------------------------------------------------------------------
c     calculate contribution only during depeening of the 
c       mixed layer
c --------------------------------------------------------------------
c  
c
c --------------------------------------------------------------------
c     determine the C and c13 at time t_pre
c --------------------------------------------------------------------
c
         C_pre = harmonic(t_pre,c_sdic,5) 
     $        + spline(t_pre,'sdic',error)
c
         dC13_tmp = harmonic(t_pre,c_dc13,5) 
     $        + spline(t_pre,'dc13',error)
c
         c13_pre = delta_to_c13(dc13_tmp,c_pre)
         
c
c --------------------------------------------------------------------
c     determine the difference in mixed layer depth during the 
c       entrainment recurrence period, and the average concentration
c       at this depth
c --------------------------------------------------------------------
c     
        dh_ent = h_pst - h_pre
        C_lb    = C_pre + 0.5d0 * dCdz * dh_ent
        C13_lb  = C13_pre + 0.5d0 * dC13dz * dh_ent
c
c --------------------------------------------------------------------
c     the mixed layer concentration after the mixing
c --------------------------------------------------------------------
c
        c_pst = (c_pre*h_pre + c_lb*dh_ent)/h_pst
        c13_pst = (c13_pre*h_pre + c13_lb*dh_ent)/h_pst
c
c --------------------------------------------------------------------
c     calculate the average rate of change 
c --------------------------------------------------------------------
c
        D_c(3,k) = ent_fact * (C_pst - C_pre) * (ts/dt_ent_re)
        fluxes(3,k) = rates_to_fluxes(D_c(3,k),ts,h)
        D_c13(3,k) = ent_fact * (c13_pst - c13_pre) * (ts/dt_ent_re)
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
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

