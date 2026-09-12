c =====================================================================
c
c     SUBROUTINE GASEXCHANGE
c 
c     PURPOSE: this routine calculates the contribution of air-sea
c                gasexchange to the mixed layer budget
c
c     VARIABLES:
c 
c     IN : k            : time loop counter
c          ts           : time step
c          pv_relship   : piston velocity relationship
c          D_fco2_corr  : correction term for fugacity
c          gasex_fact   : multiplication factor for air-sea exchange
c
c          C            : obs C in the mixed layer
c          C13          : obs C13 in the mixed layer
c          h            : obs mixed layer depth
c          temp         : obs temperature
c          sal          : obs salinity
c          fco2_o       : obs fCO2 in the mixed layer
c          pco2_a       : obs xco2 in the atmosphere
c          dC13_a       : obs dC13 in the atmosphere
c          ws           : obs wind speed
c          sum_kex      : sum of kex to calculate average
c
c     OUT: D_C          : modified array holding the rates of change for C
c          D_C13        : modified array holding the rates of change for C13
c          fluxes       : modified array holding the fluxes
c          sum_kex      : modified sum_kex 
c          fco2_a_obs   : observed atmospheric fco2 [uatm]
c          pco2_a_dry   : observed atmospheric pco2 [uatm]
c 
c     REVISION:
c
c     date     author  remarks
c   
c      1.08.96   ng    first implementation
c     22.08.96   ng    included estimates of relative humidity and
c                         total pressure from St.Davids head on
c                         the island of Bermuda
c     16.12.96   ng    found error in calculating fugacity, permit
c                         to use old,wrong formulation with a cpp 
c                         option (-Dwith_error_fugacity)
c     28.05.97   ng    found error in wanninkhof wind speed dependency of
c                         pv (**(-0.5) instead of **(0.5))
c     15.10.99   ng    adapted for interannual variability run, make
c                         sure units are compatible (now using years)
c     26.10.02   ng    added pco2_a_obs_dry
c
c =====================================================================
c
      subroutine gasexchange(k,ts,
     $     pv_relship,D_fco2_corr,gasex_fact,
     $     C,C13,h,temp,sal,fco2_o,pco2_a,dC13_a,ws,
     $     sum_kex,fco2_a_obs,pco2_a_dry,
     $     D_C,D_C13,fluxes)
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
# include "gasex_params.h"
# include "averages.h"
# include "c13_params.h"
# include "obs_values.h"
      integer k
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision fluxes_to_rates,delta_to_ratio ! functions
      double precision solco2m3, schmidt_number, vp_seawater ! functions
      double precision c13_to_delta,fco2 ! function
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision sc,sc_old,scnfact,sol,pv, tk,kex,vp_sw
      double precision pco2_a_vp,fco2_a,fco2_a2
      double precision alpha_eq,alpha_am
      double precision flux_C13_ex
      double precision dc13_o
      double precision dc13_ex,rstar,p_dry
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      tk = temp+273.15d0
c
c --------------------------------------------------------------------
c     calculate Schmidt-number (using Wanninkhof, 1992)
c         no units
c     sc_old is the old relationship from Heimann, pers. comm.
c
c --------------------------------------------------------------------
c
#ifdef same_as_vs13    
      sc = 600.d0*10.d0**(1966.d0/293.15d0*(293.15d0/tk-1.0d0))
      Scnfact = sc/600.d0
#else
      sc = schmidt_number('CO2',temp)
      ScnFact = sc/660.d0
#endif
c
c      write(*,*) 'sc (norm to 660) and sc_old (norm to 600) : ',
c     $     sc/660.d0,sc_old/600.d0
c
c --------------------------------------------------------------------
c     calculate CO2 solubility (using Weiss, 1974)
c       units are mol/(m3*atm) then change directly to units of
c       mol m-3 uatm-1
c --------------------------------------------------------------------
c        
      sol = solco2m3(temp,sal)
      sol = sol * 1.d-6
c
c --------------------------------------------------------------------
c     calculate gas transfer velocity (also called piston velocity) pv
c        for Wanninkhof use formula for long term average winds
c        his formula (1)
c     the relationships are given for wind-speed in m s-1 and
c        resulting pv in cm hr-1
c --------------------------------------------------------------------
c        
      if (pv_relship .eq. 'lm') then
c
         pv = 0.17d0*(ScnFact)**(-2.0d0/3.0d0) * ws + 
     $        2.68d0*(ws-3.6d0)*(ScNFact)**(-0.5d0) ! cm hr-1
c
      elseif (pv_relship .eq. 'wa') then
c
         pv = 0.39d0 * ws**2.d0 / dsqrt(ScnFact) ! cm hr-1
c
      endif
c
      pv = pv / (100.d0 * 3600.d0) ! convert to m s-1
c
c
c --------------------------------------------------------------------
c     calculate gas transfer coefficient in units of 
c         mol m-2 s-1 uatm-1
c     include directly also the air-sea gasexchange factor
c --------------------------------------------------------------------
c        
      kex = gasex_fact * pv * Sol
c
      Sum_kex = Sum_kex + kex
c
c      write(*,*) 'sol, pv, kex : ',sol,pv,kex
c
c --------------------------------------------------------------------
c     calculate atmospheric fugacity
c        calculate first water vapour pressure at the temperature
c        of the water and assuming 100% saturation
c        then calculate pco2_a_vp which includes effect of 
c        vapor pressure of seawater, assume 1.03 atm total pressure
c --------------------------------------------------------------------
c
#ifdef same_as_vs13
      fco2_a = pco2_a
      fco2_a_obs(k) = fco2_a
#else
      vp_sw = vp_seawater(temp,sal) ! units [kPa]
      vp_sw = vp_sw / 101.325d0 ! units [atm]
c
      p_dry = mean_air_press - mean_rel_humid*vp_sw
c
c      write(*,*) 'k, p_dry : ',k,p_dry
c
      fco2_a = fco2(temp,p_dry,pco2_a) ! units [uatm]
      pco2_a_dry(k) = pco2_a*p_dry
      fco2_a_obs(k) = fco2_a

c
c --------------------------------------------------------------------
c     old stuff
c --------------------------------------------------------------------
c  
c
c#  ifdef with_error_fugacity
c      pco2_a_vp = pco2_a * p_dry ! units [ppm*atm = uatm]
c      fco2_a = fco2(temp,1.d0,pco2_a_vp) ! units [uatm]
c#  else
c
c      fco2_a = fco2(temp,p_dry,pco2_a) ! units [uatm]
c#  endif
c
c      write(*,*) 'pco2_a,pco2_a_vp,fco2_a,fco2_a2 : ',
c     $     pco2_a,pco2_a_vp,fco2_a,fco2_a2
#endif
c
c --------------------------------------------------------------------
c     apply correction for fCO2_o
c --------------------------------------------------------------------
c  
      fCO2_o = fCO2_o + D_fco2_corr
c
c --------------------------------------------------------------------
c     calculate finally the flux of CO2 and rate of change
c --------------------------------------------------------------------
c        
      fluxes(1,k) = kex * (fCO2_a - fCO2_o) ! units [mol/m2/s]
c
      D_C(1,k) = fluxes_to_rates(fluxes(1,k), ts, h) ! units [umol/kg]
c
c --------------------------------------------------------------------
c     for debugging purposes
c --------------------------------------------------------------------
c        
c      write(*,*) 'k : ',k
c      write(*,*) 'kex, fCO2_a, fCO2_o : ',kex,fco2_a,fco2_o
c      write(*,*) 'flux, D_C : ',fluxes(1,k),D_C(1,k)
c
c --------------------------------------------------------------------
c     calculate the equilibrium fractionation factor for dC13
c       from Heimann and Keeling, 1989
c     define also the kinetic fractionation factor alpha_am
c       also from Heimann and Keeling, 1989
c --------------------------------------------------------------------
c        
      Alpha_eq = 1.02389d0 - 9.483d0/TK
      Alpha_am = 0.99820d0
c
c --------------------------------------------------------------------
c     calculate the air-sea exchange flux for C13
c       units are [mol 13C m-2 s-1] and the rate of change
c --------------------------------------------------------------------
c   
      dc13_o = c13_to_delta(c13,c) ! calculate the dC13 value in the ocean
c
      flux_c13_ex = 
     $     kex * alpha_am *(rstar_std/r_std) 
     $     * (fCO2_a * delta_to_ratio(dC13_a)
     $     - fCO2_o * delta_to_ratio(dc13_o) * alpha_eq) ! units [mol/m2/s]
c
      D_c13(1,k) = fluxes_to_rates(flux_c13_ex, ts, h) ! units [umol/kg]
c
c --------------------------------------------------------------------
c     for debugging purposes
c --------------------------------------------------------------------
c
#ifdef debugging
      if (mod(k,90) .eq. 0) then
         write(*,*)
         write(*,*) '--GASEXCHANGE: '
         write(*,*) '  k : ',k
         write(*,*) '  alpha_am : ',alpha_am
         write(*,*) '  alpha_eq : ',alpha_eq
         write(*,*) '  flux_c13_ex : ',flux_c13_ex
         write(*,*) '  D_c13_ex    : ',D_c13(1,k)
c
         if (D_c(1,k) .ne. 0.d0) then
            rstar = D_c13(1,k) / D_c(1,k)
            dc13_ex = 1.d3 * (1.d0/r_std*(rstar/(1.d0-rstar)) - 1.d0)
         else
            dc13_ex = 0.d0
         endif
         write(*,*) '  dc13_ex : ',dc13_ex
      endif
#endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

