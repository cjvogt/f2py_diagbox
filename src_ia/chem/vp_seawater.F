c =====================================================================
c
c     DOUBLE PRECISION FUNCTION VP_SEAWATER
c
c     PURPOSE: this function calculates the vapor pressure of
c                seawater. The vapor pressure of water in sea water
c                is related to that of pure water
c
c                vp(H20,sw) = VP(H20)*exp(-0.018*phi*sum(mb/m0))
c
c     REFERENCES: DOE (1994). Handbook of methods for the analysis
c                    of the various parameters of the carbon dioxide
c                    system in sea water; version 2, A.G. Dickson and
c                    C. Goyet, eds. ORNL/CDIAC-74. pages 5.22ff
c
c     SUBROUTINES: needs function vp_purewater
c
c     CHECK_VALUE: for 25 deg C, vp(H20,sw) = 3.1094 kPa
c
c     VARIABLES:
c 
c     IN:   temp        : temperature [deg C]
c           sal         : sallinity [psu]
c     OUT:  vp_seqwater : vapor pressure of pure-water [kPa]
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     1.08.96   ng      first implementation 
c 
c =====================================================================
c
      double precision function vp_seawater(temp,sal)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none

      double precision temp,sal
c
c --------------------------------------------------------------------
c     external functions
c --------------------------------------------------------------------
c   
      double precision vp_purewater
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      double precision phi,sum_mb_m0,vp_water
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     calculate vapor pressure for pure water
c --------------------------------------------------------------------
c        
      vp_water = vp_purewater(temp)
c
c --------------------------------------------------------------------
c     evaluate sum(mb/m0)
c --------------------------------------------------------------------
c   
      sum_mb_m0 = (31.998d0*sal)/(1.d3 - 1.005d0*sal)
c
c --------------------------------------------------------------------
c     evaluate osmotic coefficient phi
c --------------------------------------------------------------------
c   
      phi = 0.90799d0 - 0.08992*(0.5d0*sum_mb_m0) 
     $     + 0.18458d0*(0.5d0*sum_mb_m0)**2
     $     - 0.07395d0*(0.5d0*sum_mb_m0)**3
     $     - 0.00221d0*(0.5d0*sum_mb_m0)**4
c
c --------------------------------------------------------------------
c     calculate water vapor pressure
c --------------------------------------------------------------------
c        
      vp_seawater = vp_water*dexp(-0.018d0*phi*sum_mb_m0)
c
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      end
