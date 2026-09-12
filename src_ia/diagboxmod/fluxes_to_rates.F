c =====================================================================
c
c     DOUBLE PRECISION FUNCTION FLUXES_TO_RATES
c 
c     PURPOSE: this function calculates the rates for a given
c                 flux. 
c 
c           Flux = dDIC * h / (ts * secPday * daysPyear)
c
c     VARIABLES:
c 
c     IN : fluxes          : fluxes [mol m-2 s-1]
c          ts              : timestep [years]
c          h               : mixed layer depth
c     OUT: fluxes_to_rates : rate due to fluxes [umol kg-1]
c 
c     REVISION:
c
c     date     author  remarks
c   
c     31.07.96   ng    first implementation
c     15.10.99   ng    adapted for interannual variability runs, units
c                        are now years
c
c =====================================================================
c
      double precision function fluxes_to_rates(fluxes,ts,h)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
      double precision fluxes,ts,h
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision rho,sec_day, day_year
      parameter (rho = 1026.2d0,sec_day = 24.d0 * 3600.d0)
      parameter (day_year = 365.d0)
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      fluxes_to_rates = fluxes * (ts * sec_day * day_year) 
     $                 / (1.d-6 * rho * h)
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

