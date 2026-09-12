c =====================================================================
c
c     DOUBLE PRECISION FUNCTION RATES_TO_FLUXES
c 
c     PURPOSE: this function calculates the fluxes for a given
c                 rate of change. It also converts the fluxes
c                 directly to a mol per square meter and second base.
c
c           Flux = dDIC * h / (ts * secPday * dayPyear)
c
c     VARIABLES:
c 
c     IN : rate            : rate of change of DIC [umol kg-1]
c          ts              : timestep [years]
c          h               : mixed layer depth
c     OUT: rates_to_fluxes : flux due to rate [mol m-2 s-1]
c 
c     REVISION:
c
c     date     author  remarks
c   
c     31.07.96   ng    first implementation
c     15.10.99   ng    adapted for interannual variability runs, units
c                        are now years
cc
c =====================================================================
c
      double precision function rates_to_fluxes(rate,ts,h)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
      double precision rate,ts,h
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision rho,sec_day,day_year
      parameter (rho = 1026.2d0,sec_day = 24.d0 * 3600.d0)
      parameter (day_year = 365.d0)
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      rates_to_fluxes = rate * 1.d-6 * rho * h 
     $                      / (ts * sec_day * day_year)
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

