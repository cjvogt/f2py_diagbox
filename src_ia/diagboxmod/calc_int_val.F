c =====================================================================
c
c     SUBROUTINE CALC_INT_VAL
c 
c     PURPOSE: this routine calculates the integrated rates of change
c              and fluxes over the different periods
c
c     VARIABLES:
c
c     IN : D_c        : rates of change of C
c          D_c13      : rates of change of C13
c          fluxes     : fluxes of C
c          tstart     : start time of integration
c          tend       : end time of integration
c          ts         : time step
c          nstep      : number of steps
c 
c     OUT: Int_Fluxes : two dim array holding the integrated Fluxes of
c                          the different processes during the
c                          different periods
c          Int_Rates  : two dim array holding the integrated Rates of
c                          change of the different processes during the
c                          different periods
c          nyears     : number of analysis years
c
c     REVISION:
c
c     date     author  remarks
c   
c      3.08.96   ng    first implementation
c     17.10.99   ng    adapting for interannual var. runs, changed time
c                        time from days to years.
cc
c =====================================================================
c
      subroutine calc_int_val(D_c,D_c13,fluxes,tstart,tend,ts,nstep,
     $     int_fluxes,int_rates,nyears)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "simulation.h"
# include "periods.h"
# include "processes.h"
# include "slabs.h"
# include "intval.h"
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer ny,i,k,iyr

      double precision rho,sec_day,day_year
      parameter (rho = 1026.2d0,sec_day = 24.d0 * 3600.d0)
      parameter (day_year = 365.d0)
c
      integer nstep_yr(nyearmax)
c
      double precision eps
      parameter(eps = 1.d-8)

c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--CALC_INT_VAL : calculating integrated values...'
c
c --------------------------------------------------------------------
c     make sure that all integrated values are set to zero
c      at the beginning
c --------------------------------------------------------------------
c        
      do ny = 1,nyearmax
         do i = 1,nproc_fluxes
            int_fluxes(ny,i) = 0.0d0
         enddo
         do i = 1,nproc_rates
            int_rates(ny,i) = 0.0d0
         enddo
      enddo
c
c --------------------------------------------------------------------
c     determine the number of analysis years 
c --------------------------------------------------------------------
c        
      nyears = int(tend - tstart)
c
c --------------------------------------------------------------------
c     set the numbers of time steps in an analysis year to zero
c --------------------------------------------------------------------
c        
      do ny = 1,nyearmax
         nstep_yr(ny) = 0
      enddo
c
c --------------------------------------------------------------------
c     loop over all timesteps
c --------------------------------------------------------------------
c        
      time = tstart
c
      do k = 1,nstep
c
c --------------------------------------------------------------------
c     determine analysis year
c      add a little eps (small number) because of inexact addition
c      operation
c      keep track of how many time steps there are in an analysis year
c --------------------------------------------------------------------
c        
         iyr = int(time + eps - tstart - year_start) + 1
         nstep_yr(iyr) = nstep_yr(iyr) + 1
c         if (abs(time - nint(time)) .le. 0.0025) then
c            write(*,*) '--CALC_INT_VAL, time = ',time
c            write(*,*) '  iyr,nstep_yr(iyr) = ',iyr,nstep_yr(iyr)
c         endif
c
c --------------------------------------------------------------------
c     calculate integrated rates of change (no unit conversion yet)
c       here (umol/kg)
c     no need to multiply with ts, because D_C denotes just Delta_C
c       and not Delta_C/Delta_t
c --------------------------------------------------------------------
c        
         do i = 1,nproc_rates
            int_rates(iyr,i) = int_rates(iyr,i) + D_C(i,k)
         enddo
c
c --------------------------------------------------------------------
c     calculate integrated fluxes
c     here ts must be included to get integral, because fluxes
c      are in mol m-2 s-1
c --------------------------------------------------------------------
c        
         do i = 1,nproc_fluxes
            int_fluxes(iyr,i) = int_fluxes(iyr,i) 
     $           + fluxes(i,k) * ts * sec_day * day_year
         enddo
c
c --------------------------------------------------------------------
c     increase time
c --------------------------------------------------------------------
c        
         time = time + ts
c
c --------------------------------------------------------------------
c     end of loop over all timesteps
c --------------------------------------------------------------------
c        
      enddo
c
c --------------------------------------------------------------------
c     give info
c --------------------------------------------------------------------
c        
c      write(*,*)
c      write(*,*) '--CALC_INT_VAL: number of analysis years : ',nyears
c      write(*,*) '        number of steps in each of these years: '
c      do ny = 1,nyears
c         write(*,*) '            year = ',ny, '  # steps = ',
c     $        nstep_yr(ny)
c      enddo
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

