c =====================================================================
c
c     SUBROUTINE compute_anomalies
c
c     PURPOSE: computes the anomalies of the rates and fluxes
c     
c
c     VARIABLES:
c 
c     IN:  tstart          : start time of integration
c          tend            : end time of integration
c          ts              : time step
c          nstep           : number of steps
c          D_C             : rates of change
c          fluxes          : fluxes
c
c
c
c     OUT: D_C_avg         : average timeseries of dC/dt
c          D_C_anom        : anomaly timeseries of dC/dt
c          fluxes_avg      : average timeseries of fluxes
c          fluxes_anom     : anomaly timeseries of fluxes
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     08.06.00   ng     first implementation
c 
c =====================================================================
c
      subroutine compute_anomalies(tstart,tend,ts,nstep,
     $        D_C,fluxes,
     $        D_C_avg,D_C_anom,fluxes_avg,fluxes_anom)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
# include "anomalies.h"
c
      logical error
      integer unit
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      integer k,i,k2
      integer nstep_per_year
      integer nsum(nstepmax)
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     determine the number of timesteps per year
c --------------------------------------------------------------------
c     
      nstep_per_year = int(1.d0/ts)
c
      write(*,*) '--COMPUTE_ANOMALIES: number of ts per year : ',
     $     nstep_per_year
c
c --------------------------------------------------------------------
c     intitialize all averages to zero
c --------------------------------------------------------------------
c     
      do k = 1,nstep
         do i = 1,nproc_rates
            D_C_avg(i,k) = 0.0d0
         enddo
c
         do i = 1,nproc_fluxes
            fluxes_avg(i,k) = 0.0d0
         enddo
c
         nsum(k) = 0
      enddo
c
c --------------------------------------------------------------------
c     sum values, store them in the first year,
c       keep track of number of samples
c --------------------------------------------------------------------
c     
      do k = 1,nstep
         k2 = mod(k-1,nstep_per_year)+1
         nsum(k2) = nsum(k2) + 1
c
         do i = 1,nproc_rates
            D_C_avg(i,k2) = D_C_avg(i,k2) + D_C(i,k)
         enddo
         do i = 1,nproc_fluxes
            fluxes_avg(i,k2) = fluxes_avg(i,k2) + fluxes(i,k)
         enddo
      enddo
c
c --------------------------------------------------------------------
c     compute averages
c --------------------------------------------------------------------
c
      do k = 1,nstep
         k2 = mod(k-1,nstep_per_year)+1
c
         if (k .eq. k2) then
            do i = 1,nproc_rates
               D_C_avg(i,k) = D_C_avg(i,k2)/nsum(k2)
            enddo
c     
            do i = 1,nproc_fluxes
               fluxes_avg(i,k) = fluxes_avg(i,k2)/nsum(k2)
            enddo
         else
            do i = 1,nproc_rates
               D_C_avg(i,k) = D_C_avg(i,k2)
            enddo
c     
            do i = 1,nproc_fluxes
               fluxes_avg(i,k) = fluxes_avg(i,k2)
            enddo
c            
         endif
      enddo
c
c --------------------------------------------------------------------
c     compute anomalies
c --------------------------------------------------------------------
c
      do k = 1,nstep
         do i = 1,nproc_rates
            D_C_anom(i,k) = D_C(i,k) - D_C_avg(i,k)
         enddo
         do i = 1,nproc_fluxes
            fluxes_anom(i,k) = fluxes(i,k) - fluxes_avg(i,k)
         enddo
      enddo
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c
      write(*,*) '--COMPUTE_ANOMALIES: terminated.'
c    
      return
c
c --------------------------------------------------------------------
c     error handling
c --------------------------------------------------------------------
c     
 9900 write(*,*) '--COMPUTE_ANOMALIES: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











