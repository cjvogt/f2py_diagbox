c =====================================================================
c
c     SUBROUTINE WRITE_RATEANOMALIES_TOFILE
c
c     PURPOSE: writes the flux anomalies to a file
c     
c
c     VARIABLES:
c 
c     IN:  unit            : unit number to write to
c          tstart          : start time of integration
c          ts              : time step
c          nstep           : number of steps
c          D_C_avg         : average timeseries of dC/dt
c          D_C_anom        : anomaly timeseries of dC/dt
c          D_C_runavg      : running average timeseries of dC/dt
c
c
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     08.06.00   ng     first implementation
c     29.08.01   ng     added running averge
c 
c =====================================================================
c
      subroutine write_rateanomalies_tofile(unit,tstart,ts,nstep,
     $        D_C_avg,D_C_anom,D_C_runavg,
     $        error)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "anomalies.h"
# include "runavg.h"
c
      logical error
      integer unit
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      integer k,i
      double precision dcdt_avg(nproc_rates)
      double precision dcdt_anom(nproc_rates)
      double precision dcdt_runavg(nproc_rates)
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      error = .false.
      time = tstart
c
c --------------------------------------------------------------------
c     loop over all simulated values
c --------------------------------------------------------------------
c     
      do k = 1,nstep
c
c --------------------------------------------------------------------
c     change units of D_C (now it is umol kg-1) to umol kg-1 yr-1
c --------------------------------------------------------------------
c
         do i = 1,nproc_rates
            dcdt_avg(i) = d_C_avg(i,k) / ts
            dcdt_anom(i) = d_C_anom(i,k) / ts
            dcdt_runavg(i) = d_C_runavg(i,k) / ts
         enddo
c
c --------------------------------------------------------------------
c     write to file
c --------------------------------------------------------------------
c
         write(unit,9000) 
     $        time,
     $        (dcdt_avg(i),i=1,nproc_rates),
     $        (dcdt_anom(i),i=1,nproc_rates),
     $        (dcdt_runavg(i),i=1,nproc_rates)
c
 9000    format(f8.3,1x,21(f8.3,1x))
c
c --------------------------------------------------------------------
c     increase time
c --------------------------------------------------------------------
c
         time = time + ts
c
c --------------------------------------------------------------------
c     end of loop over all simulated values
c --------------------------------------------------------------------
c
      enddo
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c     
      return
c
c --------------------------------------------------------------------
c     error handling
c --------------------------------------------------------------------
c     
 9900 write(unit,*) '--WRITE_RATEANOMALIES_TOFILE: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











