c =====================================================================
c
c     SUBROUTINE WRITE_MC_FLUXANOMALIES_TOFILE
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
c          fluxes_avg      : timeseries of average fluxes
c          fluxes_anom     : timeseries of flux anomalies
c          fluxes_runavg   : timeseries of running averaged fluxes
c
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     08.06.00   ng     first implementation
c     29.08.01   ng     added fluxes_runavg
c 
c =====================================================================
c
      subroutine write_mc_fluxanomalies_tofile(
     $     unit,tstart,ts,nstep,nstep2,
     $     fluxes_std_runavg, mean_fluxes_runavg,sg_fluxes_runavg,error)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "diagboxmod.h"
# include "simulation.h"
# include "monte_carlo.h"
# include "processes.h"
# include "anomalies.h"
# include "mc_plot_slabs.h"
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
      double precision fluxes_std_gd(nproc_fluxes)
      double precision fluxes_mean_gd(nproc_fluxes)
      double precision fluxes_sg_gd(nproc_fluxes)
       double precision sec_day,g_mol,day_year
      parameter (sec_day = 24.d0 * 3600.d0, g_mol = 12.01d0)
      parameter (day_year = 365.d0)
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
c     change units of fluxes to mol m-2 y-1
c --------------------------------------------------------------------
c
         do i = 1,nproc_fluxes
            fluxes_std_gd(i) = fluxes_std_runavg(i,k) * sec_day * 
     $           day_year 
            fluxes_mean_gd(i) = mean_fluxes_runavg(i,k) * sec_day * 
     $           day_year 
            fluxes_sg_gd(i) = sg_fluxes_runavg(i,k) * sec_day * 
     $           day_year 
         enddo
c
c --------------------------------------------------------------------
c     write to file
c --------------------------------------------------------------------
c
         write(unit,9000) 
     $        time,
     $        (fluxes_std_gd(i),i=1,nproc_fluxes),
     $        (fluxes_mean_gd(i),i=1,nproc_fluxes),
     $        (fluxes_sg_gd(i),i=1,nproc_fluxes)
c
 9000    format(f9.4,1x,3(7(f8.2,1x),3(f10.2,1x)))
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
 9900 write(unit,*) '--WRITE_MC_FLUXANOMALIES_TOFILE: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











