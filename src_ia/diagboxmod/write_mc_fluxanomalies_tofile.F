c =====================================================================
c
c     SUBROUTINE WRITE_MC_FLUXANOMALIES_TOFILE
c
c     PURPOSE: writes the flux anomalies to a file
c     
c
c     VARIABLES:
c 
c     IN:  unit                : unit number to write to
c          tstart              : start time of integration
c          ts                  : time step
c          nstep               : number of steps
c          nstep2              : number of samples in sg* arrays
c          rsmpl_step          : resampling time step of sg* arrays
c          fluxes_std_runavg   : running avg of fluxes_std
c          mean_fluxes_runavg  : running avg of mean_fluxes
c          sg_fluxes_runavg    : running avg of sg_fluxes
c
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     09.06.02   ng     first implementation, on the basis of write_mc_fluxes
c                                         to_file.F
c 
c =====================================================================
c
      subroutine write_mc_fluxanomalies_tofile(
     $     unit,tstart,ts,nstep,nstep2,rsmpl_step,
     $     fluxes_std_runavg, mean_fluxes_runavg,
     $     sg_fluxes_runavg,error)
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
      integer k,k2,i
      double precision fluxes_gd(nproc_fluxes)
      double precision fluxes_gd_ub(nproc_fluxes)
      double precision fluxes_gd_lb(nproc_fluxes)
      double precision sec_day,g_mol,day_year
      parameter (sec_day = 24.d0 * 3600.d0, g_mol = 12.01d0)
      parameter (day_year = 365.d0)
      character*27 tmp(nproc_fluxes)
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     statement
c --------------------------------------------------------------------
c     
      write(*,*) '--WRITE_MC_FLUXANOMALIES_TOFILE: saving values...'


      error = .false.
      time = tstart
c
c --------------------------------------------------------------------
c     loop over all simulated values
c --------------------------------------------------------------------
c     
      do k2 = 1,nstep2
c
         k = (k2-1)*rsmpl_step + 1
         
         write(*,*) ' working on step k2,k : ',k2,k
c
c --------------------------------------------------------------------
c     change units of fluxes to mol m-2 day-1
c      calculate also upper and lower boundaries from mc analysis
c --------------------------------------------------------------------
c
         do i = 1,nproc_fluxes
            fluxes_gd(i) = fluxes_std_runavg(i,k) * sec_day * day_year 
            fluxes_gd_ub(i) = (fluxes_std_runavg(i,k) 
     $           + sg_fluxes_runavg(i,k2))
     $           * sec_day * day_year 
            fluxes_gd_lb(i) = (fluxes_std_runavg(i,k) 
     $           - sg_fluxes_runavg(i,k2))
     $           * sec_day * day_year 
         enddo
c
c --------------------------------------------------------------------
c     write to temporary array for output
c --------------------------------------------------------------------
c
         do i = 1,nproc_fluxes-3
            write(tmp(i),'(3(f8.4,1x))') 
     $           fluxes_gd(i),fluxes_gd_lb(i),fluxes_gd_ub(i)
         enddo
c
c         write(*,*) 'values : ',(tmp(i),i=1,nproc_fluxes-3)
c
c --------------------------------------------------------------------
c     write to file
c --------------------------------------------------------------------
c
         write(unit,9000) 
     $        time,
     $        (tmp(i),i=1,nproc_fluxes-3)
c
 9000    format(f9.4,1x,10(a26))
c
c --------------------------------------------------------------------
c     increase time
c --------------------------------------------------------------------
c
         time = time + float(rsmpl_step)*ts
         write(*,*) ' time : ',time
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











