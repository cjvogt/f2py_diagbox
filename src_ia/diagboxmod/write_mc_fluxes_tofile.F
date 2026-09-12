c =====================================================================
c
c     SUBROUTINE WRITE_MC_FLUXES_TOFILE
c
c     PURPOSE: writes the rates of change to a file specified by
c                the unit number unit, this special version writes also
c                the means, weighted means, the sigma and weighted sigmas
c                calculated from the monte carlo analysis to the unit 
c     
c
c     VARIABLES:
c 
c     IN:  unit            : unit number to write to
c          tstart          : start time of integration
c          ts              : time step
c          nstep2          : number of samples in sg* arrays
c          rsmpl_step      : resampling time step of sg* arrays
c          fluxes_std      : fluxes of the different processes
c          mean_fluxes     : mean fluxes of mc runs
c          wmean_fluxes    : weighted mean fluxes of mc runs
c          sg_fluxes       : 1-sg uncertainty of fluxes
c          wsg_fluxes      : weighted 1-sg uncertainty of fluxes
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     03.08.96   ng     first implementation
c     18.10.99   ng     adapted for interannual variability
c 
c =====================================================================
c
      subroutine write_mc_fluxes_tofile(unit,tstart,ts,
     $     nstep2,rsmpl_step,
     $     fluxes_std,
     $     mean_fluxes,wmean_fluxes,sg_fluxes,wsg_fluxes,
     $     error)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "simulation.h"
# include "monte_carlo.h"
# include "periods.h"
# include "processes.h"
# include "mc_statistics.h"
# include "slabs.h"
# include "mc_plot_slabs.h"
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
      parameter (day_year = 365.)
      character*27 tmp(nproc_fluxes)

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
      do k2 = 1,nstep2
c
         k = (k2-1)*rsmpl_step + 1
c
c --------------------------------------------------------------------
c     change units of fluxes to mol m-2 day-1
c      calculate also upper and lower boundaries from mc analysis
c --------------------------------------------------------------------
c
         do i = 1,nproc_fluxes
            fluxes_gd(i) = fluxes_std(i,k) * sec_day * day_year 
            fluxes_gd_ub(i) = (fluxes_std(i,k) + wsg_fluxes(i,k2))
     $           * sec_day * day_year 
            fluxes_gd_lb(i) = (fluxes_std(i,k) - wsg_fluxes(i,k2))
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
c
c --------------------------------------------------------------------
c     end of loop over all simulated values
c --------------------------------------------------------------------
c
c
      enddo
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
 9900 write(unit,*) '--WRITE_FLUXES_TOFILE: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











