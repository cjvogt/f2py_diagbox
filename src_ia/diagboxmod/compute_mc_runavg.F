c =====================================================================
c
c     SUBROUTINE compute_mc_runavg
c
c     PURPOSE: computes the one year running averages of 
c                               the rates and fluxes from MC simulations
c     
c
c     VARIABLES:
c 
c     IN:  tstart          : start time of integration
c          tend            : end time of integration
c          ts              : time step
c          nstep           : number of steps
c          nstep2          : number of samples in sg* arrays
c          rsmpl_step      : resampling time step of sg* arrays
c          fluxes_std      : fluxes of the different processes
c          mean_fluxes     : mean fluxes of mc runs
c          sg_fluxes       : 1-sg uncertainty of fluxes
c          D_c_std         : rates of change of the different processes
c          mean_d_c        : mean d_c of mc runs
c          sg_d_c          : 1-sg uncertainty of d_c
c
c
c     OUT: D_C_std_runavg    : running avg of D_C standard
c          mean_d_c_runavg   : running avg of mean_D_C
c          sg_D_c_runavg     : running avg of sg_D_c
c          fluxes_std_runavg : running avg of fluxes_std
c          mean_fluxes_runavg: running avg of mean_fluxes
c          sg_fluxes_runavg  : running avg of sg_fluxes
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     08.06.02   ng     first implementation, on the basis of 
c                             compute_runavg.F
c 
c =====================================================================
c
      subroutine compute_mc_runavg(
     $     tstart,tend,ts,nstep, nstep2,
     $     rsmpl_step, D_c_std, mean_d_c,sg_D_c, fluxes_std,
     $     mean_fluxes,sg_fluxes, 
     $     D_C_std_runavg, mean_d_c_runavg,
     $     sg_D_c_runavg, fluxes_std_runavg, mean_fluxes_runavg,
     $     sg_fluxes_runavg,error)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "periods.h"
# include "monte_carlo.h"
# include "mc_statistics.h"
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
      integer k,i,k2,kk
      integer nstep_per_year,nstep_per_halfyr,nstep_start,nstep_stop
      double precision sum1,sum2
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c =====================================================================
c     need to do two passes, one for the stddev and mean, which
c       has been resamples, and one for the fluxes_std and D_C_std
c       which has the original time resolution
c =====================================================================
c
c
c --------------------------------------------------------------------
c     determine the number of timesteps per year, for the resampled
c       case
c --------------------------------------------------------------------
c     
      nstep_per_year = int(1.d0/ts)/rsmpl_step
c
      write(*,*) '--COMPUTE_MC_RUNAVG: number of ts per year',
     $     ' for resampled fields : ',
     $     nstep_per_year, ' resample timestep : ',rsmpl_step
c
c --------------------------------------------------------------------
c     determine the first and last steps
c --------------------------------------------------------------------
c     
      nstep_per_halfyr = int(nstep_per_year/2)
      nstep_start = nstep_per_halfyr+1
      nstep_stop = nstep2 - nstep_per_halfyr
      if (nstep_stop .le. nstep_start) then
         write(*,*) '--COMPUTE_RUNAVG: running average not possible,',
     $        ' because record is not long enough '
         write(*,*) '     nstep_start, nstep_stop : ',nstep_start,
     $        nstep_stop
         error = .true.
         goto 9900
      endif
c
      write(*,*) '--COMPUTE_MC_RUNAVG: nstep_start,nstep_stop : ',
     $     nstep_start,nstep_stop
c
c --------------------------------------------------------------------
c     pad at the beginning
c --------------------------------------------------------------------
c     
      
      do k = 1,nstep_start-1
         do i = 1,nproc_rates
            mean_d_c_runavg(i,k) = 0.0d0
            sg_D_c_runavg(i,k) = 0.0d0
         enddo
         do i = 1,nproc_fluxes
            mean_fluxes_runavg(i,k) = 0.0d0
            sg_fluxes_runavg(i,k) = 0.0d0
         enddo
      enddo
c
c --------------------------------------------------------------------
c     do the central (the real) part
c --------------------------------------------------------------------
c     
      do k = nstep_start,nstep_stop
c         write(*,*) '  working on step : ',k
         do i = 1,nproc_rates
            sum1 = 0.d0
            sum2 = 0.d0
             do kk = k-nstep_per_halfyr,k+nstep_per_halfyr
               sum1 = sum1 + mean_D_C(i,kk)
               sum2 = sum2 + sg_D_C(i,kk)
            enddo
            mean_D_C_runavg(i,k) = sum1/(2*nstep_per_halfyr+1)
            sg_D_C_runavg(i,k) = sum2/(2*nstep_per_halfyr+1)
         enddo
         do i= 1, nproc_fluxes
            sum1 = 0.0d0
            sum2 = 0.0d0
            do kk = k-nstep_per_halfyr,k+nstep_per_halfyr
               sum1 = sum1 + mean_fluxes(i,kk)
               sum2 = sum2 + sg_fluxes(i,kk)
            enddo
            mean_fluxes_runavg(i,k) = sum1/(2*nstep_per_halfyr+1)
            sg_fluxes_runavg(i,k) = sum2/(2*nstep_per_halfyr+1)
         enddo
      enddo
c
c --------------------------------------------------------------------
c     pad at the end
c --------------------------------------------------------------------
c     
      
      do k = nstep_stop+1,nstep2
         do i = 1,nproc_rates
            mean_d_c_runavg(i,k) = 0.0d0
            sg_D_c_runavg(i,k) = 0.0d0
         enddo
         do i = 1,nproc_fluxes
            mean_fluxes_runavg(i,k) = 0.0d0
            sg_fluxes_runavg(i,k) = 0.0d0
         enddo
      enddo
             
c
c =====================================================================
c     do the second pass, required for std cases
c
c =====================================================================
c
c --------------------------------------------------------------------
c     determine the number of timesteps per year
c --------------------------------------------------------------------
c     
      nstep_per_year = int(1.d0/ts)
c
      write(*,*) '--COMPUTE_MC_RUNAVG: number of ts per year for ',
     $     ' standard case : ',
     $     nstep_per_year
c
c --------------------------------------------------------------------
c     determine the first and last steps
c --------------------------------------------------------------------
c     
      nstep_per_halfyr = int(nstep_per_year/2)
      nstep_start = nstep_per_halfyr+1
      nstep_stop = nstep - nstep_per_halfyr
      if (nstep_stop .le. nstep_start) then
         write(*,*) '--COMPUTE_RUNAVG: running average not possible,',
     $        ' because record is not long enough '
         write(*,*) '     nstep_start, nstep_stop : ',nstep_start,
     $        nstep_stop
         error = .true.
         goto 9900
      endif
c
      write(*,*) '--COMPUTE_MC_RUNAVG: nstep_start,nstep_stop : ',
     $     nstep_start,nstep_stop
c
c --------------------------------------------------------------------
c     pad at the beginning
c --------------------------------------------------------------------
c     
      
      do k = 1,nstep_start-1
         do i = 1,nproc_rates
            D_C_std_runavg(i,k) = 0.0d0
         enddo
         do i = 1,nproc_fluxes
            fluxes_std_runavg(i,k) = 0.0d0
         enddo
      enddo
             
c
c --------------------------------------------------------------------
c     do the central (the real) part
c --------------------------------------------------------------------
c     
      do k = nstep_start,nstep_stop
         do i = 1,nproc_rates
            sum1 = 0.d0
             do kk = k-nstep_per_halfyr,k+nstep_per_halfyr
               sum1 = sum1 + D_C_std(i,kk)
            enddo
            D_C_std_runavg(i,k) = sum1/(2*nstep_per_halfyr+1)
         enddo
         do i= 1, nproc_fluxes
            sum1 = 0.0d0
            do kk = k-nstep_per_halfyr,k+nstep_per_halfyr
               sum1 = sum1 + fluxes_std(i,kk)
            enddo
            fluxes_std_runavg(i,k) = sum1/(2*nstep_per_halfyr+1)
         enddo
      enddo
c
c --------------------------------------------------------------------
c     pad at the end
c --------------------------------------------------------------------
c     
      
      do k = nstep_stop+1,nstep
         do i = 1,nproc_rates
            D_C_std_runavg(i,k) = 0.0d0
         enddo
         do i = 1,nproc_fluxes
            fluxes_std_runavg(i,k) = 0.0d0
         enddo
      enddo
             
cc
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c
      write(*,*) '--COMPUTE_MC_RUNAVG: terminated.'
c    
      return
c
c --------------------------------------------------------------------
c     error handling
c --------------------------------------------------------------------
c     
 9900 write(*,*) '--COMPUTE_MC_RUNAVG: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











