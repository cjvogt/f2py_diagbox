c =====================================================================
c
c     SUBROUTINE WRITE_MC_RATES_TOFILE
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
c          D_c_std         : rates of change of the different processes
c          c_sim_std       : simulated seasonal cycle of c
c          c_obs           : observed seasonal cycle of C
c
c          mean_d_c        : mean d_c of mc runs
c          mean_c_sim      : mean c_sim of mc runs
c          wmean_d_c       : weighted mean d_c of mc runs
c          wmean_c_sim     : weighted mean c_sim of mc runs
c          sg_d_c          : 1-sg uncertainty of d_c
c          sg_c_sim        : 1-sg uncertainty of c_sim
c          wsg_d_c         : weighted 1-sg uncertainty of d_c
c          wsg_c_sim       : weighted 1-sg uncertainty of c_sim
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     02.09.96   ng     first implementation, only weighted intervals
c                         are written to unit
c     05.05.02   ng     copied changed version over from Holger,
c                         adds another output file
c 
c =====================================================================
c
      subroutine write_mc_rates_tofile(unit,tstart,ts,
     $     nstep2,rsmpl_step,
     $     D_c_std,c_sim_std,c_obs,
     $     mean_d_c,mean_c_sim,
     $     wmean_d_c,wmean_c_sim,
     $     sg_D_c,sg_c_sim,
     $     wsg_D_c,wsg_c_sim,
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
cBX      integer k,k2,i
      integer k,k2,i,unit2
      double precision dcdt(nproc_rates)
      double precision dcdt_lb(nproc_rates),dcdt_ub(nproc_rates)
cBX      character*27 tmp(nproc_rates)
      character*36 tmp(nproc_rates)
cBX-a
      character*30 filename
      character*24 tmp2(nproc_rates)
cBX-e
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      error = .false.
      time = tstart
cBX-a
c
c --------------------------------------------------------------------
c      add another output file with differences instead of upper
c        and lower limits
c --------------------------------------------------------------------
c
         unit2 = 20
         filename = 'mc_rates_diff.mc'
         open (unit2,file = filename, form = 'formatted', 
     $        status = 'unknown')

cBX-e
c
c --------------------------------------------------------------------
c     loop over all simulated values ! here only resampled values
c --------------------------------------------------------------------
c     
      do k2 = 1,nstep2
c
         k = (k2-1)*rsmpl_step + 1
c
c --------------------------------------------------------------------
c     change units of D_C (now it is umol kg-1) to umol kg-1 day-1
c      calculate also upper and lower bounds of weighted intervals
c --------------------------------------------------------------------
c
         do i = 1,nproc_rates
            dcdt(i) = d_C_std(i,k) / ts
            dcdt_lb(i) = (d_C_std(i,k)-wsg_d_c(i,k2))/ts
            dcdt_ub(i) = (d_C_std(i,k)+wsg_d_c(i,k2))/ts
         enddo
c
c --------------------------------------------------------------------
c     write to temporary array for output
c --------------------------------------------------------------------
c
         do i = 1,nproc_rates
cBX            write(tmp(i),'(3(f8.5,1x))') 
cBX            write(tmp(i),'(3(f10.5,1x))') 
            write(tmp(i),'(3(f11.4,1x))') 
     $           dcdt(i),dcdt_lb(i),dcdt_ub(i)
cBX-a
            write(tmp2(i),'(2(f11.4,1x))') 
     $           dcdt(i),wsg_d_c(i,k2)/ts
         enddo
cBX-e
c
c --------------------------------------------------------------------
c     write to unit
c --------------------------------------------------------------------
c
         write(unit,9000) 
     $        time,
     $        (tmp(i),i=1,nproc_rates),
     $        c_sim_std(k),
     $        c_sim_std(k)-wsg_c_sim(nstep2),
     $        c_sim_std(k)+wsg_c_sim(nstep2),
     $        c_obs(k)
c
cBX 9000    format(f9.4,1x,7(a26),4(f7.2,1x))
 9000    format(f9.4,1x,7(a36),1x,4(f7.2,1x))
cBX-a
c
c --------------------------------------------------------------------
c     write to 2nd unit for 2nd file
c --------------------------------------------------------------------
c
         write(unit2,9001) 
     $        time,
     $        (tmp2(i),i=1,nproc_rates),
     $        c_sim_std(k),
     $        wsg_c_sim(nstep2),
     $        c_obs(k)
c
cBX 9000    format(f9.4,1x,7(a26),4(f7.2,1x))
 9001    format(f9.4,1x,7(a24),1x,3(f7.2,1x))
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
 9900 write(unit,*) '--WRITE_MC_RATES_TOFILE: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











