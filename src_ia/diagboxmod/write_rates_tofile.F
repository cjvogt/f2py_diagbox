c =====================================================================
c
c     SUBROUTINE WRITE_RATES_TOFILE
c
c     PURPOSE: writes the rates of change to a file specified by
c                the unit number unit
c     
c
c     VARIABLES:
c 
c     IN:  unit            : unit number to write to
c          tstart          : start time of integration
c          ts              : time step
c          nstep           : number of steps
c          D_C             : array holding the rates of change
c          c_obs           : array holding the observed seasonal cycle of C
c          c_sim           : array holding the simulated seasonal cycle of C
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     03.08.96   ng     first implementation
c      5.08.96   ng     included c_obs and c_sim
c 
c =====================================================================
c
      subroutine write_rates_tofile(unit,tstart,ts,nstep,
     $        D_c,c_obs,c_sim,error)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
c
      logical error
      integer unit
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      integer k,i
      double precision dcdt(nproc_rates)
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
            dcdt(i) = d_C(i,k) / ts
         enddo
c
c --------------------------------------------------------------------
c     write to file
c --------------------------------------------------------------------
c
         write(unit,9000) 
     $        time,
     $        (dcdt(i),i=1,nproc_rates),
     $        c_sim(k),c_obs(k)
c
 9000    format(f8.3,1x,7(f8.3,1x),2(f7.2,1x))
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
 9900 write(unit,*) '--WRITE_RATES_TOFILE: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











