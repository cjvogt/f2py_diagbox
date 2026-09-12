c =====================================================================
c
c     SUBROUTINE WRITE_FLUXES_TOFILE
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
c          fluxes          : array holding the fluxes
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     03.08.96   ng     first implementation
c     07.06.00   ng     adapted for interannual variability run, ATTENTION
c                        time units have changed from d-1 to yr-1
c 
c =====================================================================
c
      subroutine write_fluxes_tofile(unit,tstart,ts,nstep,
     $        fluxes,error)
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
      double precision fluxes_gd(nproc_fluxes)
      double precision sec_day,g_mol,day_year
      parameter (sec_day = 24.d0 * 3600.d0, g_mol = 12.01d0)
      parameter (day_year = 365.d0)

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
            fluxes_gd(i) = fluxes(i,k) * sec_day * day_year 
         enddo
c
c --------------------------------------------------------------------
c     write to file
c --------------------------------------------------------------------
c
         write(unit,9000) 
     $        time,
     $        (fluxes_gd(i),i=1,nproc_fluxes)
c
 9000    format(f8.3,1x,7(f8.4,1x),3(f10.4,1x))
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
 9900 write(unit,*) '--WRITE_FLUXES_TOFILE: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











