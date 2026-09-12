c =====================================================================
c
c     SUBROUTINE WRITE_SOLUTION_TOFILE
c
c     PURPOSE: writes all other solutions to the model to a file
c     
c
c     VARIABLES:
c 
c     IN:  unit            : unit number to write to
c          tstart          : start time of integration
c          ts              : time step
c          nstep           : number of steps
c          sol_scheme      : solution scheme
c          u_var           : computed u 
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     19.10.99   ng     first implementation
c 
c =====================================================================
c
      subroutine write_solution_tofile(unit,tstart,ts,nstep,
     $        sol_scheme,u_var,error)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "diagboxmod.h"
# include "simulation.h"
# include "solution.h"
c
      logical error
      integer unit
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      integer k,i
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
c     write to file
c --------------------------------------------------------------------
c
         write(unit,9000,err=9900) 
     $        time,
     $        u_var(k)
c
 9000    format(f9.3,1x,f8.4)
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
 9900 write(unit,*) '--WRITE_SOLUTION_TOFILE: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











