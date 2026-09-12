c =====================================================================
c
c     SUBROUTINE INIT_VAL
c 
c     PURPOSE: this routine initializes the values needed for
c                doing the simulation
c
c     VARIABLES:
c 
c     IN : time         : actual time
c          dc13_coeff      : character specifying which dC13 coefficients
c                              are used       
c          h_coeff      : character specifying which coefficients for
c
c     OUT: C_old        : mixed layer C concentration at time t
c          C13_old      : mixed layer C13 concentration at time t
c          h_old        : mixed layer depth at time t
c          C_sim        : mixed layer C concentration at start time
c          error        : logical
c 
c     REVISION:
c
c     date     author  remarks
c   
c     31.07.96   ng    first implementation
c      5.09.96   ng    included possibility to do backwards integration
c     15.10.99   ng    adapted for interannual variability run, time
c                         is now in years,added error
c     07.06.00   ng    avoid h less than 15
c
c =====================================================================
c
      subroutine init_val(time,dc13_coeff,h_coeff,
     $     C_old,c13_old,h_old,c_sim,error)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
# include "harm_coeff.h"
# include "diffent_params.h"
      logical error
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision harmonic,delta_to_c13,spline
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision dc13_tmp
       
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      C_old = harmonic(time,c_sdic,7) + spline(time,'sdic ',error)
c
      dc13_tmp = harmonic(time,c_dc13,7) + spline(time,'dc13 ',error)
      c13_old = delta_to_c13(dc13_tmp,C_old)
c
      h_old = harmonic(time,c_h,7) + spline(time,'h    ',error)
c
      if (h_old .lt. 15.) then
         h_old = 15.
      endif

c
      C_sim(1) = C_old
c
      if (error) goto 9900
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
 9900 write(*,*) '--INIT_VAL: an error occured in spline.'
      error = .true.
      return      
      end

