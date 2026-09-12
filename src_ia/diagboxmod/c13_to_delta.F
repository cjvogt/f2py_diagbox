c =====================================================================
c
c     DOUBLE PRECISION FUNCTION C13_TO_DELTA
c 
c     PURPOSE: this function calculates the dC13C value from the
c                 13C concentration following
c          
c                dC13 = 1/rstd * 13C/(C-13C) - rstd
c     VARIABLES:
c 
c     IN : C13            : C13 concentration [umol/kg]
c          C              : total C (12C + 13C) [umol/kg]
c     OUT: 13C_to_delta   : delta 13C value [per mil]
c 
c     REVISION:
c
c     date     author  remarks
c   
c      2.08.96   ng    first implementation
c
c =====================================================================
c
      double precision function C13_to_delta(C13,C)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "c13_params.h"
c
      double precision C13,C
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      c13_to_delta = 1.d3 * (C13/(r_std*(C-C13)) - 1.d0)
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

