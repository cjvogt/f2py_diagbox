c =====================================================================
c
c     DOUBLE PRECISION FUNCTION DELTA_TO_C13
c 
c     PURPOSE: this function calculates the 13C concentration from the
c                 dC13 value following
c          
c                C13 = C * rStd * (dC13 + 1) /
c                      (1 + rStd* dC13 + rs)
c     VARIABLES:
c 
c     IN : dC13           : delta 13C value [per mil]
c          C              : total C (12C + 13C) [umol/kg]
c     OUT: delta_to_C13   : C13 concentration [umol/kg]
c 
c     REVISION:
c
c     date     author  remarks
c   
c      1.08.96   ng    first implementation
c
c =====================================================================
c
      double precision function delta_to_C13(dC13,C)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "c13_params.h"
c
      double precision dC13,C
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
      delta_to_c13 = C * r_std * (dC13 * 1.d-3 + 1.d0)/
     $     (1.d0 + r_std * dC13 * 1.d-3 + r_std)
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

