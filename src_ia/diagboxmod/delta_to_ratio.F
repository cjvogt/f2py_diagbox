c =====================================================================
c
c     DOUBLE PRECISION FUNCTION DELTA_TO_RATIO
c 
c     PURPOSE: this function calculates the 13C/12C ratio from the
c                 dC13 value following
c          
c                ratio = rStd * (dC13 * 0.001d0 + 1.d0)
c     VARIABLES:
c 
c     IN : dC13           : delta 13C value [per mil]
c     OUT: delta_to_ratio : 13C/12C ratio [no units]
c 
c     REVISION:
c
c     date     author  remarks
c   
c      1.08.96   ng    first implementation
c
c =====================================================================
c
      double precision function delta_to_ratio(dC13)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "c13_params.h"
c
      double precision dC13
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
      delta_to_ratio = r_std * (dC13 * 0.001d0 + 1.d0)
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

