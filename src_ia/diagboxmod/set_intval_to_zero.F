c =====================================================================
c
c     SUBROUTINE SET_INTVAL_TO_ZERO
c 
c     PURPOSE: this routine sets all variables that are used for
c                 doing the integrations to zero
c
c     VARIABLES:
c 
c     OUT: Int_Fluxes : two dim array holding the integrated Fluxes of
c                          the different processes during the
c                          different periods
c          Int_Rates  : two dim array holding the integrated Rates of
c                          change of the different processes during
c
c     REVISION:
c
c     date     author  remarks
c   
c     31.07.96   ng    first implementation
c     18.10.99   ng    adapted for interannual variability runs
c
c =====================================================================
c
      subroutine set_intval_to_zero(Int_Fluxes,Int_Rates)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "periods.h"
# include "processes.h"
# include "intval.h"
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer ny,i
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      do ny = 1,nyearmax
         do i = 1,9
            int_fluxes(ny,i) = 0.0d0
         enddo
         do i = 1,6
            int_rates(ny,i) = 0.0d0
         enddo
      enddo
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

