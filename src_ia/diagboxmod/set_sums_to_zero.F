c =====================================================================
c
c     SUBROUTINE SET_SUMS_TO_ZERO
c 
c     PURPOSE: this routine sets all variables that are used for
c                 calculating averages to zero
c
c     VARIABLES:
c 
c     OUT: sum_fract    : sum for calculating average fractionation factor
c          sum_kex      : sum for calculating average kex
c          sum_dCdz     : sum for calculating average vertical dic gradient
c          sum_ddC13dz  : sum for calculating average vertical dC13 gradient
c          sum_dc13_prc : sum for calculating sum of dc13 value of prc
c          nsum_dc13_prc: number of values used in sum_dc13_prc
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
      subroutine set_sums_to_zero(
     $     sum_fract,sum_kex,
     $     sum_dcdz, sum_ddC13dz,
     $     sum_dc13_prc,nsum_dc13_prc)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "averages.h"
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer i
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      sum_fract = 0.0d0
      sum_kex = 0.0d0
      sum_dcdz = 0.0d0
      sum_ddC13dz = 0.0d0
      do i = 1,4
         sum_dc13_prc(i) = 0.d0
         nsum_dc13_prc(i) = 0
      enddo
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

