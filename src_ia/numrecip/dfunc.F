c =====================================================================
c
c     SUBROUTINE DFUNC
c
c     PURPOSE: calculate the derivative of the
c                 function value at the vector p, calculates
c                 the derivative numerically
c     
c     VARIABLES:
c 
c     IN:   p            : input array of parameters
c     OUT:  xi           : vector of partial derivatives
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     13.08.96   ng     first implementation
c 
c =====================================================================
c
      subroutine dfunc(p,xi)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
#include "conj_grad.h"
      double precision xi(nparams)
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      integer i,np
      double precision delta
c
      double precision p_l(nparams),p_r(nparams)
      double precision x_l,x_r
c
c --------------------------------------------------------------------
c     external functions
c --------------------------------------------------------------------
c
      double precision func
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*)
c      write(*,*) '--DFUNC: calculating derivative'
c      write(*,*) '           p : ',p
c
c --------------------------------------------------------------------
c     calculate the partial derivatives (use centered differences)
c       2.nd order accurate
c     loop over all parameters
c --------------------------------------------------------------------
c     
      do np = 1,nparams
c
c --------------------------------------------------------------------
c     estimate size of delta
c --------------------------------------------------------------------
c     
         if (p(np) .ne. 0.d0) then
            delta = 0.0001 * p(np)
         else
            delta = 0.0001
         endif
c
c --------------------------------------------------------------------
c     calculate the derivative
c --------------------------------------------------------------------
c     
         do i = 1,nparams
            p_l(i) = p(i)
            p_r(i) = p(i)
         enddo
c
         p_l(np) = p_l(np) - delta
         p_r(np) = p_r(np) + delta
c
         x_r = func(p_r)
         x_l = func(p_l)
c
         xi(np) = 0.5d0 / delta * (x_r - x_l)
c
c --------------------------------------------------------------------
c     give info
c --------------------------------------------------------------------
c     
c         write(*,*) '--DFUNC: values of np,xi(np) : ',np,xi(np)
c         write(*,*) '                   x_r,x_l   : ',x_r,x_l
c
c --------------------------------------------------------------------
c     end of loop over parameters
c --------------------------------------------------------------------
c     
      enddo
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c     
c      write(*,*) '--DFUNC: finished'
      return

      end











