c =====================================================================
c
c     SUBROUTINE CALC_HESSIAN
c 
c     PURPOSE: this routine calculates the hessian matrix at the
c                parameter optimum. The hessian matrix is defined
c                as :
c                        hessian_ij = d2costfn/(dx_i dx_j)
c
c              attention ! it is necessary to compute the
c                partial derivatives with the original parameter
c                values and not with the transformed ones !
c
c     VARIABLES:
c
c     IN : p_opt      : optimum parameter set
c          p_init     : inital parameter set
c          sg_p       : uncertainty of parameters
c 
c     OUT: hessian    : hessian matrix
c
c     REVISION:
c
c     date     author  remarks
c   
c     15.08.96   ng    first implementation
c
c =====================================================================
c
      subroutine calc_hessian(p_opt,p_init,sg_p,hessian)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "conj_grad.h"
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer i,j,np
      double precision delta_i,delta_j
      double precision p_l(nparams),p_r(nparams)
      double precision x_r,x_l,dx_i_r,dx_i_l
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision func
c      
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--CALC_HESSIAN : calculating the hessian matrix...'
c
c --------------------------------------------------------------------
c     loop over the upper diagonal matrix to calculate its values
c --------------------------------------------------------------------
c        
      do i = 1,nparams
c
c --------------------------------------------------------------------
c     estimate size of delta-i
c --------------------------------------------------------------------
c     
            if (p_opt(i) .ne. 0.d0) then
               delta_i = 0.0001 * p_opt(i)
            else
               delta_i = 0.0001
            endif

         do j = 1,i
c
c --------------------------------------------------------------------
c     estimate size of delta-j
c --------------------------------------------------------------------
c     
            if (p_opt(j) .ne. 0.d0) then
               delta_j = 0.0001 * p_opt(j)
            else
               delta_j = 0.0001
            endif
c
c --------------------------------------------------------------------
c     calculate first derivative at pj-delta_j : centered differences
c --------------------------------------------------------------------
c     
            do np = 1,nparams
               p_l(np) = p_opt(np)
               p_r(np) = p_opt(np)
            enddo
            p_l(j) = p_l(j) - delta_j
            p_r(j) = p_r(j) - delta_j
            p_l(i) = p_l(i) - delta_i
            p_r(i) = p_r(i) + delta_i
c
            x_r = func(p_r)
            x_l = func(p_l)
c
c --------------------------------------------------------------------
c     attention: remove normalization for having the true derivative !
c --------------------------------------------------------------------
c            
            dx_i_l = 0.5d0 /(sg_p(i)*delta_i) 
     $           * (x_r - x_l)
c
c --------------------------------------------------------------------
c     calculate first derivative at pj+delta_j : centered differences
c --------------------------------------------------------------------
c     
            do np = 1,nparams
               p_l(np) = p_opt(np)
               p_r(np) = p_opt(np)
            enddo
            p_l(j) = p_l(j) + delta_j
            p_r(j) = p_r(j) + delta_j
            p_l(i) = p_l(i) - delta_i
            p_r(i) = p_r(i) + delta_i
c
            x_r = func(p_r)
            x_l = func(p_l)
c
c --------------------------------------------------------------------
c     attention: remove normalization for having the true derivative !
c --------------------------------------------------------------------
c            
            dx_i_r = 0.5d0 / (sg_p(i)*delta_i) 
     $           * (x_r - x_l)
c
c --------------------------------------------------------------------
c     calculate second derivative by centered differences
c      take also normalization into account !
c --------------------------------------------------------------------
c  
            hessian(i,j) = 0.5d0 / (sg_p(j)*delta_j)
     $           * (dx_i_r - dx_i_l)
c
c --------------------------------------------------------------------
c     use symmetry for lower diagnonal matrix
c --------------------------------------------------------------------
c
            if (j .ne. i) then
               hessian(j,i) = hessian(i,j)
            endif
c
c --------------------------------------------------------------------
c     end of loops over i and j
c --------------------------------------------------------------------
c
         enddo
      enddo
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

