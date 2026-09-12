c =====================================================================
c
c     SUBROUTINE INVERT_HESSIAN
c 
c     PURPOSE: this routine calculates the inverse of the
c                hessian matrix using singular value decomposition.
c              it is important that the singular values (stored in w)
c                are edited before the inverse is calculated
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 51ff
c
c
c     SUBROUTINES: needs routine svdcmp
c
c
c     VARIABLES:
c
c     IN : hessian    : hessian matrix
c          nparams    : dimension of hessian matrix (is a square matrix)
c 
c     OUT: covm       : inverse of the hessian matrix is equal to
c                          the covariance matrix
c
c     REVISION:
c
c     date     author  remarks
c   
c     15.08.96   ng    first implementation
c
c =====================================================================
c
      subroutine invert_hessian(hessian,covm)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "conj_grad.h"
      double precision w(nparams), v(nparams,nparams)
      double precision a(nparams,nparams)

c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer i,j,k
      double precision tol,thresh,wmax,s
      parameter (tol=1.e-6)
      double precision tmp(nparams,nparams),check(nparams,nparams)
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
c      
c =====================================================================
c     begin of executable code
c =====================================================================
c
      write(*,*) '--INVERT_HESSIAN : inverting the hessian matrix...'
c
c --------------------------------------------------------------------
c     performs singular value decomposition
c       in order to keep original hessian matrix, assign it to
c       matrix a, the return matrices are a,w,v, so that the
c       singular value decomposition definition is satisfied
c                   a = u * w * vt 
c       where vt is the transposed matrix of v.
c --------------------------------------------------------------------
c
      do i = 1,nparams
         do j = 1,nparams
            a(i,j) = hessian(i,j)
         enddo
      enddo
c
c      write(*,*) '--INVERT_HESSIAN: calling svdcmp...'
      call svdcmp(a,nparams,nparams,nparams,nparams,w,v)
c      write(*,*) '--INVERT_HESSIAN: back from svdcmp...'
c
c --------------------------------------------------------------------
c     info for debugging
c --------------------------------------------------------------------
c        
c      write(*,*) '     values of w right after svdcmp'
c      do j = 1,nparams
c         write(*,*) 'w(j) : ',w(j)
c      enddo
c
c --------------------------------------------------------------------
c     edit the singular values, given tol from the parameter statement
c       this is necessary
c --------------------------------------------------------------------
c
      wmax = 0.d0
      do j=1,nparams
         if (w(j) .gt. wmax) wmax = w(j)
      enddo
      thresh = tol*wmax
      do j=1,nparams
         if (w(j) .lt. thresh) w(j)=0.
      enddo
c
c --------------------------------------------------------------------
c     info for debugging
c --------------------------------------------------------------------
c        
c      write(*,*) '     values of w right after editing'
c      do j = 1,nparams
c         write(*,*) 'w(j) : ',w(j)
c      enddo
c
c --------------------------------------------------------------------
c     calculate the inverse of a by
c       a-1 = v * diag(1/w) * ut
c --------------------------------------------------------------------
c
c --------------------------------------------------------------------
c     calculate first diag(1/w) * ut
c       i : index for rows ; k or j : index for columns
c --------------------------------------------------------------------
c        
      do i = 1,nparams
         do k = 1,nparams
            if (w(i) .ne. 0.d0) then
               tmp(i,k) = a(k,i) / w(i) 
            else
               tmp(i,k) = 0.d0
            endif
         enddo
      enddo
c
c --------------------------------------------------------------------
c     calculate covariance matrix by multiplication with v
c --------------------------------------------------------------------
c        
      do j=1,nparams
         do i=1,nparams
            covm(i,j) = 0.0d0
         enddo
         do k=1,nparams
            do i=1,nparams
               covm(i,j) = covm(i,j)+v(i,k)*tmp(k,j)
            enddo
         enddo
      enddo

c
c --------------------------------------------------------------------
c     make a check: covm * hessian = 1
c --------------------------------------------------------------------
c        
c      do j=1,nparams
c         do i=1,nparams
c            check(i,j)=0.0d0
c         enddo
c         do k=1,nparams
c            do i=1,nparams
c               check(i,j)=check(i,j) + covm(i,k) * hessian(k,j)
c            enddo
c         enddo
c      enddo
cc
c      write(*,*) '--INVERT_HESSIAN: check for inversion...'
c      s = 0.0d0
c      do i = 1,nparams
c         write(*,8000) i,(check(i,j),j=1,nparams)
c 8000    format(' row : ',i1,2x,8(f8.4,x))
c         do j = 1,nparams
c            s = s + check(i,j)
c         enddo
c      enddo
c      s = s / (float(nparams))
cc
c      if (s .ne. 1.d0) then
c         write(*,*) '--INVERT_HESSIAN: problems !!'
c         write(*,*) '     value of check matrix : ',s
c      endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

