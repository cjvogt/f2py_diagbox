c =====================================================================
c
c     DOUBLE PRECISION FUNCTION f1dim
c
c     PURPOSE:  routine which constructs an artificial function of
c               one variable. needed by linmin to communicate with
c               mnbrak, brent and dbrent
c
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 413ff
c
c     SUBROUTINES: need subroutine func
c                    func : evaluates the value at point x
c
c     VARIABLES:
c 
c     IN:      x       : input value x
c  
c     OUT:     f1dim   : value of func at x, communicates with linmin
c                         through a common block
c
c     REVISIONS:
c
c     date     author   remarks
c
c     12. 8.96   ng     first implementation
c 
c =====================================================================
c
      double precision function f1dim(x)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      implicit none
      integer nmax
      double precision func,x
      parameter (nmax=50)
      integer j,ncom
      double precision pcom(nmax),xicom(nmax),xt(nmax)
      common /f1com/ pcom,xicom,ncom
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
c
c --------------------------------------------------------------------
c     external routines
c --------------------------------------------------------------------
c
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
c
c --------------------------------------------------------------------
c     info
c --------------------------------------------------------------------
c 
c      write(*,*) '--F1DIM:  parameter values : '
c      do j = 1,ncom
c         write(*,*) '                j, x : ',j,x
c         write(*,*) '   pcom(j), xicom(j) : ',pcom(j),xicom(j)
c         write(*,*) '               xt(j) : ',xt(j)
c      enddo
c
c --------------------------------------------------------------------
c     call function
c --------------------------------------------------------------------
c 
      f1dim=func(xt)
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software +50.
