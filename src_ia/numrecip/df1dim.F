c =====================================================================
c
c     DOUBLE PRECISION FUNCTION df1dim
c
c     PURPOSE:  routine which constructs an artificial derivative 
c               function of
c               one variable. needed by linmin to communicate with
c               mnbrak, brent and dbrent
c
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 417ff
c
c     SUBROUTINES: need subroutine dfunc
c                    dfunc : evaluates the derivative of the func at
c                               point x
c
c     VARIABLES:
c 
c     IN:      x       : input value x
c  
c     OUT:     df1dim  : value of derivative of func at x, 
c                         communicates with linmin
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
      double precision function df1dim(x)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      implicit none
      integer nmax
      double precision x
      parameter (nmax=50)
      integer j,ncom
      double precision df(nmax),pcom(nmax),xicom(nmax),xt(nmax)
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
      call dfunc(xt,df)
      df1dim=0.
      do 12 j=1,ncom
        df1dim=df1dim+df(j)*xicom(j)
12    continue
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software +50.
