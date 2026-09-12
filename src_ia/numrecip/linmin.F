c =====================================================================
c
c     SUBROUTINE linmin
c
c     PURPOSE:  calculate minimum along a line
c
c               Given an n-dimensional point p(1:n) and an n-dimensional
c               direction xi(1:n), moves and resets p to where the
c               function func(p) takes on a minimum along the direction
c               xi from p, and replaces xi by the actual vector 
c               diplacement that p was moved. Also returns as fret the
c               value of func at the returned location p. This is
c               actually all accomplished by calling the routine mnbrak and 
c               brent
c
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 412
c
c     SUBROUTINES: needs routines dbrent,f1dim,df1dim,mnbrak
c                      dbrent  : 
c                      f1dim   : 
c                      df1dim  : 
c                      mnbrak  :
c
c     VARIABLES:
c 
c     IN:      p       : n-dimensional point
c              xi      : n-dimensional direction
c              n       : dimension
c  
c     OUT:     p       : modified n-dimensional point
c              xi      : modified direction (how much p was moved)
c              fret    : value of func at new location p
c
c     REVISIONS:
c
c     date     author   remarks
c
c     12. 8.96   ng     first implementation, note that this is the
c                        version where the call to brent has been
c                        replaced by a call to dbrent
c 
c =====================================================================
c
      subroutine linmin(p,xi,n,fret)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      implicit none

      integer n,nmax
      double precision fret,p(n),xi(n),tol
      parameter (nmax=50,tol=1.e-4)
      integer ncom
      double precision ax,bx,fa,fb,fx,xmin,xx,pcom(nmax),xicom(nmax)
      common /f1com/ pcom,xicom,ncom
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      integer j
c
c --------------------------------------------------------------------
c     external routines
c --------------------------------------------------------------------
c
      double precision dbrent,f1dim,df1dim
      external f1dim,df1dim
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--LINMIN: starting to find minimum on a line',
c     $     ' segment...'
c
c --------------------------------------------------------------------
c     set up the common block
c --------------------------------------------------------------------
c
      ncom=n
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
11    continue
c
c --------------------------------------------------------------------
c     initial guess for brackets
c --------------------------------------------------------------------
c
      ax = 0.d0
      xx = 1.d0
c
c --------------------------------------------------------------------
c     get brackets
c --------------------------------------------------------------------
c
c      write(*,*) '--LINMIN: calling mnbrak...'
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
c
c --------------------------------------------------------------------
c     call dbrent instead of brent here
c --------------------------------------------------------------------
c      fret=brent(ax,xx,bx,f1dim,tol,xmin)
c      write(*,*) '--LINMIN: calling dbrent...'
      fret=dbrent(ax,xx,bx,f1dim,df1dim,tol,xmin)

c
c --------------------------------------------------------------------
c     give info about minimum found
c --------------------------------------------------------------------
c
c      write(*,*) '--LINMIN: back from dbrent: xmin : ',xmin
c
c --------------------------------------------------------------------
c     construct the vector results for return
c --------------------------------------------------------------------
c
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software +50.
