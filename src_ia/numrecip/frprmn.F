c =====================================================================
c
c     SUBROUTINE frprmn
c
c     PURPOSE:  calculate minimum of a function using the conjugate 
c                 gradient technique in multi-dimensions
c
c               Given a starting point p that is a vector of length n,
c               Fletcher-Reeves-Polak-Ribiere minimization is performed
c               on a function func, using its gradient as calculated
c               by a routine dfunc. The convergence tolerance on the
c               function value is input as ftol. Returned quantities
c               are p (the location of the minimum), iter (the number
c               of iterations that were performed), and fret (the
c               minimum value of the function). The routine linmin is
c               called to perform line minimization.
c               parameters: nmas is the maximum anticipated value of
c                  n; itmax is the maximum allowed number of iterations
c                  eps is a small number to rectify special case of
c                  converging exactly zero function value.
c
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 416ff
c
c     SUBROUTINES: needs routines dfunc,func,linmin
c                      dfunc  : evaluates derivative at point p
c                      func   : evaluate function at point p
c                      linmin : calculates minimum along a line
c
c     VARIABLES:
c 
c     IN:     p    : start vector (dimension n)
c             n    : dimension of p
c             ftol : convergence tolerance  
c  
c     OUT:    iter : number of iterations
c             fret : value of func at minimum
c             p    : minimum vector
c
c     REVISIONS:
c
c     date     author   remarks
c
c     12.08.96   ng     first implementation
c 
c =====================================================================
c
      subroutine frprmn(p,n,ftol,iter,fret)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      implicit none
      integer iter,n,nmax,itmax
      double precision fret,ftol,p(n),eps,func
      parameter (nmax=50,itmax=200,eps=1.e-10)
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      integer its,j
      double precision dgg,fp,gam,gg,g(nmax),h(nmax),xi(nmax)
c
c --------------------------------------------------------------------
c     external routines
c --------------------------------------------------------------------
c
      external func
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--FRPRMN: finding minimum...'
c
c --------------------------------------------------------------------
c     initializations
c --------------------------------------------------------------------
c
      fp=func(p)
      call dfunc(p,xi)
      do 11 j=1,n
         g(j)=-xi(j)
         h(j)=g(j)
         xi(j)=h(j)
 11   continue
c
c --------------------------------------------------------------------
c     loop over iterations
c --------------------------------------------------------------------
c
      do 14 its=1,itmax
c
         write(*,*) '  iteration # ',its
c         do j=1,n
c            write(*,*) '    p(j) : ',p(j)
c            write(*,*) '   xi(j) : ',xi(j)
c         enddo
c
         iter=its
         call linmin(p,xi,n,fret)
c
c --------------------------------------------------------------------
c     return if within tolerance
c --------------------------------------------------------------------
c
         if (2.d0*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+eps))
     $        return
c
c --------------------------------------------------------------------
c     else continue with the iteration
c --------------------------------------------------------------------
c
         fp=func(p)
         call dfunc(p,xi)
         gg=0.
         dgg=0.
         do 12 j=1,n
            gg=gg+g(j)**2
c            dgg=dgg+xi(j)**2              ! for Fletcher-Reeves
            dgg=dgg+(xi(j)+g(j))*xi(j)     ! for Polak-Ribiere 
 12      continue
c
c --------------------------------------------------------------------
c     return if gradient gg is zero. unlikely but possible
c --------------------------------------------------------------------
c
         if(gg.eq.0.)return
         gam=dgg/gg
         do 13 j=1,n
            g(j)=-xi(j)
            h(j)=g(j)+gam*h(j)
            xi(j)=h(j)
 13      continue
 14   continue
c
c --------------------------------------------------------------------
c     maximum iteration exceeded
c --------------------------------------------------------------------
c
      pause 'frprmn maximum iterations exceeded'
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software +50.
