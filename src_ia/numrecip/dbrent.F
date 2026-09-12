c =====================================================================
c
c     SUBROUTINE dbrent
c
c     PURPOSE:  isolates the minimum to a fractional precision og about 
c                tol
c
c               given a function f and its derivative df, and given
c               a bracketing triplet of abscissas ax,bx,cx (such that
c               bx is between ax and cx, and f(bx) is less than both
c               f(ax) and f(cx)), this routine isolates the minimum
c               to a fractional precision of about tol using a 
c               modification of Brent's method that uses derivatives.
c               The abscissa of the minimum is returned as xmin, and
c               the minimum function value is returned as dbrent, the
c               returned function value
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 400ff
c
c     SUBROUTINES: needs routines f and df
c                      f    : evaluates function value at x
c                      df   : evaluates value of derivative at x
c
c     VARIABLES:
c 
c     IN:      ax      : initial abscissa bracket point a
c              bx      : initial abscissa bracket point b
c              cx      : initial abscissa bracket point c
c              f       : name of function f
c              df      : name of derivative function df
c              tol     : tolerance
c  
c     OUT:     xmin    : abscissa of minimum
c              dbrent  : returned function value
c
c     REVISIONS:
c
c     date     author   remarks
c
c     12. 8.96   ng     first implementation
c 
c =====================================================================
c
      double precision function dbrent(ax,bx,cx,f,df,tol,xmin)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      double precision ax,bx,cx,tol,xmin,df,f
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      double precision zeps
      integer itmax
      parameter (itmax=100,zeps=1.0d-10)
      integer iter
      double precision a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw,fx
      double precision olde,tol1,tol2,u,u1,u2,v,w,x,xm
      logical ok1,ok2
c
c --------------------------------------------------------------------
c     external routines
c --------------------------------------------------------------------
c
      external df,f
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--DBRENT: finding minimum ...'

      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      dx=df(x)
      dv=dx
      dw=dx
      do 11 iter=1,itmax
        xm=0.5d0*(a+b)
        tol1=tol*abs(x)+zeps
        tol2=2.d0*tol1
        if(abs(x-xm).le.(tol2-.5d0*(b-a))) goto 3
        if(abs(e).gt.tol1) then
          d1=2.d0*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0.d0).and.(dx*d1.le.0.d0)
          ok2=((a-u2)*(u2-b).gt.0.d0).and.(dx*d2.le.0.d0)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(abs(d1).lt.abs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(abs(d).gt.abs(0.5d0*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0.) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5*e
2       if(abs(d).ge.tol1) then
          u=x+d
          fu=f(u)
        else
          u=x+sign(tol1,d)
          fu=f(u)
          if(fu.gt.fx)goto 3
        endif
        du=df(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software +50.
