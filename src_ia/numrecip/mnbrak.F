c =====================================================================
c
c     SUBROUTINE mnbrak
c
c     PURPOSE:  routine for initially bracketing a minimum
c
c               Given a function func2, and given distinct initial points
c               ax and bx, this routine searches in the downhill
c               direction (defined by the function as evaluated by
c               the initial points) and returns new points ax,bx,cx that
c               bracket a minimum of the function. Also returned are the
c               function values at the three points, fa,fb and fc.
c               parameters: gold is the default ratio by which 
c                    successive intervals are magnified; glimit
c                    is the maximum magnification allowed for a parabolic
c                    fit step.
c
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 393ff
c
c     SUBROUTINES: need subroutine func
c                    func2 : evaluates the value at point x
c
c     VARIABLES:
c 
c     IN:      ax      : initial point a
c              bx      : initial point b
c              cx      : initial point c
c              func2    : external double precision function
c  
c     OUT:     ax      : end point a
c              bx      : end point b
c              cx      : end point c
c              fa      : value of func2 at new point a
c              fb      : value of func2 at new point b
c              fc      : value of func2 at new point c
c
c     REVISIONS:
c
c     date     author   remarks
c
c     12. 8.96   ng     first implementation
c     13.08.96   ng     modified if condition to avoid explosion
c 
c =====================================================================
c
      subroutine mnbrak(ax,bx,cx,fa,fb,fc,func2)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      implicit none
      double precision ax,bx,cx,fa,fb,fc,func2,gold,glimit,tiny
      parameter (gold=1.618034d0, glimit=100.d0, tiny=1.d-20)
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      double precision dum,fu,q,r,u,ulim,tmp
c
c --------------------------------------------------------------------
c     external routines
c --------------------------------------------------------------------
c
      external func2
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--MNBRAK : trying to bracket minimum...'
c
c      write(*,*) '      values of ax,bx : ',ax,bx
c
c --------------------------------------------------------------------
c     initialization
c --------------------------------------------------------------------
c
      fa=func2(ax)
      fb=func2(bx)
c
c --------------------------------------------------------------------
c     switch role of a and b so that we can go downhill in the direction
c        from a to b
c --------------------------------------------------------------------
c
      if(fb.gt.fa)then
         dum=ax
         ax=bx
         bx=dum
         dum=fb
         fb=fa
         fa=dum
      endif
c
c --------------------------------------------------------------------
c     first guess for c
c --------------------------------------------------------------------
c
      cx=bx+gold*(bx-ax)
      fc=func2(cx)
c
c --------------------------------------------------------------------
c     for debugging
c --------------------------------------------------------------------
c
c      write(*,*) 'first guess values : ax,bx,cx : ',ax,bx,cx
c      write(*,*) '                   : fa,fb,fc : ',fa,fb,fc
c
c --------------------------------------------------------------------
c     do while: keep returning here until we bracket
c --------------------------------------------------------------------
c
 1    if ((fb.ge.fc) .and. (fb .ge. tiny) .and.
     $     ((fb .ne. fc) .and. (fb .ne. fa))) then
c 1    if (fb.ge.fc) then          ! original
c
c         write(*,*) 'ax,bx,cx : ',ax,bx,cx
c         write(*,*) 'fa,fb,fc : ',fa,fb,fc
c
c --------------------------------------------------------------------
c     calculate parabolic u
c --------------------------------------------------------------------
c
         r=(bx-ax)*(fb-fc)
         q=(bx-cx)*(fb-fa)
         u=bx-((bx-cx)*q-(bx-ax)*r)/
     $        (2.d0*sign(max(abs(q-r),tiny),q-r))
         ulim=bx+glimit*(cx-bx)
c
c --------------------------------------------------------------------
c     parabolic u is between b and c : try it
c --------------------------------------------------------------------
c
         if((bx-u)*(u-cx) .gt. 0.d0)then
c
            fu=func2(u)
c
c --------------------------------------------------------------------
c     got a minimum beteen b and c
c --------------------------------------------------------------------
c
            if (fu.lt.fc) then
               ax=bx
               fa=fb
               bx=u
               fb=fu
               return
c
c --------------------------------------------------------------------
c     got a minimum beteen a and u
c --------------------------------------------------------------------
c
            elseif (fu.gt.fb) then
               cx=u
               fc=fu
               return
            endif
c
c --------------------------------------------------------------------
c     parabolic fit was of no use. use default magnification
c --------------------------------------------------------------------
c
            u=cx+gold*(cx-bx)
            fu=func2(u)
c
c --------------------------------------------------------------------
c     parabolic fit between c and its allowed limit
c --------------------------------------------------------------------
c
         else if ((cx-u)*(u-ulim).gt. 0.d0) then
c
            fu=func2(u)
c
            if(fu.lt.fc)then
               bx=cx
               cx=u
               u=cx+gold*(cx-bx)
               fb=fc
               fc=fu
               fu=func2(u)
            endif
c
c --------------------------------------------------------------------
c     limit parabolic u to maximum allowed value
c --------------------------------------------------------------------
c
         else if((u-ulim)*(ulim-cx) .ge. 0.d0)then
            u=ulim
            fu=func2(u)
c
c --------------------------------------------------------------------
c     reject parabolic u, use default magnification
c --------------------------------------------------------------------
c
          else
c
             u=cx+gold*(cx-bx)
             fu=func2(u)
         endif
c
c --------------------------------------------------------------------
c     eliminate oldest point and continue
c --------------------------------------------------------------------
c
         ax=bx
         bx=cx
         cx=u
         fa=fb
         fb=fc
         fc=fu
         goto 1
c
      endif
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software +50.
