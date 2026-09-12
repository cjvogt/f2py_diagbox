c =====================================================================
c
c     FUNCTION iplin
c
c     PURPOSE: calculates linear interpolated values from input
c                data. assumes that data are sorted in increasing order
c                of x
c
c     VARIABLES:
c 
c     IN:   x0     : x for which y value is calculated                     
c           x      : array of x-values, used for interpolation             
c           y      : array of y-values, used for interpolation             
c           n      : dimension of the x, xsort, y and ysort arrays          
c
c     OUT:  iplin  : interpolated value y0=y(x0)
c           error  : error = 0: no error                               
c                    error = 1: extrapolated value                          
c     CHANGES:
c
c     date     author   remarks
c
c     15.10.99   ng     adapted from standard iplin   
c 
c =====================================================================
c
      double precision function iplin(x0,x,y,n,error)
c
      implicit none
c
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      integer n
      double precision  x0, x(n), y(n)
      integer error
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      integer i, k, index
c
      error=0
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c     
c --------------------------------------------------------------------
c     evaluate for x0 the left index of x, make also
c        check for extrapolation
c --------------------------------------------------------------------
c    
      if (x(1) .gt. x0) then
         error=1
         index=1
      elseif (x(n) .lt. x0) then
         error=1
         index=n-1
      else         
         do i=1,n-1
            if (x0 .ge. x(i) .and. x0 .le. x(i+1)) then
               index=i
            endif
         enddo
      endif
c     
c --------------------------------------------------------------------
c     calculate interpolated value (or extrapolated)
c --------------------------------------------------------------------
c         
      iplin=(y(index)-y(index+1))
     $      /(x(index)-x(index+1))
      iplin=y(index)+(x0-x(index))*iplin
c     
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c        
      return
c     
      end
