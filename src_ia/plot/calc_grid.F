c =====================================================================
c
c     SUBROUTINE calc_grid
c
c     PURPOSE: calculates the minimum and maximum values for an
c                grid axis and the number of tickmarks. Intended
c                to give nice labels and limits for automatically
c                scaled axis
c
c     SOURCE: partially taken from johannes sander, kup
c
c     SUBROUTINES: digits
c              
c
c     VARIABLES:
c 
c     IN:      x          : array of real values
c              ndata      : number of values in x 
c                         
c     OUT:     axismin    : minimum of axis
c              axismax    : maximum of axis
c              nmajticks  : number of major tickmarks
c              nminticks  : number of minor tickmarks
c
c     REVISIONS:
c
c     date     author   remarks
c
c     9. 4.96   ng     first implementation
c                       
c
c 
c =====================================================================
c
      subroutine calc_grid(x,ndata,
     $     axismin,axismax,nmajticks,nminticks)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      implicit none
      real x,axismin,axismax
      integer ndata,nmajticks,nminticks
      dimension x(ndata)
c
c --------------------------------------------------------------------
c     external functions
c --------------------------------------------------------------------
c
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      integer j,nmajdef
      parameter(nmajdef = 10)
      real xmin,xmax,dx,delta,fac
      real xl,xr
      integer iz,nvalid
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--CALC_GRID: calculating number of tickmarks...'
c      write(*,*) '      ndata, x : ',ndata,x
c
c --------------------------------------------------------------------
c     find minimum and maximum
c --------------------------------------------------------------------
c
      xmin = 9.99e10
      xmax = -9.99e10
      nvalid=0
c     
CJML      do j=1,ndata
CJML         if (x(j) .gt. xmax) then
CJML            xmax = x(j)
CJML         endif
CJML         if (x(j) .lt. xmin) then
CJML            xmin = x(j)
CJML         endif
CJML      enddo
      do j = 1, ndata
         if (x(j) .eq. x(j)) then        ! ignore NaNs (NaN .ne. NaN)
            if (x(j) .gt. xmax) xmax = x(j)
            if (x(j) .lt. xmin) xmin = x(j)
            nvalid = nvalid + 1
         endif
      enddo

      if (nvalid .le. 1) then
         if (nvalid .eq. 1) then
            delta = 0.5*max(1.0, abs(xmin))
            axismin = xmin - delta
            axismax = xmin + delta
         else
            axismin = 0.0
            axismax = 1.0
         endif
         nmajticks = 5
         nminticks = 4
         return
      endif
c
c --------------------------------------------------------------------
c     add 5% at each side
c      removed npg, 19.10.99
c --------------------------------------------------------------------
c
c      delta = xmax - xmin
c      if (delta .eq. 0.0) then
c         delta = 0.05 * xmin
c      endif
c      xmin = xmin - 0.05*delta
c      xmax = xmax + 0.05*delta
c
c --------------------------------------------------------------------
c     find optimum delta between xmin and xmax
c --------------------------------------------------------------------
c
CJML      dx = (xmax-xmin)/(float(nmajdef))
      delta = xmax - xmin
      if (delta .le. 0.0) then
         delta = 0.05*max(1.0, abs(xmin))
         if (delta .eq. 0.0) delta = 1.0
         xmin = xmin - delta
         xmax = xmax + delta
         delta = xmax - xmin
      endif

      dx = delta/float(nmajdef)

      if (abs(dx) .lt. 1)then
        iz = int(abs(log10(abs(dx)))+1)
      else
        iz = -1*int(log10(abs(dx)))
      endif
c
      fac = (10.**float(iz))
      dx = float(int(dx*fac+sign(0.5,dx)))/fac
      xl = float(int(xmin/dx))*dx
      xr = float(int(xmax/dx))*dx
c

      if (xl .gt. xmin) xl = xl-dx
      if (xr .lt. xmax) xr = xr+dx
c
c --------------------------------------------------------------------
c     assign values
c --------------------------------------------------------------------
c
      axismin = xl
      axismax = xr
CJML      nmajticks = nint((xr-xl)/dx)
      if (abs(dx) .lt. 1.0e-30) dx = 1.0
      nmajticks = nint((xr - xl)/dx)      
      nminticks = 5
c
c      write(*,*) 'xmin : ',xmin
c      write(*,*) 'xmax : ',xmax
c      write(*,*) 'xl,xr : ',xl,xr
c      write(*,*) 'dx : ',dx
c      write(*,*) 'nmajticks : ',nmajticks
c
c --------------------------------------------------------------------
c     return to calling program
c --------------------------------------------------------------------
c
c
      return
      end
