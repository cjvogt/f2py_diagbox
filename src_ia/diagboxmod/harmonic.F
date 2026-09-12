c =====================================================================
c
c     DOUBLE PRECISION FUNCTION HARMONIC
c 
c     PURPOSE: this function calculates the value of the function
c
c                harmonic = c1 + c2*sin(2*pi*t/365.) 
c                              + c3*cos(2*pi*t/365.)
c                              + c4*sin(4*pi*t/365.)
c                              + c5*cos(4*pi*t/365.) ...
c
c     VARIABLES:
c 
c     IN : time         : time (years)
c          c            : array holding the coefficients
c          nc           : number of coefficients
c     OUT: harmonic     : value of above equation
c 
c     REVISION:
c
c     date     author  remarks
c   
c     31.07.96   ng    first implementation
c     15.10.99   ng    changed time units to years
c
c =====================================================================
c
      double precision function harmonic(time,c,nc)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
      double precision time
      integer nc
      double precision c(nc)
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer i
      double precision sum,period,pi,i_dp,my_time
      parameter (pi = 3.141592654d0, period = 1.d0)
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--HARMONIC: at the beginning...'
c      write(*,*) '    coefficients : '
c      do i = 1,nc
c         write(*,*) 'i,c(i) : ' ,i,c(i)
c      enddo

      my_time = time - 1900.d0
      sum = c(1)
      do i = 2,nc,2
         i_dp = dble(float(i))
         sum = sum + 
     $        c(i)   * dsin(i_dp*pi*my_time/period) + 
     $        c(i+1) * dcos(i_dp*pi*my_time/period)
      enddo
c
      harmonic = sum
c      write(*,*) '--HARMONIC: at the end. : ',sum
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

