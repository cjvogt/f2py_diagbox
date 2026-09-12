c =====================================================================
c
c     FUNCTION ran1
c
c     PURPOSE:  calculates random number deviate between 0.0 and 1.0
c
c               minimal random number generator of Park and Miller with
c               Bays-Durham shuffle and added safeguards. Returns a uniform
c               random deviate between 0.0 and 1.0 (exclusive of the
c               end-point values). Call with idum a negative integer to
c               initialize; thereafter, do not alter idum between 
c               successive deviates in a sequence. 
c
c
c     SOURCE:  code from Michele Bruno, 
c               is different than that from 
c                    Numerical Recipes, 2.ed, p. 271ff
c                
c
c     SUBROUTINES: none
c
c     VARIABLES:
c 
c     IN:     idum : negative integer 
c  
c     OUT:    ran1 : random number deviate between 0.0 and 1.0
c
c     REVISIONS:
c
c     date     author   remarks
c
c     12.08.96   ng     first implementation
c 
c =====================================================================
c
      function ran1(idum)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      implicit none
c
      integer idum
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      integer ia,im,iq,ir,ntab,ndiv
      real ran1,am,eps,rnmx
      parameter (ia=16807,im=2147483647,am=1./im,iq=127773,ir=2836)
      parameter (ntab=32,ndiv=1+(im-1)/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer j,k,iv(ntab),iy
      save iv,iy
c
      data iv /ntab*0/, iy /0/
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      if ((idum.le.0).or.(iy.eq.0)) then
c         write(*,*) '--RAN1 : initializing'
         idum = max(-idum,1)
c
         do j = ntab+8,1,-1
            k = idum/iq
            idum = ia*(idum-k*iq)-ir*k
            if (idum .lt. 0) idum = idum + im
            if (j .le. ntab) iv(j) = idum
         enddo
         iy = iv(1)
      endif
c
c      write(*,*) '--RAN1 : computing...'
c
      k = idum/iq
      idum = ia*(idum-k*iq)-ir*k
      if (idum .lt. 0) idum = idum + im
      j = 1+iy/ndiv
      iy = iv(j)
      iv(j) = idum
      ran1 = min(am*iy,rnmx)
c
c      write(*,*) '--RAN1 : at end '
c      write(*,*) '       ran1  : ',ran1
c
c --------------------------------------------------------------------
c     return 
c --------------------------------------------------------------------
c
c
      return
c
      END
