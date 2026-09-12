c =====================================================================
c
c     FUNCTION gasdev
c
c     PURPOSE:  calculates a random gaussian deviate
c
c               returns a normally distributed deviated with zero mean
c               and unit variance, using ran1(idum) as the source of
c               uniform deviates.
c
c
c     SOURCE:  code from Numerical Recipes, 2.ed, p. 280
c                
c
c     SUBROUTINES: none
c
c     VARIABLES:
c 
c     IN:     idum   : negative integer 
c  
c     OUT:    gasdev : normally distributed random deviate
c
c     REVISIONS:
c
c     date     author   remarks
c
c     16.08.96   ng     first implementation
c 
c =====================================================================
c
      function gasdev(idum)
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
      real gasdev
      integer iset
      real  fac,gset,rsq,v1,v2
c
      save iset,gset
      data iset/0/
c
c --------------------------------------------------------------------
c     external functions
c --------------------------------------------------------------------
c
      real ran1
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--GASDEV: iset,idum : ',iset,idum
      if (iset.eq.0) then
 1       v1 = 2.0*ran1(idum)-1.0
         v2 = 2.0*ran1(idum)-1.0
         rsq = v1**2+v2**2
c
         if((rsq.ge.1.0).or.(rsq.eq.0.0)) goto 1
c
         fac = sqrt(-2.0*log(rsq)/rsq)
         gset = v1*fac
         gasdev = v2*fac
         iset = 1
      else
         gasdev = gset
         iset = 0
      endif
c      write(*,*) '--GASDEV: end gasdev = ',gasdev
c
c --------------------------------------------------------------------
c     return 
c --------------------------------------------------------------------
c
c
      return
c
      END
