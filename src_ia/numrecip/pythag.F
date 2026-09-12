c =====================================================================
c
c     FUNCTION pythag
c
c     PURPOSE:  computes (a**2 + b**2)^(1/2) without destructive
c                 underflow or overflow
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 62ff
c
c
c     VARIABLES:
c 
c     IN:      a       : array of x-values
c              b       : array of y-values
c  
c     OUT:     pythag  : (a**2 + b**2)^(1/2)
c
c     REVISIONS:
c
c     date     author   remarks
c
c     29. 7.96   ng     first implementation
c     15. 8.96   ng     changed real to double precision
c 
c =====================================================================
c
      double precision function pythag(a,b)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      double precision a,b
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      double precision absa,absb
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      absa=abs(a)
      absb=abs(b)
      if(absa.gt.absb)then
         pythag=absa*sqrt(1.d0+(absb/absa)**2)
      else
         if(absb.eq.0.d0)then
            pythag=0.d0
         else
            pythag=absb*sqrt(1.d0+(absa/absb)**2)
         endif
      endif
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
       return
      end
C  (C) Copr. 1986-92 Numerical Recipes Software +50.
