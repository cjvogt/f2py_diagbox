c =====================================================================
c
c     DOUBLE PRECISION FUNCTION VP_PUREWATER
c
c     PURPOSE: this function calculates the vapor pressure of
c                pure water over the temperature range
c                273 K - 648 K represented by the Chebyshev polynomial
c                (Ambrose and Lawrenson, 1972)
c
c     REFERENCES: DOE (1994). Handbook of methods for the analysis
c                    of the various parameters of the carbon dioxide
c                    system in sea water; version 2, A.G. Dickson and
c                    C. Goyet, eds. ORNL/CDIAC-74. pages 5.22ff
c
c     SUBROUTINES: needs function chebev to evaluate Chebyshev poly-
c                  nomials
c
c     CHECK_VALUE: for 25 deg C, vp(H20) = 3.1686 kPa
c
c     VARIABLES:
c 
c     IN:   temp         : temperature [deg C]
c     OUT:  vp_purewater : vapor pressure of pure-water [kPa]
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     1.08.96   ng      first implementation 
c 
c =====================================================================
c
      double precision function vp_purewater(temp)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none

      double precision temp
c
c --------------------------------------------------------------------
c     external functions
c --------------------------------------------------------------------
c   
      double precision chebev
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      double precision a(12),rhs
      double precision tk,x
c
      data a /2794.0144d0, 1430.6181d0, -18.2465d0, 7.6875d0, 
     $          -0.0328d0,    0.2728d0,   0.1371d0, 0.0629d0,
     $           0.0261d0,    0.0200d0,   0.0117d0, 0.0067d0 /
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     evaluate chebyshev polynomial
c --------------------------------------------------------------------
c        
      tk = temp + 273.15d0
      x = (2.d0 * tk - 921.d0)/375.d0
      rhs = chebev(-1.d0,1.d0,a,12,x)
c
c --------------------------------------------------------------------
c     calculate water vapor pressure
c --------------------------------------------------------------------
c        
      vp_purewater = 10.d0**(rhs / tk)
c
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      end
c
c
c
c
c =====================================================================
c
c     DOUBLE PRECISION FUNCTION CHEBEV
c
c     PURPOSE: this function evaluates the chebyshev polynomials
c              all arguments are input. c(1:m) is an array of
c              Chebyshev coefficients. The Chebyshev polynomial
c              sum(k=1,m)(ck*Tk-1(y) - c1/2)
c              is evaluated at a point y = [x-(b+a)/2]/[(b-a)/2],
c              and the result is returned as the function value
c
c     SOURCE:  code taken from Numerical Recipes, 2.ed, p. 673ff
c
c
c     VARIABLES:
c 
c     IN:   a      : lower limit of range
c           b      : upper limit of range
c           c      : array of coefficients
c           m      : number of coefficients
c           x      : where the chebyshev polynomial should be
c                      evaluated.
c     OUT:  chebev : value of chebyshev polynomial
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     1.08.96   ng      first implementation 
c 
c =====================================================================
c
      DOUBLE PRECISION FUNCTION chebev(a,b,c,m,x)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
      integer m
      double precision a,b,x,c(m)
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      INTEGER j
      double precision d,dd,sv,y,y2
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     check range
c --------------------------------------------------------------------
c     
      if ((x-a)*(x-b).gt.0.) then
         pause 'x not in range in chebev'
         return
      endif
c
c --------------------------------------------------------------------
c     evaluate polynomial
c --------------------------------------------------------------------
c     
      d=0.
      dd=0.
      y=(2.*x-a-b)/(b-a)
      y2=2.*y
      do 11 j=m,2,-1
        sv=d
        d=y2*d-dd+c(j)
        dd=sv
11    continue
      chebev=y*d-dd+0.5*c(1)
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software +50.
