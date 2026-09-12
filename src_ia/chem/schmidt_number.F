c =====================================================================
c
c     DOUBLE PRECISION FUNCTION SCHMIDT_NUMBER
c
c     PURPOSE: this function calculates the Schmidt number versus
c                temperature for seawater (35 psu) ranging from zero
c                to 30 deg C.
c              the dependence is of the form
c
c             Sc = A - B*t + C*t**2 - D*t**3 
c
c     REFERENCES: Wanninkhof, R. 1992. Relationship between Wind speed
c           and gas-exchange over the ocean. J. Geophys. Res., 97, C5, 
c           7373-7382
c
c     CREDITS: most of the subroutine taken from Fortunat Joos
c
c     VARIABLES:
c 
c     IN:   gas_type         : character denoting type of gas
c                                 <He>  : Helium
c                                 <Ne>  : Neon
c                                 <Ar>  : Argon
c                                 <O2>  : oxygen
c                                 <CH4> : methane
c                                 <CO2> : carbon dioxide
c                                 <SF6> : sulfur hexa-fluoride
c                                 <Rn>  : Radon
c                                 <F12> : CCl2F2
c                                 <F11> : CCl3F
c
c           temp             : temperature [deg C]
c     OUT:  schmidt_number   : Schmidt-number of this gas (no units)
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     1.08.96   ng      first implementation 
c 
c =====================================================================
c
      double precision function schmidt_number(gas_type,temp)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none

      double precision temp
      character gas_type*3
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      double precision a,b,c,d
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     get coefficients according to gas_type
c --------------------------------------------------------------------
c
      if (gas_type .eq. 'He') then
c
         a = 410.14d0
         b = 20.503d0
         c = 0.53175d0
         d = 0.0060111
c
      elseif (gas_type .eq. 'Ne') then
c
         a = 855.1d0
         b = 46.299d0
         c = 1.254d0
         d = 0.01449d0
c
      elseif (gas_type .eq. 'Ar') then
c
         a = 1909.1d0
         b = 125.09d0
         c = 3.9012d0
         d = 0.048953d0
c
      elseif (gas_type .eq. 'O2') then
c
         a = 1953.4d0
         b = 128.0d0
         c = 3.9918d0
         d = 0.050091d0
c
      elseif (gas_type .eq. 'CH4') then
c
         a = 2039.2d0
         b = 120.31d0
         c = 3.4209d0
         d = 0.040437d0
c
      elseif (gas_type .eq. 'CO2') then
c
         a = 2073.1d0
         b = 125.62d0
         c = 3.6276d0
         d = 0.043219d0
c
      elseif (gas_type .eq. 'Rn') then
c
         a = 3412.8d0
         b = 224.30d0
         c = 6.7954d0
         d = 0.08300d0
c
      elseif (gas_type .eq. 'SF6') then
c
         a = 3531.6d0
         b = 231.40d0
         c = 7.2168d0
         d = 0.090558d0
c
      elseif (gas_type .eq. 'F12') then
c
         a = 3713.2d0
         b = 243.30d0
         c = 7.5879d0
         d = 0.095215d0
c
      elseif (gas_type .eq. 'F11') then
c
         a = 4039.8d0
         b = 264.70d0
         c = 8.2552d0
         d = 0.10359d0
c
      else 
c
         write(*,*)'--SCHMIDT_NUMBER: incorrect type ',
     $        'specification'
         return
c     
      endif
c
c --------------------------------------------------------------------
c     calculate schmidt - number
c --------------------------------------------------------------------
c        
      schmidt_number = a - b*temp + C*temp**2 - D*temp**3
c
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c
      return
      end
