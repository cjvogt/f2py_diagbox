      double precision function CO2smpl(
     $     Alk,
     $     TCO2,
     $     K2divK1,
     $     AlkCdAlk)
c     c-------------------------------------------------------------------c
c     c- Calculating CO2 concentration value:                            -c
c     c- For the iterative CO2 search an estimated value for the solution-c
c     c- has to be known, which can be calculated with this subroutine   -c
c     c- here.                                                           -c
c     c-                                                                 -c
c     c- Following simplification is used (and therefore an explicit     -c
c     c- formula can be given here):                                     -c
c     c- Total Alkalinity is set to the carbon alkalinity.               -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- Input:                                                          -c
c     c- ------                                                          -c
c     c-                                                                 -c
c     c- Alk: Total Alkalinity                                           -c
c     c-                                                                 -c
c     c- TCO2: total dissolved inorganic carbon                          -c
c     c-                                                                 -c
c     c- AlkCdAlk: Estimated ratio of carbonate to total alkalinity.     -c
c     c-           note: 0<< AlkCdAlk:=AlkC/TAlk <1                      -c
c     c-           A good value for this ratio is about 0.95.            -c
c     c-           The number of loops used for the iteration is reduced,-c
c     c-           if a good value of this ratio has been choosen.       -c
c     c-           [AlkCdAlk]=1                                          -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- output:                                                         -c
c     c- -------                                                         -c
c     c- CO2 concentration                                               -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- units:                                                          -c
c     c- ------                                                          -c
c     c- "eq" in the alkalinity term should correspond to "mol" of the   -c
c     c- total dissolved inorganic carbonate concentration.              -c
c     c-                                                                 -c
c     c- e.g.                                                            -c
c     c- [Alk]=eq/(m3 of solution) and [TCO2]=mol/(m3 of solution)       -c
c     c- or                                                              -c
c     c- [Alk]=eq/(kg of solution) and [TCO2]=mol/(kg of solution)       -c
c     c-                                                                 -c
c     c- The unit of CO2 corresponds to the unit of the input            -c
c     c- concentrations.                                                 -c
c     c-                                                                 -c
c     c-------------------------------------------------------------------c
c     
c     
      implicit none
c     
c     
c     global variables:
      double precision Alk, TCO2, K2divK1, AlkCdAlk
c
c     local variables:
      double precision AlkC, Diskr, h1, h2, CO3
c
c     
c     
      AlkC=AlkCdAlk*Alk
      h1=1.d0-4.d0*K2divK1
      h2=h1*AlkC-TCO2
      Diskr=h2*h2+4.d0*K2divK1*h1*AlkC*AlkC
      CO3=(h2+dsqrt(Diskr))/(2.d0*h1)
      CO2smpl=K2divK1*((AlkC-2.d0*CO3)**2.d0)/CO3
c
      end
c     co2smpl
c
c
c
c
c
c
c
c
      double precision function TCO2smpl(
     $     Alk,
     $     CO2,
     $     K2divK1,
     $     AlkCdAlk)
c     c-------------------------------------------------------------------c
c     c- Calculating TCO2 concentration value:                           -c
c     c- For the iterative TCO2 search an estimated value for the soluti--c
c     c- on has to be known, which can be calculated with this subroutine-c
c     c- here.                                                           -c
c     c-                                                                 -c
c     c- Following simplification is used (and therefore an explicit     -c
c     c- formula can be given here):                                     -c
c     c- Total Alkalinity is set to the carbon alkalinity.               -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- Input:                                                          -c
c     c- ------                                                          -c
c     c-                                                                 -c
c     c- Alk: Total Alkalinity                                           -c
c     c-                                                                 -c
c     c- CO2: CO2 concentration                                          -c
c     c-                                                                 -c
c     c- AlkCdAlk: Estimated ratio of carbonate to total alkalinity.     -c
c     c-           note: 0<< AlkCdAlk:=AlkC/TAlk <1                      -c
c     c-           A good value for this ratio is about 0.95.            -c
c     c-           The number of loops used for the iteration is reduced,-c
c     c-           if a good value of this ratio has been choosen.       -c
c     c-           [AlkCdAlk]=1                                          -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- output:                                                         -c
c     c- -------                                                         -c
c     c- TCO2 concentration                                              -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- units:                                                          -c
c     c- ------                                                          -c
c     c- "eq" in the alkalinity term should correspond to "mol" of the   -c
c     c- disolved co2 gas concentration.                                 -c
c     c-                                                                 -c
c     c- e.g.                                                            -c
c     c- [Alk]=eq/(m3 of solution) and [CO2]=mol/(m3 of solution)        -c
c     c- or                                                              -c
c     c- [Alk]=eq/(kg of solution) and [CO2]=mol/(kg of solution)        -c
c     c-                                                                 -c
c     c- The unit of TCO2 corresponds to the unit of the input           -c
c     c- concentrations.                                                 -c
c     c-                                                                 -c
c     c-------------------------------------------------------------------c
c     
      implicit none
c     
c     
c     global variables:
      double precision Alk, CO2, K2divK1, AlkCdAlk
c     
c     local variables:
      double precision AlkC, Diskr, hlp, CO3
c
c
c     calculate CO3=: 
      AlkC=AlkCdAlk*Alk
      hlp=CO2/(8.d0*K2divK1)
      Diskr=hlp*(AlkC+hlp)
      CO3=AlkC/2.d0+hlp-dsqrt(Diskr)

c     calculate TCO2:
      TCO2smpl=AlkC-CO3+CO2
c
      end
c     tco2smpl
c
c
c
c
c
c
      double precision function CO3smpl(
     $     Alk,
     $     TCO2,
     $     K2divK1,
     $     AlkCdAlk)
c     c-------------------------------------------------------------------c
c     c- Calculating CO3 concentration value:                            -c
c     c- For the iterative CO3 search an estimated value for the solution-c
c     c- has to be known, which can be calculated with this subroutine   -c
c     c- here.                                                           -c
c     c-                                                                 -c
c     c- Following simplification is used (and therefore an explicit     -c
c     c- formula can be given here):                                     -c
c     c- Total Alkalinity is set to the carbon alkalinity.               -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- Input:                                                          -c
c     c- ------                                                          -c
c     c-                                                                 -c
c     c- Alk: Total Alkalinity                                           -c
c     c-                                                                 -c
c     c- TCO2: total dissolved inorganic carbon                          -c
c     c-                                                                 -c
c     c- AlkCdAlk: Estimated ratio of carbonate to total alkalinity.     -c
c     c-           note: 0<< AlkCdAlk:=AlkC/TAlk <1                      -c
c     c-           A good value for this ratio is about 0.95.            -c
c     c-           The number of loops used for the iteration is reduced,-c
c     c-           if a good value of this ratio has been choosen.       -c
c     c-           [AlkCdAlk]=1                                          -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- output:                                                         -c
c     c- -------                                                         -c
c     c- CO3 concentration                                               -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- units:                                                          -c
c     c- ------                                                          -c
c     c- "eq" in the alkalinity term should correspond to "mol" of the   -c
c     c- total dissolved inorganic carbonate concentration.              -c
c     c-                                                                 -c
c     c- e.g.                                                            -c
c     c- [Alk]=eq/(m3 of solution) and [TCO2]=mol/(m3 of solution)       -c
c     c- or                                                              -c
c     c- [Alk]=eq/(kg of solution) and [TCO2]=mol/(kg of solution)       -c
c     c-                                                                 -c
c     c- The unit of CO3 corresponds to the unit of the input            -c
c     c- concentrations.                                                 -c
c     c-                                                                 -c
c     c-------------------------------------------------------------------c
c     
c     
      implicit none
c     
c     
c     global variables:
      double precision Alk, TCO2, K2divK1, AlkCdAlk
c
c     local variables:
      double precision AlkC, Diskr, h1, h2
c
c     
c     
      AlkC=AlkCdAlk*Alk
      h1=1.d0-4.d0*K2divK1
      h2=h1*AlkC-TCO2
      Diskr=h2*h2+4.d0*K2divK1*h1*AlkC*AlkC
      CO3smpl=(h2+dsqrt(Diskr))/(2.d0*h1)
c
      end
c     co3smpl
c
c
c
      double precision function TAsmpl(
     $     TCO2,
     $     CO2,
     $     K2divK1,
     $     AlkCdAlk)
c     c-------------------------------------------------------------------c
c     c- Calculating TAlk concentration value:                           -c
c     c- For the iterative TAlk search an estimated value for the        -c
c     c- solution                                                        -c
c     c- has to be known, which can be calculated with this subroutine   -c
c     c- here.                                                           -c
c     c-                                                                 -c
c     c- Following simplification is used (and therefore an explicit     -c
c     c- formula can be given here):                                     -c
c     c- Total Alkalinity is set to the carbon alkalinity.               -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- Input:                                                          -c
c     c- ------                                                          -c
c     c-                                                                 -c
c     c- TCO2: total dissolved inorganic carbon                          -c
c     c-                                                                 -c
c     c- CO2: CO2 concentration                                          -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- AlkCdAlk: Estimated ratio of carbonate to total alkalinity.     -c
c     c-           note: 0<< AlkCdAlk:=AlkC/TAlk <1                      -c
c     c-           A good value for this ratio is about 0.95.            -c
c     c-           The number of loops used for the iteration is reduced,-c
c     c-           if a good value of this ratio has been choosen.       -c
c     c-           [AlkCdAlk]=1                                          -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- output:                                                         -c
c     c- -------                                                         -c
c     c- TAlk concentration                                              -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- units:                                                          -c
c     c- ------                                                          -c
c     c- "eq" in the alkalinity term should correspond to "mol" of the   -c
c     c- total dissolved inorganic carbonate concentration.              -c
c     c-                                                                 -c
c     c- e.g.                                                            -c
c     c- [Alk]=eq/(m3 of solution) and [TCO2]=mol/(m3 of solution)       -c
c     c- or                                                              -c
c     c- [Alk]=eq/(kg of solution) and [TCO2]=mol/(kg of solution)       -c
c     c-                                                                 -c
c     c- The unit of TAlk corresponds to the unit of the input           -c
c     c- concentrations.                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- revisions:                                                      -c
c     c- ---------                                                       -c
c     c-                                                                 -c
c     c- this part added by Nicolas Gruber, 12.7.96                      -c
c     c-                                                                 -c
c     c-------------------------------------------------------------------c
c     
c     
      implicit none
c     
c     
c     global variables:
      double precision CO2, TCO2, K2divK1, AlkCdAlk
c
c     local variables:
      double precision diskr, a,b,c,co3,AlkC
c
c 
      a = 1.d0
      b = CO2 - 2.d0*TCO2 - CO2/K2divK1
      c = TCO2**2 + CO2**2 - 2*TCO2*CO2
c
      co3 = (-b - dsqrt(b**2 - 4*a*c)) / (2*a)
c     
      AlkC = dsqrt(co3*co2/K2divK1) + 2.d0*co3
c
      TAsmpl = AlkC / AlkCdAlk
c
c ---------------------------------------------------------------------
c     for debugging purposes
c ---------------------------------------------------------------------
c
c      write(*,*) '--TAsmpl : debugging...'
c      write(*,*) 'a,b,c : ',a,b,c
c      write(*,*) 'co3 : ',co3
c      write(*,*) 'alkC : ',alkC
c      write(*,*) 'TAsmpl : ',TAsmpl
      end

c     TAsmpl
c
c
c
