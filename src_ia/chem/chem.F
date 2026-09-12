c--------------------------------------------------------------------------c
c-                                                                        -c
c- Author :  Roger Fink,                                                  -c
c-           Physics Institute,                                           -c
c-           Climate and Environmental Physics,                           -c
c-           CH-3012 Bern                                                 -c
c-                                                                        -c
c-           Nicolas Gruber                                               -c
c-           AOS Program, Princeton University                            -c
c-                                                                        -c
c-           addition of k1_dm(tc,s) and k2_dm(tc,s)                      -c
c-                by Nicolas Gruber (1996/07/26)                          -c
c-           addition of virial coefficient for pure CO2                  -c
c-                       virial coefficient for CO2 in air                -c
c-                       fco2 (calc of fco2 from xco2                     -c
c-           addition of Mehrbach et al. (1973) constants                 -c
c                        k1_mb(tc,s) and k2_mb(tc,s)                      -c
c-           addition of Mucci calcite and aragonite solubilities         -c
c-                  kspcalc_mu(tc,s,dp) and ksparag_mu(tc,s,dp)           -c
c-           addition of new pressure dependencies, taken from pmolvol.r  -c
c-                  given to me by Chris Sabine (March 1999)              -c
c-                                                                        -c
c- Date   :  1996/07/16                                                   -c
c-           revised 1996/07/26                                           -c
c-           revised 1998/03/03                                           -c
c-           revised 1999/04/01                                           -c
c-           revised 1999/04/19                                           -c
c- Version:  2.02                                                         -c
c-           2.03                                                         -c
c-           2.04                                                         -c
c-           2.05                                                         -c
c-                                                                        -c
c--------------------------------------------------------------------------c
c
c     This file contains functions, which are used in order
c     to compute chemical properties of seawater. 
c
c
c     double precision function solco2m3(tc,s)
c     (CO2 solubility after Weiss referring to one m3 of solution)
c     
c     double precision function solco2kg(tc,s)
c     (CO2 solubility after Weiss referring to one kg of solution)
c
c     double precision function vir_coeff_pure_co2(tc)  (after Weiss, 74)
c
c     double precision function vir_coeff_co2_air(tc)   (after Weiss, 74)
c
c     double precision function fco2(tc,p,xco2)  (after Handbook)
c     
c     double precision function solo2(tc,s)      (O2  solubility after Weiss)
c
c     double precision function k1_COM(tc,s)
c                                         (Millero 95, COMbined Roy93, G&P89)
c     double precision function k2_COM(tc,s)
c                                         (Millero 95, COMbined Roy93, G&P89)
c     double precision function k1_gp(tc,s)      (K1 of Goyet & Poisson)
c     double precision function k2_gp(tc,s)      (K2 of Goyet & Poisson)
c
c     double precision function k1_dm(tc,s) (Dickson and Millero, 1987)
c     double precision function k2_dm(tc,s) (Dickson and Millero, 1987)
c
c     double precision function k1_mb(tc,s) (Mehrbach et al., 1973) 
c     double precision function k2_mb(tc,s) (Mehrbach et al., 1973) 
c
c     double precision function kb_dicks(tc,s)
c
c     double precision function Kw_m95drcph(tc,s)
c     double precision function kw_DiRi(tc,s)
c
c     double precision function k1si_ym95(tc,s)
c     double precision function k1si_sni
c
c     function kp1_ym95(tc,s)
c     function kp2_ym95(tc,s)
c     function kp3_ym95(tc,s)
c     double precision function kp1_DiRi(tc,s)
c     double precision function kp2_DiRi(tc,s)
c     double precision function kp3_DiRi(tc,s)
c
c     double precision tb_uppst(s)
c
c     double precision function ca2(s)           (Ca2+ ion concentration)
c
c     double precision function kspcalc(tc,s,dp) (solubility product of calcit)
c     double precision function kspcalcFB(tc,s,dp)
c
c     double precision funciton kspcalc_mu(tc,s,dp) solub. product after
c     double precision function ksparag_mu(tc,s,dp)   Mucci, 1983
c     
c     double precision function
c     co3sat(tc,s,dp,ca2plus)     (saturation concentration of CaCO3)
c     
c     double precision function
c     getsatho(tc,s,co3,ca2plus,density,g) 
c                                       (depth of saturation horizon of caco3)
c
c     double precision function putponk1(dp,tc,s) ("put press on disso k1")
c     double precision function putponk2(dp,tc,s) ("put press on disso k2")
c     double precision function putponkbOLD(dp,tc)("put press on disso kb")
c     double precision function putponkb(dp,tc,s) ("put press on disso kw")
c     double precision function putponkw(dp,tc,s) ("put press on disso kw")
c
c     double precision function putponk1_mi(dp,tc,s) ("put press on disso k1")
c     double precision function putponk2_mi(dp,tc,s) ("put press on disso k2")
c     double precision function putponkb_mi(dp,tc,s) ("put press on disso kw")
c     double precision function putponkw_mi(dp,tc,s) ("put press on disso kw")
c
c
c----------------------------------------------------------------------------c
c
c
c
c
      double precision function solco2m3(tc,s)
c--------------------------------------------------------c
c     input:                                            -c
c     ------                                            -c
c     tc: temperature [tc] = oC                         -c
c     s : salinity [s] = o/oo                           -c
c                                                       -c
c                                                       -c
c     output:                                           -c
c     --------                                          -c
c     - co2 solubility after Weiss 1974, c.f. UNESCO 87 -c
c                                                       -c
c     - Units:                                          -c
c     - mol/(m3*atm) (referring to one m3 of solution)  -c
c--------------------------------------------------------c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc, s
c
c
c     constants:
c     ----------
      double precision a1,a2,a3,b1,b2,b3
      parameter(a1=-58.0931d0, a2=90.5069d0, a3=22.2940d0)
      parameter(b1=0.027766d0, b2=-0.025888d0, b3=0.0050578d0)
c
c     Zur Umrechnung von (dm)**(-3) auf m**(-3):
      double precision factor
      parameter(factor=1000.0d0)
c
c
c     local variables:
c     ----------------
      double precision h
c
c
      h=(tc+2.7315d2)/1.d2
      solco2m3=factor*dexp(a1+a2/h+a3*dlog(h)+s*(b1+(b2+b3*h)*h))
      end
c
c end of function "solco2m3"
c
c
c
c
      double precision function solco2kg(tc,s)
c--------------------------------------------------------c
c     input:                                            -c
c     ------                                            -c
c     tc: temperature [tc] = oC                         -c
c     s : salinity [s] = o/oo                           -c
c                                                       -c
c                                                       -c
c     output:                                           -c
c     --------                                          -c
c     - co2 solubility after Weiss 1974, c.f. UNESCO 87 -c
c     - K0 = [CO2*]/fCO2                                -c
c                                                       -c
c     - Units:                                          -c
c     - mol/(kg*atm) (referring to one kg of solution)  -c
c                                                       -c
c     check values:                                     -c
c     -----------                                       -c
c     T = 298.15 K; S = 35 permil;                      -c
c     ln(solco2kg/k0) = -3.5617 (k0 = 1 mol kg soln-1)  -c  
c--------------------------------------------------------c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc, s
c
c
c     constants:
c     ----------
      double precision a1,a2,a3,b1,b2,b3
      parameter(a1=-60.2409d0, a2=93.4517d0, a3=23.3585d0)
      parameter(b1=0.023517d0, b2=-0.023656d0, b3=0.0047036d0)
c
c
c     local variables:
c     ----------------
      double precision h
c
c
      h=(tc+2.7315d2)/1.d2
      solco2kg=dexp(a1+a2/h+a3*dlog(h)+s*(b1+(b2+b3*h)*h))
c      write(*,*) 'ln(K0/k0) = ',dlog(solco2kg)
      end
c
c end of function "solco2kg"
c
c
c
c
      double precision function vir_coeff_pure_co2(tc)
c--------------------------------------------------------c
c     input:                                            -c
c     ------                                            -c
c     tc: temperature [tc] = oC                         -c
c                                                       -c
c                                                       -c
c     output:                                           -c
c     --------                                          -c
c     - first virial coefficient of pure CO2 B(T)       -c
c         after Weiss, 1974                             -c
c                                                       -c
c     - Units:                                          -c
c     - cm3/mol                                         -c
c                                                       -c
c     references:                                       -c
c     -----------                                       -c
c                                                       -c
c     R. Weiss, 1974. Carbon dioxide in water and sea-  -c
c       water: the solubility of a non-ideal gas. Mar.  -c
c       Chem., 2, 203-215                               -c
c                                                       -c
c                                                       -c
c     check values:                                     -c
c     -----------                                       -c
c     T = 298.15 K; B(T) = -123.2 cm3/mol               -c
c                                                       -c
c--------------------------------------------------------c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc
c
c
c     constants:
c     ----------
      double precision a1,a2,a3,a4
      parameter(a1 = -1636.75d0, a2 = 12.0408d0)
      parameter(a3 = -3.27957d-2, a4 = 3.16528d-5)
c
c
c     local variables:
c     ----------------
      double precision h
c
c
      h=(tc+2.7315d2)
      vir_coeff_pure_co2 = a1 + a2*h + a3*h**2 + a4*h**3
      end
c
c end of function "vir_coeff_pure_co2"
c
c
c
c
      double precision function vir_coeff_co2_air(tc)
c--------------------------------------------------------c
c     input:                                            -c
c     ------                                            -c
c     tc: temperature [tc] = oC                         -c
c                                                       -c
c                                                       -c
c     output:                                           -c
c     --------                                          -c
c     - first virial coefficient of CO2 in air          -c
c         after Weiss, 1974                             -c
c                                                       -c
c     - Units:                                          -c
c     - cm3/mol                                         -c
c                                                       -c
c     references:                                       -c
c     -----------                                       -c
c                                                       -c
c     R. Weiss, 1974. Carbon dioxide in water and sea-  -c
c       water: the solubility of a non-ideal gas. Mar.  -c
c       Chem., 2, 203-215                               -c
c                                                       -c
c                                                       -c
c     check values:                                     -c
c     -----------                                       -c
c     T = 298.15 K; del(T) = 22.5 cm3/mol               -c
c                                                       -c
c--------------------------------------------------------c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc
c
c
c     constants:
c     ----------
      double precision a1,a2
      parameter(a1 = 57.7d0, a2 = -0.118d0)
c
c
c     local variables:
c     ----------------
      double precision h
c
c
      h=(tc+2.7315d2)
      vir_coeff_co2_air = a1 + a2*h 
      end
c
c end of function "vir_coeff_co2_air"
cc
c
c
c
c
      double precision function fco2(tc,p,xco2)
c---------------------------------------------------------c
c     input:                                             -c
c     ------                                             -c
c     tc   : temperature [tc] = oC                       -c
c     p    : total pressure [atm]                        -c
c     xco2 : mole fraction of CO2 in air [ppm]           -c
c                                                        -c
c                                                        -c
c     output:                                            -c
c     --------                                           -c
c     - fugacity of CO2 in air                           -c
c                                                        -c
c     - Units:                                           -c
c     - uatm                                             -c
c                                                        -c
c     references:                                        -c
c     -----------                                        -c
c                                                        -c
c     R. Weiss, 1974. Carbon dioxide in water and sea-   -c
c       water: the solubility of a non-ideal gas. Mar.   -c
c       Chem., 2, 203-215                                -c
c     DOE, 1995. Handbook of methods for the analysis    -c
c        of the various parameters of the carbon dioxide -c
c        system in sea water; version 2. A.G. Dickson    -c
c        and C. Goyet, eds. ORNL/CDIAC-74.               -c
c                                                        -c
c     special:                                           -c
c     --------                                           -c
c                                                        -c
c     need functions vir_coeff_pure_co2 and              -c
c                    vir_coeff_co2_air                   -c
c                                                        -c
c     check values:                                      -c
c     -----------                                        -c
c     T= 298.15 K, p = 1atm, xCO2 = 350 ppm,             -c
c        fco2 = 348.9 uatm                               -c
c                                                        -c
c---------------------------------------------------------c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc,p,xco2
c
c
c     constants:
c     ----------
c
c     external functions:
c     -------------------
      double precision vir_coeff_pure_co2,vir_coeff_co2_air
c
c
c     local variables:
c     ----------------
      double precision tk,bt,del,pp,xco2_1,fco2_tmp
c
c
      tk=(tc+2.7315d2)          ! temperature Kelvin
      pp = p * 101325d0         ! pressure in Pa
      xco2_1 = xco2 * 1.d-6
c
      bt = vir_coeff_pure_co2(tc)
      del = vir_coeff_co2_air(tc)
      
      fco2_tmp = xco2_1 * pp * 
     $     dexp((bt*1.d-6+2*(1.d0-xco2_1)**2 * del * 1.d-6)*pp/
     $     (8.31451d0*tk))      ! units Pa
c
      fco2 = fco2_tmp / .101325d0 ! units uatm
      end
c
c end of function "fco2"
cc
c
c
c
c
c
c
      double precision function solo2(tc,s)
c-----------------------------------------------------------------------------c
c     - R.F. Weiss, Deep Sea Research                                        -c
c       Vol. 17. pp. 721 to 735                                              -c
c                                                                            -c
c     "air solubility"                                                       -c
c                                                                            -c
c     - Units:                                                               -c
c     --------                                                               -c
c     tc   : temperature in degree celcius                                   -c
c     s    : salinity in o/oo                                                -c
c     solo2: m3(STP) of o2 / m3 of seawater for a total atmospheric pressure -c
c            of 1 atm.                                                       -c
c                                                                            -c
c     rem.:                                                                  -c
c     The effect of water vapore pressure and the mole                       -c
c     fraction of o2 (about 0.209) is taken into account                     -c
c     in the "air solubility". It is related to a total                      -c
c     atmospheric air-pressure of one atmosphere.                            -c
c-----------------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision  tc, s
c
c
c     constants:
c     ----------
      double precision a1,a2,a3,a4,b1,b2,b3
      parameter(a1=-173.4292d0, a2=249.6339d0,
     $          a3=143.3483d0, a4=-21.8492d0)
      parameter(b1=-3.3096d-2, b2=1.4259d-2, b3=-1.7d-3)
c
c     unit conversion factor:
      double precision lperml
      parameter(lperml=1.d-3)
c
c     local variables:
c     ----------------
      double precision h
c
c
      h=(tc+2.7315d2)/1.d2
      solo2=lperml*dexp(a1+a2/h+a3*dlog(h)+a4*h+s*(b1+(b2+b3*h)*h))
      end
c
c end of function "solo2"
c
c
c
c
      double precision function k1_COM(tc,s)
c----------------------------------------------------------------------------c
c- Thermodynamics of the dioxide system in the oceans                       -c
c- Frank J. Millero                                                         -c
c- Goechimica et Cosmochimica Acta, Vol. 59, No 4. pp. 661-677, 1995        -c
c-                                                                          -c
c- Millero (1995) has COMbined the results of                               -c
c-                                                                          -c
c- 1) The dissociation constants of carbonic acid in seawater at salinities -c
c-    5 to 45 and temperatures 0 to 45 oC                                   -c
c-    Rabindra N. Roy et al.                                                -c
c-    Marine Chemistry, 44 (1993) 249-267                                   -c
c-                                                                          -c
c- with                                                                     -c
c-                                                                          -c
c- 2) Deep-Sea Research, Vol. 36, N0.11, pp. 1635-1654, 1989                -c
c-    New determination of carbonic acid dissociation constants in seawater -c
c-    a function of temperature and salinity                                -c
c-    Catherine Goyet and Alain Poisson                                     -c
c-    SWS pH-scale including fluorid                                        -c
c-    range of temperature: -1 C - 40 C                                     -c
c-    range of salinity   : 10 - 50                                         -c
c-                                                                          -c
c-                                                                          -c
c- These constants are related to the pH(SWS) scale.                        -c
c- Concentration unit: mol/(kg of solution)                                 -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc,s
c
c
c     constants:
c     ----------
      double precision a,b,c,d,e,f,g
      parameter(
     $     a=2.18867d0,
     $     b=-2.2750360d3, c=-1.468591d0,
     $     d=-0.138681d0,  e=-9.33291d0,
     $     f=7.26483d-2,  g=-5.74938d-3)
c
c     
c     local variables:
c     ----------------
      double precision tk
c
c
      tk=2.7315d2+tc
c
      k1_COM=
     $     dexp(
     $     a+b/tk+c*dlog(tk)+(d+e/tk)*dsqrt(s)+f*s+g*s**(1.5d0)
     $     )
c
      end
c end of function "k1_COM"
c
c
c
      double precision function k2_COM(tc,s)
c----------------------------------------------------------------------------c
c- Thermodynamics of the dioxide system in the oceans                       -c
c- Frank J. Millero                                                         -c
c- Goechimica et Cosmochimica Acta, Vol. 59, No 4. pp. 661-677, 1995        -c
c-                                                                          -c
c- Millero (1995) has COMbined the results of                               -c
c-                                                                          -c
c- 1) The dissociation constants of carbonic acid in seawater at salinities -c
c-    5 to 45 and temperatures 0 to 45 oC                                   -c
c-    Rabindra N. Roy et al.                                                -c
c-    Marine Chemistry, 44 (1993) 249-267                                   -c
c-                                                                          -c
c- with                                                                     -c
c-                                                                          -c
c- 2) Deep-Sea Research, Vol. 36, N0.11, pp. 1635-1654, 1989                -c
c-    New determination of carbonic acid dissociation constants in seawater -c
c-    a function of temperature and salinity                                -c
c-    Catherine Goyet and Alain Poisson                                     -c
c-    SWS pH-scale including fluorid                                        -c
c-    range of temperature: -1 C - 40 C                                     -c
c-    range of salinity   : 10 - 50                                         -c
c-                                                                          -c
c-                                                                          -c
c- These constants are related to the pH(SWS) scale.                        -c
c- Concentration unit: mol/(kg of solution)                                 -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc,s
c
c
c     constants:
c     ----------
      double precision a,b,c,d,e,f,g
      parameter(
     $     a=-0.84226d0,
     $     b=-3.7411288d3, c=-1.437139d0,
     $     d=-0.128417d0,  e=-24.41239d0,
     $     f=0.1195308d0,  g=-9.12840d-3)
c     
c     local variables:
c     ----------------
      double precision tk
c
c
      tk=2.7315d2+tc
c
      k2_COM=
     $     dexp(
     $     a+b/tk+c*dlog(tk)+(d+e/tk)*dsqrt(s)+f*s+g*s**(1.5d0)
     $     )
c
      end
c end of function "k2_COM"
c
c
c
c
c
c
      double precision function k1_gp(tc,s)
c--------------------------------------------------------------------------c
c     input:                                                              -c
c     ------                                                              -c
c     tc: temperature [tc] = oC                                           -c
c     s : salinity [s] = o/oo                                             -c
c                                                                         -c
c                                                                         -c
c                                                                         -c
c     output:                                                             -c
c     -------                                                             -c
c     - Deep-Sea Research, Vol. 36, N0.11, pp. 1635-1654, 1989            -c
c     - New determination of carbonic acid dissociation constants in      -c
c                   seawater                                              -c
c     - a function of temperature and salinity                            -c
c     - Catherine Goyet and Alain Poisson                                 -c
c     - SWS pH-scale including fluorid                                    -c
c     - range of temperature: -1 degree C to 40 degree C                  -c
c     - range of salinity   : 10 - 50                                     -c
c                                                                         -c
c     - Output units:                                                     -c
c       [k1_gp]=mol/(kg of solution)                                      -c
c--------------------------------------------------------------------------c
c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc,s
c
c
c     constants:
c     ----------
      double precision a,b,c,d
      parameter(a=8.1227d2, b=3.356d0, c=-1.71d-3, d=9.1d-5)
c
c     
c     local variables:
c     ----------------
      double precision pk1, tk
c
c
      tk=2.7315d2+tc
      pk1=a/tk+b+c*S*dlog(tk)+d*s*s
      k1_gp=1.d1**(-pk1)
c
      end
c end of function "k1_gp"
c
c
c
c
c
      double precision function k2_gp(tc,s)
c---------------------------------------------------------------------------c
c     input:                                                               -c
c     ------                                                               -c
c     tc: temperature [tc] = oC                                            -c
c     s : salinity [s] = o/oo                                              -c
c                                                                          -c
c                                                                          -c
c                                                                          -c
c     output:                                                              -c
c     -------                                                              -c
c     - Deep-Sea Research, Vol. 36, N0.11, pp. 1635-1654, 1989             -c
c     - New determination of carbonic acid dissociation constants in       -c
c                    seawater                                              -c
c     - a function of temperature and salinity                             -c
c     - Catherine Goyet and Alain Poisson                                  -c
c     - SWS pH-scale including fluorid                                     -c
c     - range of temperature: -1 degree C to 40 degree C                   -c
c     - range of salinity   : 10 - 50                                      -c
c                                                                          -c
c     - Output units:                                                      -c
c       [k1_gp]=mol/(kg of solution)                                       -c
c---------------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc,s
c
c
c     constants:
c     ----------
      double precision a,b,c,d
      parameter(a=1450.87d0, b=4.604d0, c=-3.85d-3, d=1.82d-4)
c
c     
c     local variables:
c     ----------------
      double precision pk2, tk
c
c
      tk=2.7315d2+tc
      pk2=a/tk+b+c*s*dlog(tk)+d*s*s
      k2_gp=1.d1**(-pk2)
c
      end
c end of function "k2_gp"
c
c
c
c
c
c
      double precision function k1_dm(tc,s)
c--------------------------------------------------------------------------c
c     input:                                                              -c
c     ------                                                              -c
c     tc: temperature [tc] = oC                                           -c
c     s : salinity [s] = o/oo                                             -c
c                                                                         -c
c     output:                                                             -c
c     -------                                                             -c
c     - Deep-Sea Research, Vol. 34, N0.10, pp. 1733-1743, 1987            -c
c     - A comparison of the equilibrium constants for the dissociation    -c
c            of carbonin acid in seawater media                           -c
c     - A.G. Dickson and F.J. Millero                                     -c
c     - SWS pH-scale including fluorid                                    -c
c     - range of temperature: 2 degree C to 35 degree C                   -c
c     - range of salinity   : 20 - 40                                     -c
c                                                                         -c
c     - Output units:                                                     -c
c       [k1_dm]=mol/(kg of solution)                                      -c
c--------------------------------------------------------------------------c
c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc,s
c
c
c     constants:
c     ----------
      double precision a,b,c,d
      parameter(a=845.0d0, b=3.248d0, c=-0.0098d0, d=0.000087d0)
c
c     
c     local variables:
c     ----------------
      double precision pk1, tk
c
c
      tk=2.7315d2+tc
      pk1=a/tk+b+c*S+d*s*s
      k1_dm=1.d1**(-pk1)
c
      end
c end of function "k1_dm"
c
c
c
c
c
c
      double precision function k2_dm(tc,s)
c--------------------------------------------------------------------------c
c     input:                                                              -c
c     ------                                                              -c
c     tc: temperature [tc] = oC                                           -c
c     s : salinity [s] = o/oo                                             -c
c                                                                         -c
c     output:                                                             -c
c     -------                                                             -c
c     - Deep-Sea Research, Vol. 34, N0.10, pp. 1733-1743, 1987            -c
c     - A comparison of the equilibrium constants for the dissociation    -c
c            of carbonin acid in seawater media                           -c
c     - A.G. Dickson and F.J. Millero                                     -c
c     - SWS pH-scale including fluorid                                    -c
c     - range of temperature: 2 degree C to 35 degree C                   -c
c     - range of salinity   : 20 - 40                                     -c
c                                                                         -c
c     - Output units:                                                     -c
c       [k2_dm]=mol/(kg of solution)                                      -c
c--------------------------------------------------------------------------c
c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc,s
c
c
c     constants:
c     ----------
      double precision a,b,c,d
      parameter(a=1377.3d0, b=4.824d0, c=-0.0185d0, d=0.000122d0)
c
c     
c     local variables:
c     ----------------
      double precision pk2, tk
c
c
      tk=2.7315d2+tc
      pk2=a/tk+b+c*S+d*s*s
      k2_dm=1.d1**(-pk2)
c
      end
c end of function "k2_dm"
c
c
c
c
c
c
c
c
      double precision function k1_mb(tc,s)
c--------------------------------------------------------------------------c
c     input:                                                              -c
c     ------                                                              -c
c     tc: temperature [tc] = oC                                           -c
c     s : salinity [s] = o/oo                                             -c
c                                                                         -c
c     output:                                                             -c
c     -------                                                             -c
c     - Limnology and Oceanography, 18, 897-907                           -c
c     - Measurement of the apparent dissociation constants of carbonic    -c
c       acid in seawater at atmospheric pressure                          -c
c     - C. Mehrbach, C.H. Culberson, J.E. Hawley, and R.M. Pytkowicz      -c
c     - formula given in Millero, 1995 is used here (fitted original      -c
c         data of Mehrbach et al. (1973) to seawater pH scale             -c
c     - SWS pH-scale including fluorid                                    -c
c     - range of temperature: 2 degree C to 35 degree C                   -c
c     - range of salinity   : 20 - 40                                     -c
c                                                                         -c
c     - Output units:                                                     -c
c       [k1_mb]=mol/(kg of solution)                                      -c
c--------------------------------------------------------------------------c
c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc,s
c
c
c     constants:
c     ----------
      double precision a,b,c,d,e
      parameter(a=3670.7d0, b=-62.008d0, c=9.7944d0, d=-0.0118d0,
     $     e=0.000116d0)
c
c     
c     local variables:
c     ----------------
      double precision pk1, tk
c
c
      tk=2.7315d2+tc
      pk1 = a/tk + b + c*dlog(tk) + d*S + e*S*S
      k1_mb=1.d1**(-pk1)
c
      end
c end of function "k1_mb"
c
c
c
c
      double precision function k2_mb(tc,s)
c--------------------------------------------------------------------------c
c     input:                                                              -c
c     ------                                                              -c
c     tc: temperature [tc] = oC                                           -c
c     s : salinity [s] = o/oo                                             -c
c                                                                         -c
c     output:                                                             -c
c     -------                                                             -c
c     - Limnology and Oceanography, 18, 897-907                           -c
c     - Measurement of the apparent dissociation constants of carbonic    -c
c       acid in seawater at atmospheric pressure                          -c
c     - C. Mehrbach, C.H. Culberson, J.E. Hawley, and R.M. Pytkowicz      -c
c     - formula given in Millero, 1995 is used here (fitted original      -c
c         data of Mehrbach et al. (1973) to seawater pH scale             -c
c     - SWS pH-scale including fluorid                                    -c
c     - range of temperature: 2 degree C to 35 degree C                   -c
c     - range of salinity   : 20 - 40                                     -c
c                                                                         -c
c     - Output units:                                                     -c
c       [k2_mb]=mol/(kg of solution)                                      -c
c--------------------------------------------------------------------------c
c
c
      implicit none
c
c     global variables:
c     -----------------
      double precision tc,s
c
c
c     constants:
c     ----------
      double precision a,b,c,d
      parameter(a=1394.7d0, b=4.777d0, c=-0.0184d0, d=0.000118d0)
c
c     
c     local variables:
c     ----------------
      double precision pk2, tk
c
c
      tk=2.7315d2+tc
      pk2 = a/tk + b + c*S + d*S*S
      k2_mb=1.d1**(-pk2)
c
      end
c end of function "k2_mb"
c
c
c
c
c
c
      double precision function kb_dicks(tc,s)
c     c-----------------------------------------------------------------c
c     c- input:                                                        -c
c     c- ------                                                        -c
c     c- tc: temperature [tc] = oC                                     -c
c     c- s : salinity [s] = o/oo                                       -c
c     c-                                                               -c
c     c-                                                               -c
c     c-                                                               -c
c     c- output:                                                       -c
c     c- -------                                                       -c
c     c- [kb_dicks] = mol/(kg of solution)                             -c
c     c-                                                               -c
c     c-                                                               -c
c     c- Dissociation constant of boric acid from A. G. Dickson:       -c
c     c- =======================================================       -c
c     c-                                                               -c
c     c- Deep-Sea Research, Vol. 37, No.5, pp. 755-766, 1990           -c
c     c- Thermodynamics of the dissociation of boric acid              -c
c     c- in synthetic seawater from 273.15 to 318.15 K                 -c
c     c- Andrew G. Dickson                                             -c
c     c- SWS pH-scale without fluorid                                  -c
c     c- range of temperature: 0 oC - 45 oC                            -c
c     c- range of salinity   : 5 - 45                                  -c
c     c-                                                               -c
c     c-                                                               -c
c     c- This constant is corrected for the SWS pH scale, using        -c
c     c- the beta terms for fluor and sulphur:                         -c
c     c- =======================================================       -c
c     c-                                                               -c
c     c-                                                               -c
c     c-                                                               -c
c     ------------------------------------------------------------------c
c
      implicit none
c
c     global varaibles:
c     -----------------
      double precision tc,s
c
c     local variables:
c     ----------------
      double precision tk,
     $     BHSO4,
     $     BHF,
     $     mTSO4,mTSO4IM,
     $     mTF,
     $     mINatural,
     $     mIArtific
c
c
c
c
c     temperature in Kelvin:
c     ----------------------
      tk=273.15d0+tc
c
c
c
c     dissociation constant for artificial seawater (without fluor):
c     --------------------------------------------------------------
      kb_dicks=
     $     dexp
     $     (
     $     (
     $     -8966.90d0-2890.53d0*dsqrt(s)-77.942d0*s+
     $     1.728d0*(s**1.5d0)-0.0996d0*s*s
     $     )
     $     /tk
     $     +148.0248d0+137.1942d0*dsqrt(s)+1.62142d0*s
     $     -(24.4344d0+25.085d0*dsqrt(s)+0.2474d0*s)*dlog(tk)
     $     +0.053105d0*dsqrt(s)*tk
     $     )
c
c
c
c
c     correction for SWS:
c     -------------------
c
c     total concentrations:
      mTSO4=8.069d-4*s/(1.d0-1.005d-3*s)
      mTF=2.0d-6*s/(1.d0-1.005d-3*s)
c
c     Ionic strength of natural seawater:
      mINatural=1.99273d-2*s/(1.d0-1.005d-3*s)
c
c     Ionic strength of artificial seawater:
      mIArtific=(19.919d0*s)/(1000.d0-1.00198d0*s)
c
c     total sulphate concentration in artificial seawater:
      mTSO4IM=(2.927d-2/0.7225d0)*mIArtific
c
c     beta-terms:
      BHSO4=
     $     10.d0**
     $     (
     $     647.59d0/tk-6.3451d0+1.9085d-2*tk-5.208d-1*dsqrt(mINatural)
     $     )
      BHF=
     $     dexp
     $     (
     $     -1590.2d0/tk+12.641d0-1.525d0*dsqrt(mINatural)
     $     )
c
c     correction:
      kb_dicks=kb_dicks*(1+BHSO4*mTSO4+BHF*mTF)/(1+BHSO4*mTSO4IM)
c
      end
c     end of function kb_dicks
c
c
c
c
c
c
c
      double precision function Kw_m95drcph(tc,s)
c----------------------------------------------------------------------------c
c- Thermodynamics of the dioxide system in the oceans                       -c
c- Frank J. Millero                                                         -c
c- Goechimica et Cosmochimica Acta, Vol. 59, No 4. pp. 661-677, 1995        -c
c-                                                                          -c
c-                                                                          -c
c- Input: Tc: Temperature in degree celcius                                 -c
c-        S : Salinity, dimensionless                                       -c
c-                                                                          -c
c- Output: Unit of Kw: mol/(kg of solution)                                 -c
c-                                                                          -c
c----------------------------------------------------------------------------c
c-                                              -c
c- equation of Millero (1995), fitted to        -c
c- data of:                                     -c
c-                                              -c
c-                                              -c
c-                                              -c
c- 1) A. G. Dickson and J. P. Riley             -c
c-    Marine Chemistry, 7(1979) 89-99           -c
c-                                              -c
c-    The estimation of acid dissociation       -c
c-    constants in seawater media from          -c
c-    potentiometric titrations with strong base-c
c-    I. The ionic product of water-Kw          -c
c-                                              -c
c-    SWS pH-scale  (including fluorid)         -c
c-                                              -c
c-    Range of temperature:  5 C - 35 C  (?)    -c
c-    Range of salinity   :  20 - 40  (?)       -c
c-                                              -c
c-                                              -c
c-                                              -c
c- 2) Culberson and Pytcowicz (1973),           -c
c-                                              -c
c-    SWS pH-scale (without fluorid?)           -c
c-                                              -c
c-    Range of temperature:  0 C - 50 C         -c
c-    Range of salinity   :  0 - 40             -c
c-                                              -c
c-                                              -c
c- 3) Geochimica et Cosmochimica Acta           -c
c-    Vol. 43, pp. 1651 to 1661                 -c
c-                                              -c
c-    Hansson (1972),                           -c
c-                                              -c
c-    SWS pH-scale (without fluorid?)           -c
c-                                              -c
c-    Range of temperature:  0øC - 50øC   -c
c-    Range of salinity   :  0 - 40             -c
c-                                              -c
c-                                              -c
c-                                              -c
c------------------------------------------------c
      implicit none
c
c     global variables:
      double precision tc,s
c
c     local variable:
      double precision tk, hlp
c
      tk=273.15d0+tc
c
      hlp=dlog(tk)
      Kw_m95drcph=
     $     dexp(
     $     148.9802d0-13847.26d0/tk-23.6521d0*hlp
     $     -1.615d-2*s+dsqrt(s)*(-5.977d0+118.67d0/tk+1.0495d0*hlp)
     $     )
c
      end
c
c
c
c
      double precision function kw_DiRi(tc,s)
c     c---------------------------------------------c
c     c- Marine Chemistry, 7(1979) 89-99           -c
c     c-                                           -c
c     c- The estimation of acid dissociation       -c
c     c- constants in seawater media from          -c
c     c- potentiometric titrations with strong base-c
c     c- I. The ionic product of water-Kw          -c
c     c- A. G. Dickson and J. P. Riley             -c
c     c-                                           -c
c     c- SWS pH-scale  (including fluorid)         -c
c     c-                                           -c
c     c- Range of temperature:  5oC - 35oC  (?)    -c
c     c- Range of salinity   :  20 - 40  (?)       -c
c     c---------------------------------------------c
c
      implicit none
c
c     global variables:
      double precision tc,s
c
c     local variable:
      double precision tk
c
      tk=273.15d0+tc
      kw_DiRi=10.d0**(-3441.d0/tk-2.241d0+0.09415d0*dsqrt(s))
c
      end
c
c
c
c
      double precision function k1si_ym95(tc,s)
c----------------------------------------------------------------------------c
c-                                                                          -c
c-                                                                          -c
c- Equation found in                                                        -c
c-                                                                          -c
c- Thermodynamics of the dioxide system in the oceans                       -c
c- Frank J. Millero                                                         -c
c- Goechimica et Cosmochimica Acta, Vol. 59, No 4. pp. 661-677, 1995        -c
c-                                                                          -c
c-                                                                          -c
c- Equation from                                                            -c
c-                                                                          -c
c- Yao W. and Millero F.J. (1995)                                           -c
c- The chemistry of the anoxic waters in the Framvaren Fjord, Norway.       -c
c- J. Aquatic Geochem. (in press)                                           -c
c-                                                                          -c
c- fitted to data of                                                        -c
c- Busey and Mesmer (1977)                                                  -c
c- Sjoebberg et al. (1981)                                                  -c
c- Baes and Mesmer (1976)                                                   -c
c-                                                                          -c
c-                                                                          -c
c- Input:                                                                   -c
c- Tc: Temperature in degree celcius                                        -c
c- IonicStrength: Ionic strength in units of mol/(kg of solution)           -c
c-                                                                          -c
c-                                                                          -c
c- Output:                                                                  -c
c- The unit of KSi as given in Millero 1995 is Molality. But in this        -c
c- function the unit is converted to mol/(kg of solution).                  -c
c-                                                                          -c
c-                                                                          -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
      double precision tc,s
c
c     constants:
      double precision a,b,c,d,e,f,g,h,i
      parameter(
     $     a= 117.40d0,
     $     b=-8.9042d3,  c=-19.334d0,
     $     d= 3.5913d0,  e=-458.79d0,
     $     f=-1.5998d0,  g= 188.74d0,
     $     h= 0.07871d0, i=-12.1652d0)
c
c     local variables:
      double precision tk, ionicstrength

      tk=273.15d0+tc
      ionicstrength=1.d0/(1.d0-1.005d-3*s)*1.99273d-2*s
c
      k1si_ym95=(1.d0-1.005d-3*s)*
     $     dexp(
     $     a+b/Tk+c*dlog(Tk)
     $     +(d+e/Tk)*dsqrt(IonicStrength)
     $     +(f+g/Tk)*IonicStrength
     $     +(h+i/Tk)*IonicStrength**2.d0
     $     )
      end
c
c
c

      double precision function k1si_sni()
c     c-----------------------------------------------------------c
c     c- input:                                                  -c
c     c- ------                                                  -c
c     c- tc: temperature [tc] = oC                               -c
c     c- s : salinity [s] = o/oo                                 -c
c     c-                                                         -c
c     c-                                                         -c
c     c-                                                         -c
c     c- output:                                                 -c
c     c- -------                                                 -c
c     c- [k1si_sni] = mol/(kg of solution)                       -c
c     c-                                                         -c
c     c-                                                         -c
c     c- Marine Chemistry, 10 (1981) 521-532                     -c
c     c- Equilibrium and structural studies of silicon(IV) and   -c
c     c- aluminium(III) in aqueous solution                      -c
c     c- S. Sjoeberg, A. Nordin, N. Ingri (-->SNI)               -c
c     c-                                                         -c
c     c- Corrected for seawater (vgl. Rick)                      -c
c     c-----------------------------------------------------------c
c
      implicit none
c
      k1si_sni=3.753d-10
c
      end
c
c
c
c
c
c
c
      double precision function kp1_ym95(tc,s)
c----------------------------------------------------------------------------c
c-                                                                          -c
c-                                                                          -c
c- Equation found in                                                        -c
c-                                                                          -c
c- Thermodynamics of the dioxide system in the oceans                       -c
c- Frank J. Millero                                                         -c
c- Goechimica et Cosmochimica Acta, Vol. 59, No 4. pp. 661-677, 1995        -c
c-                                                                          -c
c-                                                                          -c
c- Equation from                                                            -c
c-                                                                          -c
c- Yao W. and Millero F.J. (1995)                                           -c
c- The chemistry of the anoxic waters in the Framvaren Fjord, Norway.       -c
c- J. Aquatic Geochem. (in press)                                           -c
c-                                                                          -c
c- fitted to data of                                                        -c
c-                                                                          -c
c-                                                                          -c
c- 1) Johansson and Wedborg, 1979                                           -c
c-                                                                          -c
c- 2) Marine Chemistry, 7(1979) 101-109                                     -c
c-    The estimation of acid dissociation constants in seawater media       -c
c-    from potentiometric titrations with strong base                       -c
c-    II. The dissociation of phosphoric acid                               -c
c-    A.G. Dickson and J.P. Riley                                           -c
c-                                                                          -c
c- 3) Kester and Pytkowicz (1967)                                           -c
c-                                                                          -c
c-                                                                          -c
c-                                                                          -c
c-                                                                          -c
c- Input:                                                                   -c
c- Tc: Temperature in degree celcius                                        -c
c- S: Salinity                                                              -c
c-                                                                          -c
c-                                                                          -c
c- Output:                                                                  -c
c- K: mol/(kg of solution)                                                  -c
c-                                                                          -c
c-                                                                          -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
      double precision tc,s
c
c     constants:
      double precision a,b,c,d,e,f,g
      parameter(
     $     a= 115.54d0,
     $     b=-4576.752d0, c=-18.453d0,
     $     d= 0.69171d0,  e=-106.736d0,
     $     f=-0.01844d0,  g=-0.65643d0)
c
c     local variables:
      double precision tk
c     
      tk=273.15d0+tc
c
      kp1_ym95=
     $     dexp(
     $     a+b/tk+c*dlog(tk)
     $     +(d+e/tk)*dsqrt(s)
     $     +(f+g/tk)*s
     $     )
c
      end
c
c
c
c
      double precision function kp2_ym95(tc,s)
c-----------------------------------------------------------------------c
c-                                                                     -c
c-                                                                     -c
c- Equation found in                                                   -c
c-                                                                     -c
c- Thermodynamics of the dioxide system in the oceans                  -c
c- Frank J. Millero                                                    -c
c- Goechimica et Cosmochimica Acta, Vol. 59, No 4. pp. 661-677, 1995   -c
c-                                                                     -c
c-                                                                     -c
c- Equation from                                                       -c
c-                                                                     -c
c- Yao W. and Millero F.J. (1995)                                      -c
c- The chemistry of the anoxic waters in the Framvaren Fjord, Norway.  -c
c- J. Aquatic Geochem. (in press)                                      -c
c-                                                                     -c
c- fitted to data of                                                   -c
c-                                                                     -c
c-                                                                     -c
c- 1) Johansson and Wedborg, 1979                                      -c
c-                                                                     -c
c- 2) Marine Chemistry, 7(1979) 101-109                                -c
c-    The estimation of acid dissociation constants in seawater media  -c
c-    from potentiometric titrations with strong base                  -c
c-    II. The dissociation of phosphoric acid                          -c
c-    A.G. Dickson and J.P. Riley                                      -c
c-                                                                     -c
c- 3) Kester and Pytkowicz (1967)                                      -c
c-                                                                     -c
c-                                                                     -c
c-                                                                     -c
c-                                                                     -c
c- Input:                                                              -c
c- Tc: Temperature in degree celcius                                   -c
c- S: Salinity                                                         -c
c-                                                                     -c
c-                                                                     -c
c- Output:                                                             -c
c- K: mol/(kg of solution)                                             -c
c-                                                                     -c
c-                                                                     -c
c-----------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
      double precision tc,s
c
c     constants:
      double precision a,b,c,d,e,f,g
      parameter(
     $     a= 172.1033d0,
     $     b=-8814.715d0, c=-27.927d0,
     $     d= 1.3566d0, e=-160.340d0,
     $     f=-0.05778d0, g= 0.37335d0)
c
c     local variables:
      double precision tk
c     
      tk=273.15d0+tc
c
      kp2_ym95=
     $     dexp(
     $     a+b/tk+c*dlog(tk)
     $     +(d+e/tk)*dsqrt(s)
     $     +(f+g/tk)*s
     $     )
c
      end
c
c
c
c
      double precision function kp3_ym95(tc,s)
c----------------------------------------------------------------------------c
c-                                                                          -c
c-                                                                          -c
c- Equation found in                                                        -c
c-                                                                          -c
c- Thermodynamics of the dioxide system in the oceans                       -c
c- Frank J. Millero                                                         -c
c- Goechimica et Cosmochimica Acta, Vol. 59, No 4. pp. 661-677, 1995        -c
c-                                                                          -c
c-                                                                          -c
c- Equation from                                                            -c
c-                                                                          -c
c- Yao W. and Millero F.J. (1995)                                           -c
c- The chemistry of the anoxic waters in the Framvaren Fjord, Norway.       -c
c- J. Aquatic Geochem. (in press)                                           -c
c-                                                                          -c
c- fitted to data of                                                        -c
c-                                                                          -c
c-                                                                          -c
c- 1) Johansson and Wedborg, 1979                                           -c
c-                                                                          -c
c- 2) Marine Chemistry, 7(1979) 101-109                                     -c
c-    The estimation of acid dissociation constants in seawater media       -c
c-    from potentiometric titrations with strong base                       -c
c-    II. The dissociation of phosphoric acid                               -c
c-    A.G. Dickson and J.P. Riley                                           -c
c-                                                                          -c
c- 3) Kester and Pytkowicz (1967)                                           -c
c-                                                                          -c
c-                                                                          -c
c-                                                                          -c
c-                                                                          -c
c- Input:                                                                   -c
c- Tc: Temperature in degree celcius                                        -c
c- S: Salinity                                                              -c
c-                                                                          -c
c-                                                                          -c
c- Output:                                                                  -c
c- K: mol/(kg of solution)                                                  -c
c-                                                                          -c
c-                                                                          -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
      double precision tc,s
c
c     constants:
      double precision a,b,d,e,f,g
      parameter(
     $     a=-18.126d0,  b=-3070.75d0,
     $     d= 2.81197d0, e= 17.27039d0,
     $     f=-0.09984d0, g=-44.99486d0)
c
c     local variables:
      double precision tk
c     
      tk=273.15d0+tc
c
      kp3_ym95=
     $     dexp(
     $     a+b/tk
     $     +(d+e/tk)*dsqrt(s)
     $     +(f+g/tk)*s
     $     )
c
      end
c
c
c
c
c
c
c
      double precision function kp1_DiRi(tc,s)
c     c-------------------------------------------------------------------c
c     c- input:                                                          -c
c     c- ------                                                          -c
c     c- tc: temperature [tc] = oC                                       -c
c     c- s : salinity [s] = o/oo                                         -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- output:                                                         -c
c     c- -------                                                         -c
c     c- [k1p_DiRi] = mol/(kg of solution)                               -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- Marine Chemistry, 7(1979) 101-109                               -c
c     c- The estimation of acid dissociation constants in seawater media -c
c     c- from potentiometric titrations with strong base                 -c
c     c- II. The dissociation of phosphoric acid                         -c
c     c- A.G. Dickson and J.P. Riley                                     -c
c     c-------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
      double precision tc,s
c     
c     local variables:
      double precision tk, INatural
c     
      tk=273.15d0+tc
      INatural=1.99273d-2*s
c
      kp1_DiRi=10.d0**(75.d0/tk-2.16d0+0.35d0*dsqrt(INatural))
c
      end
c
c
c
c
      double precision function kp2_DiRi(tc,s)
c     c-------------------------------------------------------------------c
c     c- input:                                                          -c
c     c- ------                                                          -c
c     c- tc: temperature [tc] = oC                                       -c
c     c- s : salinity [s] = o/oo                                         -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- output:                                                         -c
c     c- -------                                                         -c
c     c- [k2p_DiRi] = mol/(kg of solution)                               -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- Marine Chemistry, 7(1979) 101-109                               -c
c     c- The estimation of acid dissociation constants in seawater media -c
c     c- from potentiometric titrations with strong base                 -c
c     c- II. The dissociation of phosphoric acid                         -c
c     c- A.G. Dickson and J.P. Riley                                     -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- remark:                                                         -c
c     c- for the calculation of the Ionic strength of natural seawater   -c
c     c- the formula given in                                            -c
c     c-                                                                 -c
c     c-    Millero (1974), Frank J: Seawater as a multicomponent        -c
c     c-    electrolyte solution, p.18-20, in: The Sea, E.D. Goldberg,   -c
c     c-    Vol. 5, Marine Chemistry, Wiley & Sons 1974                  -c
c     c-                                                                 -c
c     c- is used.                                                        -c
c     c-                                                                 -c
c     c-------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
      double precision tc,s
c     
c     local variables:
      double precision tk, INatural
c     
      tk=273.15d0+tc
      INatural=1.99273d-2*s
c
      kp2_DiRi=10.d0**(-737.6d0/tk-4.196d0+0.851d0*dsqrt(INatural))
c
      end



      double precision function kp3_DiRi(tc,s)
c     c-------------------------------------------------------------------c
c     c- input:                                                          -c
c     c- ------                                                          -c
c     c- tc: temperature [tc] = oC                                       -c
c     c- s : salinity [s] = o/oo                                         -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- output:                                                         -c
c     c- -------                                                         -c
c     c- [k1p_DiRi] = mol/(kg of solution)                               -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- Marine Chemistry, 7(1979) 101-109                               -c
c     c- The estimation of acid dissociation constants in seawater media -c
c     c- from potentiometric titrations with strong base                 -c
c     c- II. The dissociation of phosphoric acid                         -c
c     c- A.G. Dickson and J.P. Riley                                     -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c-                                                                 -c
c     c- remark:                                                         -c
c     c- for the calculation of the Ionic strength of natural seawater   -c
c     c- the formula given in                                            -c
c     c-                                                                 -c
c     c-    Millero (1974), Frank J: Seawater as a multicomponent        -c
c     c-    electrolyte solution, p.18-20, in: The Sea, E.D. Goldberg,   -c
c     c-    Vol. 5, Marine Chemistry, Wiley & Sons 1974                  -c
c     c-                                                                 -c
c     c- is used.                                                        -c
c     c-                                                                 -c
c     c-------------------------------------------------------------------c
c
      implicit none
c
c     global variables:
      double precision tc,s
c     
c     local variables:
      double precision tk, INatural
c     
      tk=273.15d0+tc
      INatural=1.99273d-2*s
c
      kp3_DiRi=10.d0**(-2404.0d0/Tk-1.31d0+0.87d0*dsqrt(INatural))
c
      end
c
c
c
c
c
      double precision function tb_uppst(s)
c     c--------------------------------------------------------------c
c     c- Units: [S]=o/oo, [Total Bor TB]= mol/kg of seawater        -c
c     c-                                                            -c
c     c- Uppstrom L. R., 1974, The boron/chlorinity ratio of deep-  -c
c     c- sea water from the Pacific Ocean,                          -c
c     c- Deep-Sea Research, 21, 161-162.                            -c
c     c--------------------------------------------------------------c
c
      implicit none
c
c     global variable:
      double precision s
c
      tb_uppst = 1.188d-5*s
c
      end
c
c
c
c
c
      double precision function ca2(s)
c---------------------------------------------------------------c
c- Ca2+ - Konzentration in Abhaengigkeit der Salinitaet        -c
c- Gefunden in:                                                -c
c- The thermodynamics of the carbonate system in seawater,     -c
c- Frank j. Millero, 1979, Geochim Acta Vol 43, p.1660         -c
c-                                                             -c
c- Originalreferenz:                                           -c
c- Millero 1979, Thermodynamics of seawater, In Ocean Handbook -c
c-                                                             -c
c- [s] = o/oo                                                  -c
c- [ca2] = mol/(kg of solution)                                -c
c---------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision s
c
      ca2=2.938d-4*s
c
      end
c end of function "ca2"
c
c
      double precision function kspcalc_mu(tc,s,dp)
c-----------------------------------------------------------------------c
c     -                                                                -c
c     Solubility product of calcite,                                   -c
c     depending of temperature, salinity and pressure                  -c
c-----------------------------------------------------------------------c
c     -                                                                -c
c     -                                                                -c
c     - Input:                                                         -c
c     -       tc: Temp. in grad C,                                     -c
c     -       s: in o/oo,                                              -c
c     -       dp in dbar, pressure difference to surface pressure      -c
c     -            of 1atm                                             -c
c     -                                                                -c
c     - Output:                                                        -c
c     - solubility product of calcite in umol^2/(kg of solution)^2      -c
c     -                                                                -c
c     - Mucci, A. (1983) The solubility of calcite and aragonite       -c
c                  in seawater at various salinities, temperatures     -c
c                  and one atmosphere total pressure, Americ. J. Sci   -c
c                  283, 780-799.
c     the pressure dependence is from 
c 
c
c     the code is taken from Chris Sabine routine carbconst.r
c      April, 1999
c-----------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, dp
c
c     
c     local variables:
c     ----------------
      double precision lnksp0cal
      double precision t,pbar,saltrat
      double precision kcal
      double precision dvf,dkf,dv,dk, rt
c     
c     constants:
c     ----------
c
      double precision dvfresh(3),dvsalt(3),dkfresh(3),dksalt(3)
      data dvfresh / -65.28,    0.397,     -5.155d-3/
      data dvsalt  / -45.464,   0.3529,    -4.985d-3/
      data dkfresh /  18.47d-3, 0.1956d-3, -2.212d-6/
      data dksalt  / -13.70d-3, 0.1245d-3,  0.0d-6/

c     Initialization:
c     ---------------

      pbar=dp/10.0d0            !pbar = pressure in bar
      rt=(273.15d0+tc)*83.145d0 !RT in cm^3 bar mole^-1
      saltrat=dsqrt(s)/dsqrt(35.0d0)
      t=2.7315d2+tc

c
c ksp0cal = gamma*[Ca]gamma*[CO3]/[CaCO3 calcite]
c
      lnksp0cal = -395.8293d0 + 6537.773d0/t + 71.595d0 *
     $     dlog(t) - .17959d0 * t
c
c kcal = [Ca][CO3]/[CaCO3 calcite]
c
      kcal = dexp(lnksp0cal +
     $     (-1.78938d0 + 410.64d0/t + .0065453d0 * t) * dsqrt(s) -
     $     .17755d0 * s + .0094979d0 * s**1.5d0)
      
c
c include pressure dependence
c
c
      if (dp .gt. 0.) then
         dvf = dvfresh(1) + dvfresh(2)*tc + dvfresh(3)*tc*tc
         dkf = dkfresh(1) + dkfresh(2)*tc + dkfresh(3)*tc*tc
         dv = dvf + 
     $        (dvsalt(1) + dvsalt(2)*tc + dvsalt(3)*tc*tc - dvf) 
     $        * saltrat
         dk = dkf + 
     $        (dksalt(1) + dksalt(2)*tc + dksalt(3)*tc*tc - dkf) 
     $        * saltrat
         kcal = kcal*dexp(-dv*pbar/rt + 0.5d0*dk*pbar*pbar/rt)
      endif

      kspcalc_mu = kcal
      end
c end of function "kspcalc_mu"
c
c
      double precision function ksparag_mu(tc,s,dp)
c-----------------------------------------------------------------------c
c     -                                                                -c
c     Solubility product of aragonite,                                 -c
c     depending of temperature, salinity and pressure                  -c
c-----------------------------------------------------------------------c
c     -                                                                -c
c     -                                                                -c
c     - Input:                                                         -c
c     -       tc: Temp. in grad C,                                     -c
c     -       s: in o/oo,                                              -c
c     -       dp in dbar, pressure difference to surface pressure      -c
c     -            of 1atm                                             -c
c     -                                                                -c
c     - Output:                                                        -c
c     - solubility product of calcite in umol^2/(kg of solution)^2      -c
c     -                                                                -c
c     - Mucci, A. (1983) The solubility of calcite and aragonite       -c
c                  in seawater at various salinities, temperatures     -c
c                  and one atmosphere total pressure, Americ. J. Sci   -c
c                  283, 780-799.
c     the pressure dependence is from 
c 
c
c     the code is taken from Chris Sabine routine carbconst.r
c      April, 1999
c-----------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, dp
c
c     
c     local variables:
c     ----------------
      double precision lnksp0arag
      double precision t,pbar,saltrat
      double precision karag
      double precision dvf,dkf,dv,dk, rt
c     
c     constants:
c     ----------
c
      double precision dvfresh(3),dvsalt(3),dkfresh(3),dksalt(3)
      data dvfresh /-62.50,   0.397,     -5.155d-3/
      data dvsalt /-42.680,   0.3529,    -4.985d-3/
      data dkfresh /18.47d-3, 0.1956d-3, -2.212d-6/
      data dksalt /-13.70d-3, 0.1245d-3,  0.0d-6/

c     Initialization:
c     ---------------

      pbar=dp/10.0d0            !pbar = pressure in bar
      rt=(273.15d0+tc)*83.145d0 !RT in cm^3 bar mole^-1
      saltrat=dsqrt(s)/dsqrt(35.0d0)
      t=2.7315d2+tc

c
c ksp0arag = gamma*[Ca]gamma*[CO3]/[CaCO3 aragonite]
c
      lnksp0arag = -395.9180d0 + 6685.079d0/t + 71.595d0 *
     $     dlog(t) - .17959d0 * t
c
c karag = [Ca][CO3]/[CaCO3 aragonite]
c
      karag = dexp(lnksp0arag +
     $     (-.157481d0 + 202.938d0/t + .0039780d0 * t) * dsqrt(s) -
     $     .23067d0 * s + .0136808d0 * s**1.5d0)
c
c include pressure dependence
c
      if (dp .gt. 0.) then
         dvf = dvfresh(1) + dvfresh(2)*tc + dvfresh(3)*tc*tc
         dkf = dkfresh(1) + dkfresh(2)*tc + dkfresh(3)*tc*tc
         dv = dvf + 
     $        (dvsalt(1) + dvsalt(2)*tc + dvsalt(3)*tc*tc - dvf) 
     $        * saltrat
         dk = dkf + 
     $        (dksalt(1) + dksalt(2)*tc + dksalt(3)*tc*tc - dkf) 
     $        * saltrat
         karag = karag*dexp(-dv*pbar/rt + 0.5d0*dk*pbar**2.0d0/rt)
      endif

      ksparag_mu = karag
      end
c end of function "ksparag_mu"
c
c
c
      double precision function kspcalc(tc,s,dp)
c-----------------------------------------------------------------------c
c     -                                                                -c
c     Solubility product of calcit,                                    -c
c     depending of temperature, salinity and pressure                  -c
c-----------------------------------------------------------------------c
c     -                                                                -c
c     -                                                                -c
c     - Input:                                                         -c
c     -       tc: Temp. in grad C, s in o/oo,                          -c
c     -       dp in Pa                                                 -c
c     -                                                                -c
c     - dp: Druckdifferenz (Pa) zu 1 atm                               -c
c     -     (evtl. zu 1 bar. In Literatur mies beschrieben.)           -c
c     -     1 atm = 1.01325 10**5 Pa, 1 bar = 10**5 Pa                 -c
c     -                                                                -c
c     - Zu den Einheiten:                                              -c
c     - [tc]= grad C, [s]=o/oo,                                        -c
c     - [dp] = Pa,                                                     -c
c     - [dv]=(cm)**3/mol; [dk]=(cm)**3/(mol*atm)                       -c
c     - Avogadro constant r: [r]=atm*(cm)**3/(mol*K), somit:           -c
c     - r=8.31441 J/(mol*K) = 82.0568 atm (cm)**3/(mol*K)              -c
c     -                                                                -c
c     - Output:                                                        -c
c     - Loeslichkeitsprodukt von CaCO3 in mol^2/(kg of solution)^2     -c
c     -                                                                -c
c     - Referenz zur Chemie vgl.                                       -c
c     - Thermodynamics of the carbon dioxide system in seawater,       -c
c     - UNESCO technical papers in marine science, 51, 1987, p. 14-17  -c
c-----------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, dp
c
c     constants:
c     ----------
      double precision b0, b1, b2, c0, d0
      parameter(b0=-7.7712d-1, b1=2.8426d-3, b2=1.7834d2)
      parameter(c0=-7.711d-2, d0=4.1249d-3)
c
c     
c     local variables:
c     ----------------
      double precision dpinatm
      double precision t, ln10
      double precision logksp0c, logkspsc
      double precision dvc, dkc, rt
c
c
c     Initialisierungen:
c     ------------------
c     Umwandlung von Pa in atm
      dpinatm=dp/1.01325d5 
      t=2.7315d2+tc
      ln10=dlog(1.d1)
      rt=8.20568d1*t
c
c     - 10-er log des "calcit solubility products" berechnen: 
      logksp0c=-171.9065d0-0.077993d0*t+2839.319d0/t+71.595d0*dlog10(t)
      logkspsc=logksp0c+(b0+b1*t+b2/t)*dsqrt(s)+c0*s+d0*s**1.5d0
c
c     - Druckabhaengigkeit berechnen:
      dvc=5.304d-1*tc-4.876d1
      dkc=3.692d-4*tc-1.176d-2
      kspcalc=dexp(ln10*logkspsc)*
     $     dexp( (-dvc/rt)*dpinatm+(0.5d0*dkc/rt)*dpinatm*dpinatm )
c
c
      end
c end of function "kspcalc"
c
c
c
      double precision function kspcalcFB(tc,s,dp)
c-----------------------------------------------------------------------c
c     -                                                                -c
c     Solubility product of calcit,                                    -c
c     depending of temperature, salinity and pressure                  -c
c-----------------------------------------------------------------------c
c code adappted from Frank Basinot
C	(mostly from equations and constants given in Morse and McKenzie, 1990)
c-----------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, dp
c
c     constants:
c     ----------
      double precision b0, b1, b2, c0, d0
      parameter(b0=-7.7712d-1, b1=2.8426d-3, b2=1.7834d2)
      parameter(c0=-7.711d-2, d0=4.1249d-3)
c
c     
c     local variables:
c     ----------------
      double precision dpinatm
      double precision t
c     , ln10
c      double precision logksp0c, logksp0a, logkspsc
      double precision dvc, dkc, rt
      double precision kcp0
c
c
c     Initialisierungen:
c     ------------------
c     Umwandlung von Pa in atm
      dpinatm=dp/1.01325d5 
      t=2.7315d2+tc
c      ln10=dlog(1.d1)
      rt=8.20568d1*t
c
c
c
C	(mostly from equations and constants given in Morse and McKenzie, 1990)
c     
c [d ]  mco3sat : concentrattion of  co3-- at saturation [nmol/cm3]
      Kcp0= dexp(303.1308d0-(13348.09d0/t)-(48.7537d0*DLOG(t))+
     $     (DSQRT(s)*(1.6233d0-(118.64d0/t)))-(0.06999d0*s))
c     Konversion of units [nmol/cm3] into [mol/m3]:
      Kcp0=1.d-3*Kcp0
c
c
c
c
c     - 10-er log des "calcit solubility products" berechnen: 
c      logksp0a=-1.71945d2-7.7993d-2*t+2.903293d3/t+7.1595d1*dlog(t)/ln10
c      logksp0c=logksp0a-3.85d-2-6.3974d1/t
c      logkspsc=logksp0c+(b0+b1*t+b2/t)*dsqrt(s)+c0*s+d0*s**1.5d0
c
c     - Druckabhaengigkeit berechnen:
      dvc=5.304d-1*tc-4.876d1
      dkc=3.692d-4*tc-1.176d-2
      kspcalcFB=Kcp0*
     $     dexp( (-dvc/rt)*dpinatm+(0.5d0*dkc/rt)*dpinatm*dpinatm )
c
      end
c end of function "kspcalcFB"
c
c
c
c
c
      double precision function co3sat(tc,s,dp,ca2plus)
c-----------------------------------------------------------------------c
c     - Input:                                                         -c
c     -       tc: Temp. in grad C                                      -c
c     -       s : salinity in o/oo,                                    -c
c     -       dp in Pa                                                 -c
c     -       ca2plus in mol/(kg of solution)                          -c
c     -                                                                -c
c     - Bem.:                                                          -c
c     - dp: Druckdifferenz (Pa) zu 1 atm (evtl. 1 bar)                 -c
c     -     1 atm = 1.01325 10**5 Pa, 1 bar = 10**5 Pa                 -c
c     -                                                                -c
c     - Zu den Einheiten:                                              -c
c     - [tc]= grad C, [s]=o/oo,                                        -c
c     - [dp] = Pa,                                                     -c
c     - [ca2plus]=mol/(kg of solution)                                 -c
c     - [dv]=(cm)**3/mol                                               -c
c     - [dk]=(cm)**3/(mol*atm)                                         -c
c     - [r]=atm*(cm)**3/(mol*K), somit:                                -c
c     - r=8.31441 J/(mol*K) = 82.0568 atm (cm)**3/(mol*K)              -c
c     -                                                                -c
c     - Output:                                                        -c
c     - Saetigungskonzentration des CO3= in mol/(kg of solution)       -c
c     -                                                                -c
c     - Referenz zur Chemie vgl.                                       -c
c     - Thermodynamics of the carbon dioxide system in seawater,       -c
c     - UNESCO technical papers in marine science, 51, 1987, p. 14-17  -c
c-----------------------------------------------------------------------c
c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, dp, ca2plus
c
c     external function:
      double precision kspcalc
      external kspcalc
c
c     - Saettigungskonzentration des CO3= berechnen:
      co3sat=kspcalc(tc,s,dp)/ca2plus
c
      end
c end of function "co3sat"
c
c
c
c
c
c
      double precision function getsatho(tc,s,co3,ca2plus,density,g)
c----------------------------------------------------------------------c
c     - Berechnung der Tiefe der Lysocline:                           -c
c     -                                                               -c
c     - Input:  [tc]= grad C                                          -c
c     -         [s]=o/oo                                              -c
c     -         [co3]=[ca2+]=mol/(kg of solution)                     -c
c     -         [density]=kg/m^3                                      -c
c     -                                                               -c
c     - Output: [lysdepth]=m                                          -c
c     -                                                               -c
c     - Naeherung: Dichte des Wassers wird als konstant angenommen.   -c
c----------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, co3, ca2plus, density, g
c
c
c     external function:
c     ------------------
      double precision kspcalc
      external kspcalc
c
c     local variables:
c     ----------------
      double precision a, b, c, diskr
      double precision t, rt, dvc, dkc, dp
      double precision ksp0c, kspsc
c
c
c     Loeslichkeitsprodukt von CaCO3 an der Wasseroberflaeche berechnen:
c     --------------------------------------------------------------------
      ksp0c=kspcalc(tc,s,0.d0)
      kspsc=co3*ca2plus
c
c
c     - Gleichung (10.3.5), p. 17, UNESCO 87, nach P (=dp) aufloesen:
c     ---------------------------------------------------------------
      t=2.7315d2+tc
c     r=8.31441 J/(mol*K) = 82.0568 atm (cm)**3/(mol*K)
      rt=8.20568d1*t
      dvc=0.5304d0*tc-4.876d1
      dkc=3.692d-4*tc-1.176d-2
      a=0.5d0*dkc/rt
      b=-dvc/rt
      c=dlog(ksp0c/kspsc)
      diskr=b*b-4.d0*a*c
c     dp in atm (oder bar?)
      dp=(dsqrt(diskr)-b)/(2.d0*a)
c
c
c     - Aus dem hydrostatischen Druck die Tiefe der Lysocline berechnen:
c     ------------------------------------------------------------------
c     Einheitenumwandlung in Pa:
      dp=1.01325d5*dp
      getsatho=dp/(density*g)
c
      end
c end function "getsatho"
c
c
      double precision function putponk1_mi(tc,s,dp)
c-----------------------------------------------------------------------c
c     -                                                                -c
c     pressure dependence of carbonic acid diss const. k1              -c
c     depending of temperature, salinity and pressure                  -c
c-----------------------------------------------------------------------c
c     -                                                                -c
c     - putponk1: ratio "K1(dp)/K1(0)"                                 -c
c     --------------------------------                                 -c
c     -                                                                -c
c     - Input:                                                         -c
c     -       tc: Temp. in grad C,                                     -c
c     -       s: in o/oo,                                              -c
c     -       dp in dbar, pressure difference to surface pressure      -c
c     -            of 1atm                                             -c
c     -                                                                -c
c     - Output:                                                        -c
c     -       ratio of K1(dp)/K1(0)
c     -                                                                -c
c     - Millero, 198x
c 
c
c     the code is taken from Chris Sabine routine [pmolvol.r
c      April, 1999
c-----------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, dp
c
c     
c     local variables:
c     ----------------
      double precision pbar,saltrat
      double precision dvf,dkf,dv,dk, rt
c     
c     constants:
c     ----------
c
      double precision dvfresh(3),dvsalt(3),dkfresh(3),dksalt(3)
      data dvfresh /-30.54,    0.1849,     -2.3366d-3 /
      data dvsalt  /-26.69,    0.1976,     -3.1110d-3 /
      data dkfresh / -6.22d-3, 0.1368d-3,   1.2330d-6  / 
      data dksalt  / -3.90d-3, 0.1700d-3,   0.0000d-6  /

c     Initialization:
c     ---------------

      pbar=dp/10.0d0            !pbar = pressure in bar
      rt=(273.15d0+tc)*83.145d0 !RT in cm^3 bar mole^-1
      saltrat=dsqrt(s)/dsqrt(35.0d0)

      dvf = dvfresh(1) + dvfresh(2)*tc + dvfresh(3)*tc*tc
      dkf = dkfresh(1) + dkfresh(2)*tc + dkfresh(3)*tc*tc
      dv = dvf + 
     $     (dvsalt(1) + dvsalt(2)*tc + dvsalt(3)*tc*tc - dvf) 
     $     * saltrat
      dk = dkf + 
     $     (dksalt(1) + dksalt(2)*tc + dksalt(3)*tc*tc - dkf) 
     $     * saltrat
      putponk1_mi = dexp(-dv*pbar/rt + 0.5d0*dk*pbar**2.0d0/rt)
      end
c end of function "putponk1_mi"
c
c
      double precision function putponk2_mi(tc,s,dp)
c-----------------------------------------------------------------------c
c     -                                                                -c
c     pressure dependence of bicarbonate diss const. k2                -c
c     depending of temperature, salinity and pressure                  -c
c-----------------------------------------------------------------------c
c     -                                                                -c
c     - putponk1: ratio "K2(dp)/K2(0)"                                 -c
c     --------------------------------                                 -c
c     -                                                                -c
c     - Input:                                                         -c
c     -       tc: Temp. in grad C,                                     -c
c     -       s: in o/oo,                                              -c
c     -       dp in dbar, pressure difference to surface pressure      -c
c     -            of 1atm                                             -c
c     -                                                                -c
c     - Output:                                                        -c
c     -       ratio of K1(dp)/K1(0)
c     -                                                                -c
c     - Millero, 198x
c 
c
c     the code is taken from Chris Sabine routine [pmolvol.r
c      April, 1999
c-----------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, dp
c
c     
c     local variables:
c     ----------------
      double precision pbar,saltrat
      double precision dvf,dkf,dv,dk, rt
c     
c     constants:
c     ----------
c
      double precision dvfresh(3),dvsalt(3),dkfresh(3),dksalt(3)
      data dvfresh /-29.81,    0.1150,     -1.816d-3 /
      data dvsalt  /-17.32,    0.1758,     -2.647d-3 /
      data dkfresh / -5.74d-3, 0.0930d-3,  -1.896d-6 /
      data dksalt  / -4.87d-3, 0.0900d-3,   0.000d-6 /

c     Initialization:
c     ---------------

      pbar=dp/10.0d0            !pbar = pressure in bar
      rt=(273.15d0+tc)*83.145d0 !RT in cm^3 bar mole^-1
      saltrat=dsqrt(s)/dsqrt(35.0d0)

      dvf = dvfresh(1) + dvfresh(2)*tc + dvfresh(3)*tc*tc
      dkf = dkfresh(1) + dkfresh(2)*tc + dkfresh(3)*tc*tc
      dv = dvf + 
     $     (dvsalt(1) + dvsalt(2)*tc + dvsalt(3)*tc*tc - dvf) 
     $     * saltrat
      dk = dkf + 
     $     (dksalt(1) + dksalt(2)*tc + dksalt(3)*tc*tc - dkf) 
     $     * saltrat
      putponk2_mi = dexp(-dv*pbar/rt + 0.5d0*dk*pbar**2.0d0/rt)
      end
c end of function "putponk2_mi"
c
c
      double precision function putponkw_mi(tc,s,dp)
c-----------------------------------------------------------------------c
c     -                                                                -c
c     pressure dependence of water self diss const.                    -c
c     depending of temperature, salinity and pressure                  -c
c-----------------------------------------------------------------------c
c     -                                                                -c
c     - putponk1: ratio "Kw(dp)/Kw(0)"                                 -c
c     --------------------------------                                 -c
c     -                                                                -c
c     - Input:                                                         -c
c     -       tc: Temp. in grad C,                                     -c
c     -       s: in o/oo,                                              -c
c     -       dp in dbar, pressure difference to surface pressure      -c
c     -            of 1atm                                             -c
c     -                                                                -c
c     - Output:                                                        -c
c     -       ratio of K1(dp)/K1(0)
c     -                                                                -c
c     - Millero, 198x
c 
c
c     the code is taken from Chris Sabine routine [pmolvol.r
c      April, 1999
c-----------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, dp
c
c     
c     local variables:
c     ----------------
      double precision pbar,saltrat
      double precision dvf,dkf,dv,dk, rt
c     
c     constants:
c     ----------
c
      double precision dvfresh(3),dvsalt(3),dkfresh(3),dksalt(3)
      data dvfresh /-25.60,    0.2324,    -3.6246d-3 /
      data dvsalt  /-20.02,    0.1119,    -1.4090d-3 /
      data dkfresh / -7.33d-3, 0.1368d-3, -1.2330d-6 /
      data dksalt  / -5.13d-3, 0.0794d-3,  0.0000d-6 /

c     Initialization:
c     ---------------

      pbar=dp/10.0d0            !pbar = pressure in bar
      rt=(273.15d0+tc)*83.145d0 !RT in cm^3 bar mole^-1
      saltrat=dsqrt(s)/dsqrt(35.0d0)

      dvf = dvfresh(1) + dvfresh(2)*tc + dvfresh(3)*tc*tc
      dkf = dkfresh(1) + dkfresh(2)*tc + dkfresh(3)*tc*tc
      dv = dvf + 
     $     (dvsalt(1) + dvsalt(2)*tc + dvsalt(3)*tc*tc - dvf) 
     $     * saltrat
      dk = dkf + 
     $     (dksalt(1) + dksalt(2)*tc + dksalt(3)*tc*tc - dkf) 
     $     * saltrat
      putponkw_mi = dexp(-dv*pbar/rt + 0.5d0*dk*pbar**2.0d0/rt)
      end
c end of function "putponkw_mi"
c
c
c
      double precision function putponkb_mi(tc,s,dp)
c-----------------------------------------------------------------------c
c     -                                                                -c
c     pressure dependence of borate diss const.                        -c
c     depending of temperature, salinity and pressure                  -c
c-----------------------------------------------------------------------c
c     -                                                                -c
c     - putponk1: ratio "Kb(dp)/Kb(0)"                                 -c
c     --------------------------------                                 -c
c     -                                                                -c
c     - Input:                                                         -c
c     -       tc: Temp. in grad C,                                     -c
c     -       s: in o/oo,                                              -c
c     -       dp in dbar, pressure difference to surface pressure      -c
c     -            of 1atm                                             -c
c     -                                                                -c
c     - Output:                                                        -c
c     -       ratio of K1(dp)/K1(0)
c     -                                                                -c
c     - Millero, 198x
c 
c
c     the code is taken from Chris Sabine routine [pmolvol.r
c      April, 1999
c-----------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision tc, s, dp
c
c     
c     local variables:
c     ----------------
      double precision pbar,saltrat
      double precision dvf,dkf,dv,dk, rt
c     
c     constants:
c     ----------
c
      double precision dvfresh(3),dvsalt(3),dkfresh(3),dksalt(3)
      data dvfresh /-38.75,    0.1991,    -2.813d-3 /
      data dvsalt  /-28.56,    0.1211,    -0.321d-3 /
      data dkfresh /-10.13d-3, 0.2396d-3, -3.961d-6 /
      data dksalt  / -3.00d-3, 0.0427d-3,  0.000d-6 /

c     Initialization:
c     ---------------

      pbar=dp/10.0d0            !pbar = pressure in bar
      rt=(273.15d0+tc)*83.145d0 !RT in cm^3 bar mole^-1
      saltrat=dsqrt(s)/dsqrt(35.0d0)

      dvf = dvfresh(1) + dvfresh(2)*tc + dvfresh(3)*tc*tc
      dkf = dkfresh(1) + dkfresh(2)*tc + dkfresh(3)*tc*tc
      dv = dvf + 
     $     (dvsalt(1) + dvsalt(2)*tc + dvsalt(3)*tc*tc - dvf) 
     $     * saltrat
      dk = dkf + 
     $     (dksalt(1) + dksalt(2)*tc + dksalt(3)*tc*tc - dkf) 
     $     * saltrat
      putponkb_mi = dexp(-dv*pbar/rt + 0.5d0*dk*pbar**2.0d0/rt)
      end
c end of function "putponkb_mi"
c
c
c
c
c
c
c
      double precision function putponk1(tc,s,dp)
c----------------------------------------------------------------------------c
c     -                                                                     -c
c     - putponk1: ratio "K1(dp)/K1(0)"                                      -c
c     --------------------------------                                      -c
c     -                                                                     -c
c     - Referenz:                                                           -c
c     - Thermodynamics of the carbon dioxide system in seawater,            -c
c     - UNESCO technical papers in marine science, 51, 1987, p. 14          -c
c     -                                                                     -c
c     - dp: Druckdifferenz zu 1 atm (evtl. 1 bar)                           -c
c     -     1 atm = 1.01325 10**5 Pa, 1 bar = 10**5 Pa                      -c
c     -                                                                     -c
c     - Zu den Einheiten:                                                   -c
c     - [dp] = Pa,  [dv]=(cm)**3/mol; [dk]=(cm)**3/(mol*atm)                -c
c     - Avogadrokonstante: [r]=atm*(cm)**3/(mol*K), somit:                  -c
c     - r=8.31441 J/(mol*K) = 82.0568 atm (cm)**3/(mol*K)                   -c
c     -                                                                     -c
c     - Bem: Aus Literatur nicht eindeutig klar, ob atm oder bar als Einheit-c
c     -      des Druckes fuer die Groessen dv und dk.                       -c
c     -      Widerspruechliche Angaben!                                     -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision dp, tc, s
c
c     constants:
c     ----------
      double precision a0, a1, a2, b0, b1, b2
      parameter(a0=2.55d1, a1=0.151d0, a2=-0.1271d0)
      parameter(b0=3.08d0, b1=0.578d0, b2=-8.77d-2)
c
c     
c     local variables:
c     ----------------
      double precision rt, dv, dk, dpinatm
c
c
c     Umwandlung von Pa in atm:
      dpinatm=dp/1.01325d5
c
      rt=8.20568d1*(tc+2.7315d2)
      dv=-a0-a1*(s-34.8d0)-a2*tc
      dk=(b0+b1*(s-34.8d0)+b2*tc)*(-10.d-3)
      putponk1=dexp( (-dv/rt)*dpinatm+(0.5d0*dk/rt)*dpinatm*dpinatm )
c
      end
c end of function putponk1
c
c
c
c
c
      double precision function putponk2(tc,s,dp)
c----------------------------------------------------------------------------c
c     -                                                                     -c
c     - putponk2: ratio "K2(dp)/K2(0)"                                      -c
c     --------------------------------                                      -c
c     -                                                                     -c
c     - Referenz:                                                           -c
c     - Thermodynamics of the carbon dioxide system in seawater,            -c
c     - UNESCO technical papers in marine science, 51, 1987, p. 14          -c
c     -                                                                     -c
c     - dp: Druckdifferenz zu 1 atm (evtl. 1 bar)                           -c
c     -     1 atm = 1.01325 10**5 Pa, 1 bar = 10**5 Pa                      -c
c     -                                                                     -c
c     - Zu den Einheiten:                                                   -c
c     - [dp] = Pa,  [dv]=(cm)**3/mol; [dk]=(cm)**3/(mol*atm)                -c
c     - Avogadrokonstante: [r]=atm*(cm)**3/(mol*K), somit:                  -c
c     - r=8.31441 J/(mol*K) = 82.0568 atm (cm)**3/(mol*K)                   -c
c     -                                                                     -c
c     - Bem: Aus Literatur nicht eindeutig klar, ob atm oder bar als Einheit-c
c     -      des Druckes fuer die Groessen dv und dk.                       -c
c     -      Widerspruechliche Angaben!                                     -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision dp, tc, s
c
c     constants:
c     ----------
      double precision a0, a1, a2, b0, b1, b2
      parameter(a0=1.582d1, a1=-0.321d0, a2=2.19d-2)
      parameter(b0=-1.13d0, b1=0.314d0, b2=0.1475d0)
c
c     
c     local variables:
c     ----------------
      double precision rt, dv, dk, dpinatm
c
c
c     Umwandlung von Pa in atm:
      dpinatm=dp/1.01325d5
c
      rt=8.20568d1*(tc+2.7315d2)
      dv=-a0-a1*(s-34.8d0)-a2*tc
      dk=(b0+b1*(s-34.8d0)+b2*tc)*(-10.d-3)
      putponk2=dexp( (-dv/rt)*dpinatm+(0.5d0*dk/rt)*dpinatm*dpinatm )
c
      end
c end of function putponk2
c
c
c
c
      double precision function putponkb(tc,s,dp)
c----------------------------------------------------------------------------c
c     -                                                                     -c
c     - putponk2: ratio "Kb(dp)/Kb(0)"                                      -c
c     --------------------------------                                      -c
c     -                                                                     -c
c     - Referenz:                                                           -c
c     - Thermodynamics of the carbon dioxide system in seawater,            -c
c     - UNESCO technical papers in marine science, 51, 1987, p. 14          -c
c     -                                                                     -c
c     - dp: Druckdifferenz zu 1 atm (evtl. 1 bar)                           -c
c     -     1 atm = 1.01325 10**5 Pa, 1 bar = 10**5 Pa                      -c
c     -                                                                     -c
c     - Zu den Einheiten:                                                   -c
c     - [dp] = Pa,  [dv]=(cm)**3/mol; [dk]=(cm)**3/(mol*atm)                -c
c     - Avogadrokonstante: [r]=atm*(cm)**3/(mol*K), somit:                  -c
c     - r=8.31441 J/(mol*K) = 82.0568 atm (cm)**3/(mol*K)                   -c
c     -                                                                     -c
c     - Bem: Aus Literatur nicht eindeutig klar, ob atm oder bar als Einheit-c
c     -      des Druckes fuer die Groessen dv und dk.                       -c
c     -      Widerspruechliche Angaben!                                     -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision dp, tc, s
c
c     constants:
c     ----------
      double precision a0, a1, a2, a3, b0, b1
      parameter(a0=2.948d1, a1=-0.295d0, a2=-0.1622d0, a3=2.608d-3)
      parameter(b0=2.84d0, b1=-0.354d0)
c
c     
c     local variables:
c     ----------------
      double precision rt, dv, dk, dpinatm
c
c
c     Umwandlung von Pa in atm:
      dpinatm=dp/1.01325d5
c
      rt=8.20568d1*(tc+2.7315d2)
      dv=-a0-a1*(s-34.8d0)-a2*tc-a3*tc*tc
      dk=(b0+b1*(s-34.8d0))*(-10.d-3)
      putponkb=dexp( (-dv/rt)*dpinatm+(0.5d0*dk/rt)*dpinatm*dpinatm )
c
      end
c end of function putponkb
c
c
c
      double precision function putponkbOLD(tc,dp)
c----------------------------------------------------------------------------c
c     -                                                                     -c
c     - putponkb: ratio "Kb(dp)/Kb(0)"                                      -c
c     --------------------------------                                      -c
c     -                                                                     -c
c     - Referenz:                                                           -c
c     - Culberson C., D.R. Kester and R.M. Pytkowicz (1967),                -c
c     - High pressure dissociation of carbonic and boric acids in seawater, -c
c     - science 157, p. 59 - 61                                             -c
c     -                                                                     -c
c     - dp: Druckdifferenz zu 1 atm (evtl. 1 bar)                           -c
c     -     1 atm = 1.01325 10**5 Pa, 1 bar = 10**5 Pa                      -c
c     -                                                                     -c
c     - Zu den Einheiten:                                                   -c
c     - [dp] = Pa                                                           -c
c     -                                                                     -c
c     - Bem: Aus Literatur nicht eindeutig klar, ob atm oder bar als Einheit-c
c     -      des Druckes fuer die Groessen dv und dk.                       -c
c     -      Widerspruechliche Angaben!                                     -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision dp, tc
c
c     constants:
c     ----------
      double precision c1, c2
      parameter(c1=27.5d0, c2=9.5d-2)
c
c     
c     local variables:
c     ----------------
      double precision dpinatm, cp
c
c
c     Umwandlung von Pa in atm:
      dpinatm=dp/1.01325d5
c
      cp = dpinatm/(83.143d0*(tc+273.15d0))
c
      putponkbOLD=dexp((c1-c2*tc)*cp)
c
      end
c end of function putponkbOLD
c
c
c
c
      double precision function putponkw(tc,s,dp)
c----------------------------------------------------------------------------c
c     -                                                                     -c
c     - putponkw: ratio "Kw(dp)/Kw(0)"                                      -c
c     --------------------------------                                      -c
c     -                                                                     -c
c----------------------------------------------------------------------------c
c
      implicit none
c
c     global variable:
c     ----------------
      double precision dp, tc, s
c
      putponkw = 1.d0
c
      end
c end of function putponkw
