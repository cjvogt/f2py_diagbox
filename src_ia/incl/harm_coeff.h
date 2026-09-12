c===================== include file harm_coeff.h =========================
c
c coefficients of the harmonic fits to the input data
c
c REVISIONS:
c
c date     author    remarks
c
c 30.7.96  ng        first fits (not yet checked)
c                      for wind-speed : still using old fits
c 10.9.96  ng        updated coefficients, using only coeff which
c                      are significant at the 5% level (tested by
c                      Student t-test)
c                    changes include temperature,ddic_dz and dc13_atm
c 17.10.99 ng        adapted for interannual variability runs at BERM
c 06.06.00 ng        added new coefficients based on combined BERM/BATS
c                      data set
c 24.08.01 ng        added new coefficients for 1983-2001
c 03.05.02 ng        added new coefficients for 1983 to 2002
c
c variables:
c
c name        type   description
c-------------------------------------------------------------------------
c
c c_temp      dp  a   coefficients for mixed layer temperature    [deg C]
c c_sal       dp  a   coefficients for mixed layer salinity       [psu]
c c_sdic      dp  a   coefficients for mixed layer normalized C   [umol/kg]
c c_salk      dp  a   coefficients for mixed layer normalized Alk [umol/kg]
c c_dc13      dp  a   coefficients for mixed layer dC13           [permil]
c c_fco2_o    dp  a   coefficients for mixed layer fCO2           [ppm]
c c_pco2a     dp  a   coefficients for atm pCO2                   [ppm]
c c_dc13a     dp  a   coefficients for atm dC13                   [per mil]
c c_ws        dp  a   coefficients for wind wind speed            [m s-1]
c                          fit                                    [m2 s-2]
c c_h         dp  a   coefficients for mld                        [m]
c c_dCdz      dp  a   coefficients for dC/dz                      [umol/kg m]
c c_kz_var    dp  a   coefficients for seasonally varying kz      [m2 s-1]
c
c===========================================================================
c
      double precision c_temp(7),c_sal(7)
      double precision c_sdic(7),c_salk(7)
      double precision c_dc13(7),c_fco2_o(7)
      double precision c_pco2a(7),c_dc13a(7)
      double precision c_h(7)
      double precision c_dCdz(7)
      double precision c_kz_var(7),c_ws(7)
c
c     
c --------------------------------------------------------------------
c     BATS specific coefficients (1989-1993)
c --------------------------------------------------------------------
c
#ifdef BATS_89_93

      data c_temp    / 0.23506d+2,  
     $                -0.3197d+1,-0.2429d+1, 0.5994d+0, 0.4677d+0,
     $                 0.0d+0,0.0d+0/
c
      data c_sal     / 0.36586d+2,
     $                 0.8454d-1,0.2748d-1,0.0d0,0.0d0,
     $                 0.0d+0,0.0d+0/
c
      data c_sdic   / 0.20421d4,
     $                 0.1556d+2, 0.9945d+0, 0.5501d+0, 0.2533d+0,
     $                 0.0d+0,0.0d+0/
c
      data c_dc13   / 0.13496d+1,
     $                -0.1025d+0,0.59113d-2,-0.19802d-1,-0.11080d-1,
     $                 0.0d+0,0.0d+0/
c
      data c_salk   / 0.2392d+4,
     $                 0.0d+0,0.0d+0,0.0d+0,0.0d+0,
     $                 0.0d+0,0.0d+0/
c
      data c_fco2_o  / 0.34644d+3,
     $                -0.2056d+2,-0.2703d+2, 0.8252d+1, 0.6674d+1,
     $                 0.0d+0,0.0d+0/
c
      data c_pco2a  / 0.35689d+3,
     $                0.4351d+1, 0.1126d+1,-0.2002d+1, 0.7626d+0,
     $                0.3959d+0,-0.3621d+0/
c
      data c_dc13a  /-0.78753d+1,
     $                -0.2384d+0, -0.7775d-1, 0.1085d+0,-0.3815d-1,
     $                 0.0d0,0.0d0/
c
      data c_kz_var  / 0.3624d-4 , 
     $                 0.4306d-4 , 0.8956d-5 ,0.5444d-5, -0.1909d-4,
     $                -0.7567d-5,  0.8822d-6/
c
       data c_dCdz    /0.50410d+0,
     $               -0.1578d+0, 0.1728d+0, 0.0d+0,0.0d+0,
     $                0.0d+0,0.0d+0/
#elif BERM_BATS_83_98
c     
c --------------------------------------------------------------------
c     Station 'S' and BATS combined coefficients (1983-1998)
c --------------------------------------------------------------------

      data c_temp    /0.23213d+2,
     $               -0.3438d+1,-0.2158d+1, 0.6349d+0, 0.2140d+0,
     $               -0.1516d+0,-0.1382d+0/
c
      data c_sal     /0.36565d+2,
     $                0.8409d-1, 0.5096d-1,-0.1457d-1, 0.8162d-2,
     $               -0.1859d-1, 0.2561d-1/
c
      data c_sdic    /0.19518d+4,  !@salinity 35
     $                0.1469d+2, 0.3186d+0, 0.1453d+1, 0.6547d+0,
     $                0.1597d+0,-0.9035d+0/
c
      data c_dc13    /0.14155d+1,
     $               -0.9249d-1,-0.1181d-2,-0.2559d-1,-0.5344d-2,
     $                0.3496d-2, 0.3655d-2/
c
      data c_salk    /0.22858d+4,  !@salinity 35
     $               -0.3945d-1, 0.6558d+0, 0.4173d+0,-0.1268d+0,
     $               -0.1957d+0, 0.1699d+0/
c
      data c_fco2_o / 0.33319d+3,
     $               -0.2319d+2,-0.3041d+2, 0.9512d+1, 0.5077d+1,
     $               -0.2230d+1,-0.3026d+1/
c
      data c_pco2a  / 0.35689d+3,  ! same as BATS
     $                 0.4351d+1, 0.1126d+1,-0.2002d+1, 0.7626d+0,
     $                 0.3959d+0,-0.3621d+0/
c
      data c_dc13a  /-0.78753d+1, ! same as BATS
     $                -0.2384d+0, -0.7775d-1, 0.1085d+0,-0.3815d-1,
     $                 0.0d0,0.0d0/
c
      data c_kz_var  / 0.3624d-4 , ! same as BATS
     $                 0.4306d-4 , 0.8956d-5 ,0.5444d-5, -0.1909d-4,
     $                -0.7567d-5,  0.8822d-6/
c
      data c_dCdz    /0.50410d+0, ! same as BATS
     $               -0.1578d+0, 0.1728d+0, 0.0d+0,0.0d+0,
     $                0.0d+0,0.0d+0/
#else
c     
c --------------------------------------------------------------------
c     Station 'S' and BATS combined coefficients (1983-2001)
c       updated May 03, 2002
c --------------------------------------------------------------------

      data c_temp    /0.23231d+2,
     $               -0.3491d+1,-0.2171d+1, 0.6304d+0, 0.2167d+0,
     $               -0.2627d-1,-0.8370d-1/
c
      data c_sal     /0.36557d+2,
     $                0.9186d-1, 0.3998d-1,-0.1776d-1, 0.9852d-2,
     $               -0.6274d-2, 0.2737d-1/
c
      data c_sdic    /0.19523d+4,  !@salinity 35
     $                0.1437d+2, 0.1247d+1, 0.9282d+0, 0.6024d+0,
     $               -0.2430d+0,-0.5423d+0/
c
      data c_dc13    /0.13848d+1,
     $               -0.8406d-1,-0.7534d-3,-0.1856d-1,-0.1059d-2,
     $               -0.1306d-2, 0.2271d-2/
c
      data c_salk    /0.22852d+4,  !@salinity 35
     $               -0.3019d+0, 0.7049d+0, -0.2551d-1,-0.3557d+0,
     $               -0.1768d-1, 0.4185d+0/
c
      data c_fco2_o / 0.33422d+3,
     $               -0.2340d+2,-0.2835d+2, 0.1045d+2, 0.4631d+1,
     $               -0.1253d+1,-0.2027d+1/
c
      data c_pco2a  / 0.35645d+3,  ! new coeff May 03, 2002
     $                 0.4073d+1, 0.1252d+1, -0.1315d+1, 0.4781d+0,
     $                -0.2101d-1,-0.3498d+0/
c
      data c_dc13a  /-0.78625d+1, ! new coeff May 03, 2002
     $                -0.2239d+0, -0.9744d-1, 0.6619d-1,-0.1664d-1,
     $                 0.6341d-2,  0.1724d-1/
c
      data c_kz_var  / 0.3624d-4 , ! same as BATS, not updated
     $                 0.4306d-4 , 0.8956d-5 ,0.5444d-5, -0.1909d-4,
     $                -0.7567d-5,  0.8822d-6/

      data c_dCdz    /0.6661d+0, ! BATS 1988 - 2001; stations 1-156
     $               -0.1969d+0, 0.2734d+0, 0.2231d-1, 0.4882d-1,
     $               -0.6619d-2,-0.2185d-1/


#endif
c
      
c     
c --------------------------------------------------------------------
c     unspecific coefficients or user able to choose interactively
c --------------------------------------------------------------------
c
#ifdef BATS_89_93
c Isemer and Hasse coefficients
c
      data c_ws   /   7.173d+0 , 
     $                0.764d+0, 1.894d+0, 0.057d+0, -0.279d+0,
     $                0.0d+0,0.0d+0/
#else
c NCEP coefficients data from 1983-1996
c
c       data c_ws   /  0.63541d+1,
c     $                0.7353d+0, 0.1641d+1, 0.6329d-1, 0.1882d-1,
c     $               -0.9492d-1,-0.1030d+0/   
c NCEP monthly coefficients data from 1959 - 2000
c
       data c_ws   /  0.7531d+1,
     $                0.6781d+0, 0.1967d+1, 0.1809d+0, 0.1424d+0,
     $               -0.1997d+0,-0.1780d+0/   
#endif
c
#ifdef BATS_89_93
c BATS averaged ctd coefficients
c
      data c_h     / 0.74815d+2,
     $               0.4368d+2, 0.5597d+2, 0.2176d+2,-0.9656d+1,
     $              -0.3274d+0,-0.7658d+1/
#else
cc
cc Station 'S' hydrographic data coefficients
cc
c       data c_h    /  0.77770d+2,
c     $                0.5468d+2, 0.6258d+2, 0.2984d+2,-0.1631d+2,
c     $               -0.6115d+1,-0.1778d+2/
c
c BATS hydrographic data coefficients, BATS 1-156, updated May 3, 2002
c
       data c_h    /  0.7217d+2,
     $                0.4302d+2, 0.5158d+2, 0.2149d+2,-0.9637d+1,
     $               -0.1355d+1,-0.9183d+1/
#endif

