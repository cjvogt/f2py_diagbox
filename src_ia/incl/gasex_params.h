c===================== include file gasex_params.h =========================
c
c parameters of air-sea exchange
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 13.8.96  ng        changed l_* to sg_*
c 06.6.00  ng        added comp_foc2_online
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c 
c pv_relship     c  s   character indicating the relationship between
c                       piston velocity and wind-speed that is used for
c                       calculating the piston-velocity
c                           <lm> : Liss and Merlivat
c                           <wa> : Wanninkhof  
c D_fco2_corr    dp s   correction term for calculated fCO2 ocean     [uatm]
c gasex_fact     dp s   multiplication factor for air-sea exchange
c                           usually 1.0 for Wanninkhof
c                           usually 1.7447 for Liss and Merlivat
c
c ws_coeff       c  s   character specifying wind speed coefficients
c                           <ih> : Isemer and Hasse, 1985
c                           <ae> : AEROCE, 1989-1991
c
c sg_D_fco2_corr dp a   uncertainty of D_fco2_corr
c sg_gasex_fact  dp a   uncertainty of gasex_fact
c
c mean_air_press dp s   mean air pressure [atm] (from St. Davids head)
c mean_rel_humid dp s   mean relative humidity [] (from St. Davids head)
c 
c comp_fco2_online l s  logical indicating whether fco2 should be
c                         computed on or offline
c
c===========================================================================
c
      character pv_relship*2,ws_coeff*2
      double precision D_fco2_corr,gasex_fact
      double precision sg_D_fco2_corr, sg_gasex_fact
c
      double precision mean_air_press,mean_rel_humid
      parameter (mean_air_press = 1.017,mean_rel_humid = 0.87)
c
      logical comp_fco2_online
