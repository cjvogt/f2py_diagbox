c======================== include file input.h  =============================
c
c input data for the seasonal model, stemming from the harmonic 
c      regression (except C, dC13 and h which are defined in slabs.h)
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c
c variables:
c
c name       type   description                                   units
c-----------------------------------------------------------------------------
c
c temp        dp s  mixed layer temperature                       deg C
c sal         dp s  mixed layer salinity                          psu
c ws          dp s  wind-speed                                    m s-1
c fCO2_o      dp s  CO2 fugacity in the mixed layer               uatm
c pco2_a      dp s  pCO2 in the atmosphere                        ppm
c dc13_a      dp s  dC13 in the atmosphere                        per mil
c dCdz_var    dp s  vertical C gradient below the mixed layer   umol/kg/m
c kz_var      dp s  vertical diffusion coefficient                m2 s-1
c salk        dp s  salinity normalized alkalinity              umol/kg
c 
c-----------------------------------------------------------------------------
c
      double precision temp,sal,ws,fCO2_o,pco2_a,dc13_a
      double precision dCdz_var,kz_var,salk


