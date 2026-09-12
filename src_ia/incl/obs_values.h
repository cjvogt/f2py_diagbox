c==================== include file obs_values.h  =============================
c
c observed seasonal cycles of input data 
c    (except C, dC13 and h which are defined in slabs.h)
c
c REVISIONS:
c
c date     author    remarks
c
c 08.08.96  ng       first implementation
c 07.06.00  ng       added salk_obs
c
c variables:
c
c name       type   description                                   units
c-----------------------------------------------------------------------------
c
c temp_obs     dp a  mixed layer temperature                       deg C
c sal_obs      dp a  mixed layer salinity                          psu
c ws_obs       dp a  wind-speed                                    m s-1
c fco2_o_obs   dp a  CO2 fugacity in the mixed layer               uatm
c pco2_a_obs   dp a  xco2 in the atmosphere
c fco2_a_obs   dp a  fCO2 in the atmosphere                        uatm
c pco2_a_dry   dp a  pCO2 in the atmosphere                        
cc dc13_a_obs   dp a  dC13 in the atmosphere                        per mil
c dCdz_obs     dp a  vertical C gradient below the mixed layer     umol/kg/m
c kz_obs       dp a  vertical diffusivity                          m2 s-1
c salk_obs     dp a  salinity normalized alk                       umol/kg
c 
c-----------------------------------------------------------------------------
c
      double precision temp_obs(nstepmax),sal_obs(nstepmax)
      double precision ws_obs(nstepmax),fco2_o_obs(nstepmax)
      double precision pco2_a_obs(nstepmax),fco2_a_obs(nstepmax)
      double precision pco2_a_dry(nstepmax)
      double precision dc13_a_obs(nstepmax)
      double precision dCdz_obs(nstepmax),kz_obs(nstepmax)
      double precision salk_obs(nstepmax)


