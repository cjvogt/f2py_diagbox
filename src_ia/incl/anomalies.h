c======================= include file averages.h ==============================
c
c variables that store anomalies
c
c REVISIONS:
c
c date     author    remarks
c
c 08.06.00 ng        first implementation
c
c variables:
c
c name            type   description
c-----------------------------------------------------------------------------
c 
c D_C_avg         dp a   average Delta C due to the different proc.       muM/kg/ts
c D_C_anom        dp a   anomalies of Delta C                             muM/kg/ts   
c fluxes_avg      dp a   average Flux of DIC due to the different proc.   mol C/m2/s
c fluxes_anom     dp a   anomalies of Delta C                             mol C/m2/s
c D_C_a_runavg    dp a   365 days running anomalies of Delta C            muM/kg/ts   
c fluxes_a_runavg dp a   365 days running anomalies of fluxes             mol C/m2/s 
c
c-------------------------------------------------------------------------
c
      double precision D_C_avg(nproc_rates,nstepmax)
      double precision D_C_anom(nproc_rates,nstepmax)
      double precision fluxes_avg(nproc_fluxes,nstepmax)
      double precision fluxes_anom(nproc_fluxes,nstepmax)
      double precision D_C_a_runavg(nproc_rates,nstepmax)
      double precision fluxes_a_runavg(nproc_fluxes,nstepmax)
