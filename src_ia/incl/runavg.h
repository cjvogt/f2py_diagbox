c======================= include file runavg.h ==============================
c
c variables that store running averages
c
c REVISIONS:
c
c date     author    remarks
c
c 08.06.00 ng        first implementation
c
c variables:
c
c name              type   description
c-----------------------------------------------------------------------------
c 
c D_C_runavg        dp a   one year running avg Delta C due to the 
c                                 different proc.                       muM/kg/ts
c fluxes_runavg     dp a   one year running avg Flux of DIC due 
c-----------------------------------------------------------------------------
c
      double precision D_C_runavg(nproc_rates,nstepmax)
      double precision fluxes_runavg(nproc_fluxes,nstepmax)
