c================== include file sensitivity.h =========================
c
c variables for sensitivity analysis
c
c REVISIONS:
c
c date     author    remarks
c
c 21.8.96  ng        first implementation
c 18.10.99 ng        adapted for interannual variability runs
c
c variables:
c
c name            type   description
c---------------------------------------------------------------------------
c  
c d_int_rates_dp  dp a   array holding derivative of integrated rates with
c                          respect to parameter p
c d_int_fluxes_dp dp a   array holding derivative of integrated fluxes with
c                          respect to parameter p
c 
c===========================================================================
c
      double precision d_int_rates_dp(nparams,nyearmax,nproc_rates)
      double precision d_int_fluxes_dp(nparams,nyearmax,nproc_fluxes)
