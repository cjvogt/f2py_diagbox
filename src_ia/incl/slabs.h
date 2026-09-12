c======================== include file slabs.h  =============================
c
c main variables for the seasonal carbon model:
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c
c variables:
c
c name       type  description                                   units
c-----------------------------------------------------------------------------
c
c C_new      dp s  obs DIC at time step t+ts                     muM/kg
c C_old      dp s  obs DIC at time step t                        muM/kg
c C          dp s  obs DIC at time step t+0.5*t                  muM/kg
c C_obs      dp a  array of observed DIC                         muM/kg
c C_sim      dp a  array of simulated DIC                        muM/kg
c	    
c C13_new    dp s  obs. C13 at time step t+ts                    muM/kg
c C13_old    dp s  obs. C13 at time step t                       muM/kg
c C13        dp s  obs. C13 at time step t+0.5*ts                muM/kg
c dC13_obs   dp a  array of observed dC13                        per mil
c      	     
c h_new      dp s  mixed layer depth at time step t+ts           m
c h_old      dp s  mixed layer depth at time step t              m
c h          dp s  mixed layer depth at time step t+0.5*ts       m
c h_obs      dp a  array of observed mixed layer depths          m
c
c D_C        dp a   Delta C due to the different processes       muM/kg/ts
c                    DDIC(nproc_rates,ts)
c D_C13      dp a   Delta C13 due to the different processes     per mil/ts
c                    DdC13(nproc_rates,ts)
c fluxes     dp a   Flux of DIC due to the different processes   mol C/m2/s
c                    Flux(nproc_fluxes,ts)
c
c-----------------------------------------------------------------------------
c
      double precision C_old,C_new,C13_old,C13_new,h_new,h_old
      double precision C,c13,h
      double precision C_obs(nstepmax),C_sim(nstepmax)
      double precision dC13_obs(nstepmax),h_obs(nstepmax)
c
      double precision D_C(nproc_rates,nstepmax)
      double precision D_C13(nproc_rates,nstepmax)
      double precision fluxes(nproc_fluxes,nstepmax)


