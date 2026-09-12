c
c======================= include file intval.h ==============================
c
c variables for calculating the time integrated fluxes and changes
c in concentrations. 
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 17.10.99 ng        changed nperiods to nyears
c
c variables:
c
c name          type   description
c-----------------------------------------------------------------------------
c 
c Int_Fluxes    dp a   two dimensional array holding the time integrated
c                       fluxes over the mixed layer.
c                       unit: molC/m2 
c
c Int_Rates     dp a   two dimensional array holding the time integrated
c                       rates of change
c                       unit: umolC/kg 
c-------------------------------------------------------------------------
c
      double precision Int_Fluxes(nyearmax,nproc_fluxes)
      double precision Int_Rates(nyearmax,nproc_rates)



