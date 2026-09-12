c================== include file mc_statistics.h =========================
c
c variables of statistics of monte carlo runs
c
c REVISIONS:
c
c date     author    remarks
c
c 19.8.96  ng        first implementation
c 18.10.99 ng        adapted for interannual variability runs
c
c variables:
c
c name            type   description
c---------------------------------------------------------------------------
c  
c int_rates_std   dp a   array holding standard integrated rates
c int_fluxes_std  dp a   array holding standard integrated fluxes
c costfn_std      dp s   value of cost function for standard run
c rsq_std         dp s   r**2 value for standard run
c chisq_std       dp s   chi**2 value for standard run
c amplitude_std   dp s   amplitude of standard run
c c_min_std       dp s   minimum C of standard run
c c_max_std       dp s   maximum C of standard run
c
c var_int_rates   dp a   array holding variance of integrated rates 
c                            (over monte carlo runs)
c var_int_fluxes  dp a   array holding variance of integrated fluxes
c wvar_int_rates  dp a   array holding weighted variance of integrated rates 
c                                 (over monte carlo runs)
c wvar_int_fluxes dp a   array holding weighted variance of integrated fluxes
c sum_costfn      dp s   sum of cost function
c mean_costfn     dp s   mean of cost function
c costfn_min      dp s   minimum of cost function
c costfn_max      dp s   maximum of cost function
c 
c sg_int_rates    dp s   1-sigma uncertainty of int_rates
c sg_int_fluxes   dp s   1-sigma uncertainty of int_fluxes
c wsg_int_rates   dp s   1-sigma weighted uncertainty of int_rates
c wsg_int_fluxes  dp s   1-sigma weighted uncertainty of int_fluxes
c 
c===========================================================================
c
      double precision int_rates_std(nyearmax,nproc_rates)
      double precision int_fluxes_std(nyearmax,nproc_fluxes)
      double precision costfn_std
      double precision closure_std,rsq_std,chisq_std,amplitude_std
      double precision c_min_std,c_max_std
c
      double precision var_int_rates(nyearmax,nproc_rates)
      double precision var_int_fluxes(nyearmax,nproc_fluxes)
      double precision wvar_int_rates(nyearmax,nproc_rates)
      double precision wvar_int_fluxes(nyearmax,nproc_fluxes)
c   
      double precision sg_int_rates(nyearmax,nproc_rates)
      double precision sg_int_fluxes(nyearmax,nproc_fluxes)
      double precision wsg_int_rates(nyearmax,nproc_rates)
      double precision wsg_int_fluxes(nyearmax,nproc_fluxes)
c
      double precision sum_costfn,mean_costfn
      double precision costfn_min,costfn_max
c

