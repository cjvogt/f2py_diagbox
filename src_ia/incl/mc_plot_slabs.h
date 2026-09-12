c================== include file mc_plot_slabs.h =========================
c
c variables that contain plot information for monte carlo runs
c
c REVISIONS:
c
c date     author    remarks
c
c 19.8.96  ng        first implementation
c
c variables:
c
c name                type   description
c---------------------------------------------------------------------------
c  		      
c nstepmax2           i  s   maximum number of values on time axis
c nstep               i  s   number of resampled values
c rsmpl_step          i  s   every i'th sample is resampled
c mc_D_c              dp a   array of resampled D_C array for mc analysis
c mc_fluxes           dp a   array of resampled fluxes array for mc analysis
c mc_c_sim            dp a   array of resampled c_sim array for mc analysis
c mc_u_var            dp a   array of resampled inverse calc u
c mc_ddc13_dh_var     dp a   array of resampled inverse calc ddC13_dh
c mean_D_c            dp a   mean D_C of mc runs
c mean_fluxes         dp a   mean fluxes of mc runs
c mean_c_sim          dp a   mean c_sim of mc runs
c mean_u_var          dp a   mean u_var
c mean_ddc13_dh_var   dp a   mean ddc13_dh_var
c wmean_D_c           dp a   weighted mean D_C of mc runs
c wmean_fluxes        dp a   weighted mean fluxes of mc runs
c wmean_c_sim         dp a   weighted mean c_sim of mc runs
c wmean_u_var         dp a   weighted mean u_var
c wmean_ddc13_dh_var  dp a   weighted mean ddc13_dh_var
c sg_D_c              dp a   1-sg uncertainty of D_C
c sg_fluxes           dp a   1-sg uncertainty of fluxes
c sg_c_sim            dp a   1-sg uncertainty of c_sim
c sg_u_var            dp a   1-sg uncertainty of u_var
c sg_ddc13_dh_var     dp a   1-sg uncertainty of ddc13_dh_var
c wsg_D_c             dp a   weighted 1-sg uncertainty of D_C
c wsg_fluxes          dp a   weighted 1-sg uncertainty of fluxes
c wsg_c_sim           dp a   weighted 1-sg uncertainty of c_sim
c wsg_u_var           dp a   weighted 1-sg uncertainty of u_var
c wsg_ddc13_dh_var    dp a   weighted 1-sg uncertainty of ddc13_dh_var
c c_sim_std           dp a   array holding c_sim of standard simulation
c fluxes_std          dp a   array holding fluxes of standard simulation
c d_c_std             dp a   array holding d_c of standard simulation
c u_var_std           dp a   array holding inv. calc. u of std simulation
c ddc13_dh_std        dp a   array holding inv. calc. ddc13_dh of std sim.
c
c                                 to the different proc.                mol C/m2/s
c D_C_std_runavg    dp a   one year running standard Delta C due to the 
c                                 different proc.                       muM/kg/ts
c fluxes_std_runavg dp a   one year running standard Flux of DIC due 
c                                 to the different proc.                mol C/m2/s
cc 
c===========================================================================
c
      integer nstepmax2,nstep2,rsmpl_step
      parameter(nstepmax2 = 200)
c
      double precision mc_D_c(nrunmax,nproc_rates,nstepmax2)
      double precision mc_fluxes(nrunmax,nproc_fluxes,nstepmax2)
      double precision mc_c_sim(nrunmax,nstepmax2)
      double precision mean_D_c(nproc_rates,nstepmax2)
      double precision mean_fluxes(nproc_fluxes,nstepmax2)
      double precision mean_c_sim(nstepmax2)
      double precision wmean_D_c(nproc_rates,nstepmax2)
      double precision wmean_fluxes(nproc_fluxes,nstepmax2)
      double precision wmean_c_sim(nstepmax2)
      double precision sg_D_c(nproc_rates,nstepmax2)
      double precision sg_fluxes(nproc_fluxes,nstepmax2)
      double precision sg_c_sim(nstepmax2)
      double precision wsg_D_c(nproc_rates,nstepmax2)
      double precision wsg_fluxes(nproc_fluxes,nstepmax2)
      double precision wsg_c_sim(nstepmax2)
      double precision mc_weight(nrunmax)
      double precision c_sim_std(nstepmax)
      double precision d_c_std(nproc_rates,nstepmax)
      double precision fluxes_std(nproc_fluxes,nstepmax)
c
      double precision mc_u_var(nrunmax,nstepmax2)
      double precision mc_ddc13_dh_var(nrunmax,nstepmax2)
      double precision mean_u_var(nstepmax2)
      double precision mean_ddc13_dh_var(nstepmax2)
      double precision wmean_u_var(nstepmax2)
      double precision wmean_ddc13_dh_var(nstepmax2)
      double precision sg_u_var(nstepmax2)
      double precision sg_ddc13_dh_var(nstepmax2)
      double precision wsg_u_var(nstepmax2)
      double precision wsg_ddc13_dh_var(nstepmax2)
      double precision u_var_std(nstepmax)
      double precision ddc13_dh_var_std(nstepmax)
c
      double precision fluxes_std_runavg(nproc_fluxes,nstepmax), 
     $     D_C_std_runavg(nproc_rates,nstepmax)
      double precision mean_d_c_runavg(nproc_rates,nstepmax2),
     $        sg_D_c_runavg(nproc_rates,nstepmax2),       
     $        mean_fluxes_runavg(nproc_fluxes,nstepmax2),
     $        sg_fluxes_runavg(nproc_fluxes,nstepmax2)   

