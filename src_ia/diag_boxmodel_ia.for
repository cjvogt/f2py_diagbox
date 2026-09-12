# 1 "diag_boxmodel_ia.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/home/prafter/f2py-diagbox/src_ia//"
# 1 "diag_boxmodel_ia.F"
      program diag_boxmodel_ia
c
c =====================================================================
c
c     PURPOSE:
c
c     This program is the fortran version of the diagnostic seasonal 
c     box model originally developed for station 'S' data. 
c     This is a modified version of the seasonal model, intended to run
c     with interannual data. (October 1999). The calculation scheme has
c     been changed too. Instead of solving von dC/dt, the scheme solves now
c     for u or ...
c     a further change is that fCO2_o is calculated interactively using
c         the observed DIC and Alk.
c
c     It includes several new features:
c
c       - new entrainment routine
c       - possibility to include role of advection 
c       - revised calculation scheme, C13 is handled as C13 and not
c          as dC13.
c       - possibility to do Monte Carlo simulations and 
c          parameter optimizations
c       - no common blocks, all explicit transfer
c       - each subroutine has an own file
c       - solving
c
c     REFERENCES:
c
c     Gruber, N. and C.D. Keeling, 1996. Seasonal carbon cycling in
c          the Sargasso Sea near Bermuda. Bulletin of the Scripps Institution
c          of Oceanography, 33, Univ. of Calif. Press, Berkeley, 1999.
c
c     Gruber, N., C.D. Keeling and T.F. Stocker, 1998. Carbon-13 constraints 
c          on the seasonal inorganic carbon budget at the BATS site in the 
c          northwestern Sargasso Sea, Deep-Sea Research I, 45, 673-717, 1998.
c
c     REVISIONS:
c
c     Vs   Date     Author   Remarks
c
c     1.0  30.07.96   ng     first implementation based on bermuda model
c           2.08.96   ng     still under development
c           7.08.96   ng     solved problem in net_comm_prod. for
c                              details see notebook and latex document
c                              diag_boxmod.tex
c           8.08.96   ng     included more features and matear entrainment
c     2.0  12.08.96   ng     start to include parameter optimization 
c          15.08.96   ng     parameter optimization is running, it was
c                              necessary to normalize the parameters in
c                              order to get good results
c          16.08.96   ng     included monte carlo simulation
c          21.08.96   ng     included sensitivities
c           2.09.96   ng     included output of rates and fluxes to file
c                              for mc analysis
c     3.0   2.09.96   ng     included inverse calculation of u and 
c                              ddC13/dx --> run_inv_simulation.F
c                              including monte carlo simulations
c           4.09.96   ng     included several dC13 fits
c     3.1   5.09.96   ng     included possibility to run model 
c                              backward, this is done to test the
c                              numerics of the model
c                              changed definition of cost
c                              function slightly
c           9.09.96   ng     included possibility to use variable
c                              entrainment length, included dlent_dh
c                              include more exact calculation of 
c                              model equation and gradients in model
c     3.2  20.11.96   ng     discovered mistake in gasexchange.F : calculation
c                              of fCO2
c     3.3  28.05.97   ng     found error in wanninkhof wind speed 
c                              dependency of pv (**(-0.5) instead of **(0.5))
c          11.06.97   ng     included adjusted dC13 fits
c
c     4.0  15.10.99   ng     adapted seasonal model for interannual runs
c          19.10.99   ng     is working ! included plotting of u_var
c     4.1  06.06.00   ng     included new berm-bats harmonic and spline
c                              fits
c          08.06.00   ng     adding anomalies output
c     4.2  29.08.01   ng     updating splines, switched years to
c                              4 digit, i.e. 1999
c     4.3  05.05.02   ng     updated harm coeff and splines to 
c                              new data (May 2002),
c                            added a few more plots (running avg for
c                              anomalies, etc)
c     4.4  09.06.02   ng     added running average to monte carlo runs
c               
c
c =====================================================================
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/diagboxmod.h" 1
c===================== include file diagboxmod.h =========================
c
c general variables of diagnostic box model
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c  2.9.96  ng        included inverse mode
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c               
c exp_name      c  s   name of experiment 
c notes         c  s   notes to add to experiment on output
c sim_mode      c  s   character specifying simulation mode
c                       <st>  : standard simulation
c                       <mc>  : monte carlo simulation
c                       <op>  : parameter optimization
c                       <sn>  : sensitivity mode
c                       <in>  : inverse calculation of u and dC13/dx
c sol_scheme    c  s    <std> : original standard, solving for dC/dt
c               c  s    <adv> : solve for u velocity
c               c  s    <phy> : solve for Kz + u
c               c  s    <kz>  : solve for Kz
c
c dc13_coeff    c  s   character specifying which dc13 coefficients are
c                        used
c                       <o1>  : original data with harmonic fit of 1st order
c                       <o2>  : original data with harmonic fit of 2nd order
c                       <a2>  : artificial curve with 2nd order harm fit
c def_par       l  s   logical specifying use of default parameters for run
c incl_adv      l  s   logical specifying whether advection is included
c                        or not
c plot_out      l  s   logical indicating that plots should be generated
c
c===========================================================================
c
      character exp_name*10,notes*30,dc13_coeff*2
      logical def_par,incl_adv,plot_out
      character sim_mode*2,sol_scheme*3
      
# 98 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/simulation.h" 1
c======================= include file simulation.h ==========================
c
c variables for controlling the simulation
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 15.10.99 ng        adapted for interannual variability runs
c
c variables:
c
c name          type   description
c-----------------------------------------------------------------------------
c 
c tstart        dp s   start time of simulation  [decimal year]
c tend          dp s   end time of simulation    [decimal year]
c ts            dp s   time step of simulation   [decimal year]
c nstepmax      dp s   maximum number of steps  
c nstep         i  s   number of steps
c time          dp s   actual time               [decimal year] 
c
c-------------------------------------------------------------------------
c
      double precision tstart,tend,ts,time
      integer nstep,nstepmax
c
      parameter (nstepmax = 10000)



# 99 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/gasex_params.h" 1
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
# 100 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/diffent_params.h" 1
c===================== include file diffent_params.h =========================
c
c parameter of vertical diffusion and entrainment
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 13.8.96  ng        changed l_* to sg_*
c 09.9.96  ng        included const_lent
c
c variables:
c
c name            type   description
c-----------------------------------------------------------------------------
c                 
c ent_scheme      c  s   character indicating entrainment scheme to be used
c                            <in> : integrated vertical gradient scheme (DEF)
c                            <mc> : mass conservation scheme
c                            <ep> : episodic event scheme
c                            <ma> : matear scheme
c                        the <in> scheme is the new scheme developed by ng
c                                      for station BATS
c                        the <mc> scheme is the scheme of Thomas Stocker
c                                      with the assumption that the 
c                                      concentrations at a specified depth
c                                      remain constant
c                        the <ep> scheme is the original scheme developed by
c                                      C.D. Keeling used in the station 'S'
c                                      study
c                        the <ma> scheme is the scheme used by Matear (1995)
c                                      in his model (uses an entrainment
c                                      length scale)
c h_coeff         c  s   character specifying which coefficients for
c                            mixed layer depth coefficients should be taken
c                            <b> : bats          <s> : station S  
c                            <c> : bats CTD      <a> : bats avg CTD 
c                 
c dt_ent_re      dp  s   recurrence period of entrainment (only important
c                            for the <ep> scheme
c const_Kz        l  s   logical indicating whether const Kz is used or
c                            seasonally varying Kz
c kz_const       dp  s   value of constant kz [m2 s-1]
c const_vertgrad  l  s   logical indicating whether const vertical gradients
c                            are used or seasonally varying
c dCdz_const     dp  s   value of constant dCdz [umol kg-1 m-1]
c D_dcdz         dp  s   correction term for dCdz [umol kg-1 m-1]
c dd13C_dC_bml   dp  s   ratio of ddC13/dz and dC/dz below mixed layer
c                           [per mil kg umol-1]
c diff_fact      dp  s   multiplication factor for vertical diffusion
c ent_fact       dp  s   multiplication factor for entrainment
c h_th           dp  s   depth at which thermocline values for mc entrainment
c                           scheme are taken
c c_th           dp  s   C concentration at h_th [umol kg-1]
c dC13_th        dp  s   dC13 concentration at h_th [per mil]
c 
c const_lent     dp  l   logical indicating use of a constant ent_length
c lent_const     dp  s   constant entrainment length scale 
c                                (for Matear entrainment) [m]
c dlent_dh       dp  s   slope of ent length vs h [no units]
c
c sg_ddC13_dC_bml dp  s   sigma uncertainty of dd13C_dC_bml (for OP and MC)
c sg_diff_fact    dp  s   sigma uncertainty of diff_fact (for OP and MC)
c sg_ent_fact             sigma uncertainty of ent_fact (for OP and MC)
c sg_kz_const     dp  s   sigma uncertainty of kz (for OP and MC)
c sg_lent_const   dp  s   sigma uncertainty of ent_length (for OP and MC)
c sg_D_dcdz       dp  s   sigma uncertainty of dcdz
c                
c=============================================================================
c
      character ent_scheme*2,h_coeff*1
      logical const_kz,const_vertgrad,const_lent
      double precision dt_ent_re
      double precision kz_const,dCdz_const,D_dcdz
      double precision ddC13_dC_bml,diff_fact,ent_fact
      double precision h_th,C_th,dC13_th
      double precision lent_const,dlent_dh
c
      double precision sg_ddC13_dC_bml
      double precision sg_ent_fact,sg_diff_fact
      double precision sg_kz_const,sg_lent_const,sg_D_dcdz
      double precision sg_dlent_dh
      
# 101 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/adv_params.h" 1
c===================== include file adv_params.h =========================
c
c parameters of horizontal advection
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 13.8.96  ng        changed l_* to sg_*
c 17.10.99 ng        added dC13_dh
c
c variables:
c
c name            type   description
c-----------------------------------------------------------------------------
c                 
c dC_dh          dp  s   horizontal C gradient       [umol kg-1 m-1] 
c ddC13_dh       dp  s   horizontal dC13 gradient    [per mil m-1]
c dC13_dh        dp  s   horizontal C13 gradient     [umol kg-1 m-1] (calc)
c u              dp  s   horizontal velocity         [m s-1]
c
c sg_dC_dh       dp  a   1-sigma uncertainty for dC_dh (for OP and MC)
c sg_ddC13_dh    dp  a   1-sigma uncertainty for ddC13_dh (for OP and MC)
c sg_u           dp  a   1-sigma uncertainty for u (for OP and MC)
c 
c=============================================================================
c
      double precision dC_dh,ddC13_dh,dC13_dh
      double precision u
      double precision sg_dC_dh,sg_ddC13_dh
      double precision sg_u
# 102 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/ncp_params.h" 1
c===================== include file ncp_params.h =========================
c
c parameters of net community production
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 13.8.96  ng        changed l_* to sg_*
c
c variables:
c
c name            type   description
c-----------------------------------------------------------------------------
c                 
c D_dC13_org      dp s   correction term for dC13 value of organic matter
c dc13_org        dp s   dc13 value of organic matter
c                         [per mil]
c sg_D_dc13_org   dp s   1 sigma uncertainty of D_dc13_org (for OP and MC)
c 
c=============================================================================
c
      double precision D_dc13_org
      double precision sg_D_dc13_org
      double precision dc13_org
# 103 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/periods.h" 1
c======================= include file periods.h ==============================
c
c variables for the periods that are investigated 
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 17.10.99 ng        adapted for interannual runs, periods are replaced with
c                       years
c
c variables:
c
c name          type   description
c-----------------------------------------------------------------------------
c 
c nyears        i  s   number of years to analyze
c year_start    dp a   start of analysis year
c 
c  
c------------------------------------------------------------------------
c
      integer nyearmax,nyears
      parameter(nyearmax = 20)
c
      double precision year_start
c
      data year_start / 0.0 /



      

# 104 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/processes.h" 1
c======================= include file processes.h ============================
c
c variables for the processes that are investigated 
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c
c variables:
c
c name          type   description
c-----------------------------------------------------------------------------
c 
c nproc_fluxes      i  s   number of processes for the fluxes
c nproc_rates       i  s   number of processes for the rates of change
c proc_fluxes_name  c  a   array holding the names of the flux processes
c proc_rates_name   c  a   array holding the names of the rates processes
c  
c------------------------------------------------------------------------
c
      integer nproc_fluxes,nproc_rates
      parameter(nproc_fluxes = 10, nproc_rates = 7)
c
      character proc_fluxes_name(nproc_fluxes)*4
      character proc_rates_name(nproc_rates)*4
c
      data proc_rates_name  /'ex  ','diff', 'ent ', 'adv ', 'ncp ',
     $                       'sim ','obs '/
      data proc_fluxes_name /'ex  ','diff', 'ent ', 'adv ', 'ncp ',
     $                       'sim ','obs ', 'zcmp', 'tcal', 'tobs'/



      

# 105 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/slabs.h" 1
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


# 106 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/anomalies.h" 1
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
# 107 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/runavg.h" 1
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
# 108 "diag_boxmodel_ia.F" 2


# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/intval.h" 1
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



# 110 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/obs_values.h" 1
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


# 111 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/statistics.h" 1
c==================== include file statistics.h ==============================
c
c variables for the run statistics
c
c REVISIONS:
c
c date     author    remarks
c
c 08.8.96  ng        first implementation
c 05.9.96  ng        included c_min,c_max,c_min_obs,c_max_obs
c
c variables:
c
c name          type   description
c-----------------------------------------------------------------------------
c 
c chisq         dp s   chi**2 between simulated and observed rate of change
c rsq           dp s   r**2 of fit between simulated and observed rates
c closure       dp s   annual lack of closure
c amplitude     dp s   simulated annual amplitude
c amplitude_obs dp s   observed annual amplitude
c c_min         dp s   simulated minimum C
c c_max         dp s   simulated maximum C
c c_min_obs     dp s   observed minimum C
c c_max_obs     dp s   observed maximum C
c par_sum       dp s   value of parameter misfit
c costfn        dp s   value of cost function
c chisq_term    dp s   value of chi**2 term in cost function
c clos_term     dp s   value of closure term in cost function
c ampl_term     dp s   value of amplitude term in cost function
c par_term      dp s   value of parameter term in cost function
c
c-------------------------------------------------------------------------
c
      double precision chisq,rsq,closure,costfn,amplitude
      double precision amplitude_obs,par_sum
      double precision c_min,c_max,c_min_obs,c_max_obs
      double precision chisq_term,clos_term,ampl_term,par_term



# 112 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/averages.h" 1
c======================= include file averages.h ==============================
c
c variables for calculating annual averages
c
c REVISIONS:
c
c date     author    remarks
c
c 31.07.96  ng        first implementation
c 10.09.96  ng        included dC13 values for different processes
c
c variables:
c
c name           type   description
c-----------------------------------------------------------------------------
c 
c sum_fract     dp  s   sum for calculating average fractionation factor
c sum_kex       dp  s   sum for calculating average kex
c sum_dCdz      dp  s   sum for calculating average vertical DIC gradient
c sum_ddC13dz   dp  s   sum for calculating average vertical dC13 gradient
c sum_dC13_prc  dp  a   sum for calculating average dC13 values of the proc
c nsum_dc13_prc  i  a   number of values used for summing sum_dc13_prc
c
c avg_fract     dp  s   average fractionation factor [per mil]
c dc13_prc      dp  a   average dc13 value of the processes [per mil]
c avg_kex       dp  s   average kex value [mol m-2 s-1 uatm-1]
c
c-------------------------------------------------------------------------
c
      double precision sum_fract,sum_kex,sum_dCdz,sum_ddC13dz
      double precision avg_fract,avg_kex
      double precision sum_dC13_prc(4)
      integer nsum_dc13_prc(4)
      double precision dc13_prc(4)



# 113 "diag_boxmodel_ia.F" 2
c

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/conj_grad.h" 1
c===================== include file conj_grad.h =========================
c
c variables for finding minimum using conjugate gradient technique
c
c REVISIONS:
c
c date     author    remarks
c
c 12.8.96  ng        first implementation
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c     
c assign         l  s   logical for subroutine assign_parameters
c                        true : assign values to parameters
c                        false: assign parameters to values   
c p             dp  a   array of parameters
c p_opt         dp  a   array of optimum parameters
c sg_p          dp  a   sigma of parameters
c p_init        dp  a   array of initial parameters
c nparams        i  s   dimension of parameters
c npmax          i  s   maximum dimension of p
c iter           i  s   number of iterations
c fret          dp  s   value of function at minimum
c hessian       dp  a   hessian matrix
c covm          dp  a   inverse of hessian matrix, covariance matrix
c
c sg_c          dp  s   1-sigma uncertainty of C [umol kg-1]
c sg_dcdt       sp  s   1-sigma uncertainty of dcdt [umol kg-1 day-1]
c                          (estimated from figure)
c
c===========================================================================
c
      logical assign
      double precision ftol,fret
      integer nruns
      integer nparams,iter


      parameter(nparams = 10)
# 53 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/conj_grad.h"
      double precision p(nparams),p_opt(nparams),p_init(nparams)
      double precision sg_p(nparams)
      double precision hessian(nparams,nparams),covm(nparams,nparams)
c
      
# 115 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/conj_grad_common.h" 1
c================ include file conj_grad_common.h =========================
c
c special include block for subroutine func.F. You find here all
c   parameters which are necessary to run the simulation, except for
c   the parameters that are optimized by the conjugate gradient 
c   technique. for variable definitions see *_params.h file
c
c REVISIONS:
c
c date     author    remarks
c
c 12.8.96  ng        first implementation
c  9.9.96  ng        included new definition of entrainment length
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c   
c
c===========================================================================
c
      common/sim/tstart,tend,ts,nstep
c
      common/ex/D_fco2_corr,gasex_fact,pv_relship
c
      common/coeff/dc13_coeff,h_coeff,ws_coeff
c
      common/diffent/
     $     dt_ent_re,lent_const,dlent_dh,
     $     ent_fact,diff_fact,
     $     kz_const,dcdz_const,
     $     h_th,c_th,dc13_th,ddc13_dc_bml,
     $     const_vertgrad,const_kz,const_lent,
     $     ent_scheme
c
      common/adv/dc_dh,ddc13_dh,u,incl_adv
c
      common/ncp/D_dc13_org
c
      common/conjgrad/sg_p,p_init
      
# 116 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/op_statistics.h" 1
c==================== include file op_statistics.h ==========================
c
c variables for the optimization run statistics
c
c REVISIONS:
c
c date     author    remarks
c
c 08.8.96  ng        first implementation
c
c variables:
c
c name           type   description
c-----------------------------------------------------------------------------
c 
c chisq_init     dp s   chi**2 between simulated and observed rate of change at
c                            the beginning of the optimization
c rsq_init       dp s   r**2 of fit between simulated and observed rates at
c                            the beginnning of the optimization
c closure_init   dp s   annual lack of closure at the beginning of the optim.
c amplitude_init dp s   annual amplitude at the beginning of the optimization
c costfn_init    dp s   value of cost function at the beginning of the optim.
c c_min_init     dp s   value of minimum C at the beginning of the optim.
c c_max_init     dp s   value of maximum C at the beginning of the optim.
c
c-------------------------------------------------------------------------
c
      double precision chisq_init,rsq_init,closure_init,costfn_init
      double precision amplitude_init,c_min_init,c_max_init



# 117 "diag_boxmodel_ia.F" 2
c

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/monte_carlo.h" 1
c===================== include file monte_carlo.h =========================
c
c general variables for monte carlo simulations
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c  
c nrunmax       i  s   maximal number of mc runs             
c nrun          i  s   number of mc runs
c
c *_r          dp  s   parameters with a random component added. size of
c                        random component is proportional to sg of this
c                        parameter. For the parameter definitions see
c                        the *_params.h files 
c 
c
c===========================================================================
c
      integer nrunmax,nrun
      parameter(nrunmax = 2000)
      
      double precision gasex_fact_r,D_fco2_corr_r
      double precision diff_fact_r,D_dcdz_r,ddC13_dC_bml_r,lent_const_r
      double precision dlent_dh_r
      double precision D_dc13_org_r
      double precision u_r,dc_dh_r,ddc13_dh_r
c
# 119 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/mc_statistics.h" 1
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

# 120 "diag_boxmodel_ia.F" 2

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/mc_plot_slabs.h" 1
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

# 121 "diag_boxmodel_ia.F" 2
c

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/sensitivity.h" 1
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
# 123 "diag_boxmodel_ia.F" 2
c

# 1 "/home/prafter/Gruber-Diagnostic-Box-Model/src_ia/incl/inverse_calc.h" 1
c================== include file inverse_calc.h =========================
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
c u_var           dp a   inverse calculated horizontal velocity m s-1
c 
c===========================================================================
c
      double precision u_var(nstepmax)
# 125 "diag_boxmodel_ia.F" 2
c
      logical error
      integer unit
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer nr,np,np2,seed,n,l,i,j,stdout
      parameter (stdout = 6)
      character filename*30
      double precision random,fret_best,costfn_min_old
c
      double precision int_fluxes_lr(2,nyearmax,nproc_fluxes)
      double precision int_rates_lr(2,nyearmax,nproc_rates)

c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      real ran1,gasdev
c
c =================== begin of executable code ========================
c
      error = .false.
c
c ---------------------------------------------------------------------
c     call introduction
c ---------------------------------------------------------------------
c
      call intro
c
c ---------------------------------------------------------------------
c     call program to initialize parameters
c ---------------------------------------------------------------------
c
      call init_params(incl_adv,plot_out,
     $     tstart,tend,ts,
     $     pv_relship,D_fco2_corr,sg_D_fco2_corr,
     $     gasex_fact,sg_gasex_fact,ws_coeff,
     $     h_coeff,ent_scheme,const_kz,kz_const,
     $     const_vertgrad,dCdz_const,D_dcdz,sg_D_dcdz,
     $     h_th,c_th,dc13_th,
     $     const_lent,lent_const,sg_lent_const,dlent_dh,sg_dlent_dh,
     $     ddC13_dC_bml,sg_ddC13_dC_bml,
     $     ent_fact,diff_fact,sg_diff_fact,sg_ent_fact,
     $     dC_dh,sg_dC_dh,ddC13_dh,sg_ddC13_dh,u,sg_u,
     $     D_dc13_org,sg_D_dc13_org) 
c
c ---------------------------------------------------------------------
c     get userinput
c ---------------------------------------------------------------------
c
      call userinput(exp_name,notes,
     $     sim_mode,sol_scheme,nrun,def_par,
     $     tstart,tend,ts,
     $     dc13_coeff,
     $     pv_relship,D_fco2_corr,gasex_fact,ws_coeff,
     $     comp_fco2_online,
     $     ent_scheme,h_coeff,
     $     h_th,c_th,dc13_th,dt_ent_re,
     $     const_lent,lent_const,dlent_dh,ent_fact,
     $     const_Kz,Kz_const,
     $     const_vertgrad,dCdz_const,D_dcdz,
     $     ddC13_dC_bml,diff_fact,
     $     incl_adv,dC_dh,ddC13_dh,u,
     $     D_dc13_org,
     $     plot_out)
c
c      write(*,*) 'notes : ',notes
c      write(*,*) 'exp_name : ',exp_name
c
c ---------------------------------------------------------------------
c     calculate nstep
c ---------------------------------------------------------------------
c
c      nstep = int((tend-tstart)/ts) + 1
      nstep = int((tend-tstart)/ts)
      write(*,*) '--MAIN: number of steps : ',nstep
c
c =====================================================================
c     run the simulation in standard mode
c =====================================================================
c
      if (sim_mode .eq. 'st') then
c
c ---------------------------------------------------------------------
c     call the simulation 
c ---------------------------------------------------------------------
c
         write(*,*) '--MAIN: running a standard simulation...'
         call run_simulation(sol_scheme,
     $         tstart,tend,ts,nstep,
     $         dc13_coeff,
     $         pv_relship,D_fco2_corr,gasex_fact,ws_coeff,
     $         comp_fco2_online,
     $         ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $         dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $         const_Kz,Kz_const,
     $         const_vertgrad,dCdz_const,D_dcdz,
     $         ddC13_dC_bml,diff_fact,
     $         incl_adv,dC_dh,ddC13_dh,u,u_var,
     $         D_dc13_org,
     $         D_c,D_c13,fluxes,
     $         c_obs,dc13_obs,h_obs,
     $         temp_obs,sal_obs,salk_obs,
     $         fco2_o_obs,pco2_a_obs,fco2_a_obs,pco2_a_dry,
     $         dc13_a_obs,ws_obs,dcdz_obs,kz_obs,
     $         c_sim,
     $         int_rates,int_fluxes,nyears,
     $         chisq,rsq,closure,amplitude,amplitude_obs,
     $         c_min,c_max,c_min_obs,c_max_obs,costfn,
     $         avg_fract,dc13_prc,avg_kex,
     $         error)
c
            if (error) goto 9900

c
c ---------------------------------------------------------------------
c     write the parameters and integrated results to screen and 
c       to a file
c ---------------------------------------------------------------------
c     
         do n = 1,2
c
            if (n .eq. 1) then
               unit = 6         ! output to screen
            elseif (n .eq. 2) then
               unit = 10
               filename = 'std_results.'//exp_name
               open (unit,file = filename, form = 'formatted', 
     $              status = 'unknown')
            endif
                        
            call write_std_output(unit,
     $           notes,exp_name,
     $           tstart,tend,ts, 
     $           sol_scheme,
     $           dc13_coeff,
     $           pv_relship,D_fco2_corr,gasex_fact,ws_coeff,   
     $           comp_fco2_online,
     $           ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $           dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $           const_Kz,Kz_const,
     $           const_vertgrad,dCdz_const,
     $           ddC13_dC_bml,diff_fact,
     $           incl_adv,dC_dh,ddC13_dh,u,
     $           D_dc13_org,
     $           int_rates,int_fluxes,nyears,
     $           chisq,rsq,closure,costfn,avg_fract,dc13_prc,avg_kex,
     $           error)
c
            if (n .eq. 2) then
               close(unit)
            endif
            if (error) then
               write(*,*) '--MAIN: an error occured during writing.'
               goto 9900
            endif
         enddo
c
c ---------------------------------------------------------------------
c     write the rates of changes to a file
c ---------------------------------------------------------------------
c     
         unit = 12
         filename = 'rates.'//exp_name
         open (unit,file = filename, form = 'formatted', 
     $        status = 'unknown')
c
         call write_rates_tofile(unit,tstart,ts,nstep,
     $        D_c,c_obs,c_sim,error)
         close(unit)
         if (error) then
            write(*,*) '--MAIN: an error occured during writing.'
            goto 9900
         endif
c
c ---------------------------------------------------------------------
c     write the fluxes to a file
c ---------------------------------------------------------------------
c     
         unit = 14
         filename = 'fluxes.'//exp_name
         open (unit,file = filename, form = 'formatted', 
     $        status = 'unknown')
c
         call write_fluxes_tofile(unit,tstart,ts,nstep,
     $        fluxes,error)
         close(unit)
         if (error) then
            write(*,*) '--MAIN: an error occured during writing.'
            goto 9900
         endif
c
c ---------------------------------------------------------------------
c     write the observations to a file
c ---------------------------------------------------------------------
c     
         unit = 16
         filename = 'observations.'//exp_name
         open (unit,file = filename, form = 'formatted', 
     $        status = 'unknown')
c
         call write_obs_tofile(unit,tstart,ts,nstep,
     $        c_obs,dc13_obs,h_obs,
     $        temp_obs,sal_obs,
     $        salk_obs,
     $        fco2_o_obs,pco2_a_obs,fco2_a_obs,pco2_a_dry,
     $        dc13_a_obs,
     $        ws_obs,dcdz_obs,kz_obs,
     $        error)
         close(unit)
         if (error) then
            write(*,*) '--MAIN: an error occured during writing.'
            goto 9900
         endif
c
c ---------------------------------------------------------------------
c     write the solutions to a file
c ---------------------------------------------------------------------
c     
         unit = 18
         filename = 'solution.'//exp_name
         open (unit,file = filename, form = 'formatted', 
     $        status = 'unknown')
c
         call write_solution_tofile(unit,tstart,ts,nstep,
     $        sol_scheme,u_var,
     $        error)
         close(unit)
         if (error) then
            write(*,*) '--MAIN: an error occured during writing.'
            goto 9900
         endif
c
c ---------------------------------------------------------------------
c     compute the rate and flux anomalies 
c ---------------------------------------------------------------------
c     
         call compute_anomalies(tstart,tend,ts,nstep,
     $        D_C,fluxes,
     $        D_C_avg,D_C_anom,fluxes_avg,fluxes_anom)
c
c ---------------------------------------------------------------------
c     compute the annually smoothed rates and fluxes 
c         i.e. 365 day running averages
c ---------------------------------------------------------------------
c     
         call compute_runavg(tstart,tend,ts,nstep,
     $        D_C,fluxes,
     $        D_C_anom,fluxes_anom,
     $        D_C_runavg,fluxes_runavg,
     $        D_C_a_runavg,fluxes_a_runavg)
c
c ---------------------------------------------------------------------
c     write the flux anomalies to a file
c ---------------------------------------------------------------------
c     
         unit = 20
         filename = 'flux_anomalies.'//exp_name
         open (unit,file = filename, form = 'formatted', 
     $        status = 'unknown')
c
         call write_fluxanomalies_tofile(unit,tstart,ts,nstep,
     $        fluxes_avg,fluxes_anom,fluxes_runavg,
     $        error)
         close(unit)
         if (error) then
            write(*,*) '--MAIN: an error occured during writing.'
            goto 9900
         endif
c
c ---------------------------------------------------------------------
c     write the rate anomalies to a file
c ---------------------------------------------------------------------
c     
         unit = 22
         filename = 'rate_anomalies.'//exp_name
         open (unit,file = filename, form = 'formatted', 
     $        status = 'unknown')
c
         call write_rateanomalies_tofile(unit,tstart,ts,nstep,
     $        D_C_avg,D_C_anom,D_C_runavg,
     $        error)
         close(unit)
         if (error) then
            write(*,*) '--MAIN: an error occured during writing.'
            goto 9900
         endif
c
c ---------------------------------------------------------------------
c     generate NCAR plots if desired
c ---------------------------------------------------------------------
c     
         if (plot_out) then
            write(stdout,*) '--DIAG_BOXMODEL: generating NCAR plots...'
            call generate_ncar_plots(tstart,tend,ts,nstep,
     $           D_C,fluxes,
     $           D_C_avg,D_C_anom,fluxes_avg,fluxes_anom,
     $           D_C_runavg,fluxes_runavg,
     $           D_C_a_runavg,fluxes_a_runavg,
     $           c_sim,c_obs,u_var,exp_name,notes)
         endif        
c
c =====================================================================
c     run the simulation in monte carlo mode
c =====================================================================
c
      elseif (sim_mode .eq. 'mc')  then
c
         write(*,*) '--MAIN: running a Monte Carlo simulation...'
c
         seed = -23
c
c ---------------------------------------------------------------------
c     calculate the standard simulation 
c ---------------------------------------------------------------------
c
         write(*,*) '--MAIN: running the standard simulation...'
c
         call run_simulation(sol_scheme,
     $         tstart,tend,ts,nstep,
     $         dc13_coeff,
     $         pv_relship,D_fco2_corr,gasex_fact,ws_coeff,
     $         comp_fco2_online,
     $         ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $         dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $         const_Kz,Kz_const,
     $         const_vertgrad,dCdz_const,D_dcdz,
     $         ddC13_dC_bml,diff_fact,
     $         incl_adv,dC_dh,ddC13_dh,u,u_var,
     $         D_dc13_org,
     $         D_c_std,D_c13,fluxes_std,
     $         c_obs,dc13_obs,h_obs,
     $         temp_obs,sal_obs,salk_obs,
     $         fco2_o_obs,pco2_a_obs,fco2_a_obs,pco2_a_dry,
     $         dc13_a_obs,ws_obs,dcdz_obs,kz_obs,
     $         c_sim_std,
     $         int_rates_std,int_fluxes_std,nyears,
     $         chisq_std,rsq_std,closure_std,
     $         amplitude_std,amplitude_obs,
     $         c_min_std,c_max_std,c_min_obs,c_max_obs,costfn_std,
     $         avg_fract,dc13_prc,avg_kex,error)
c
            if (error) goto 9900
c
c ---------------------------------------------------------------------
c     keep the inital parameters in p_init , this
c      is necessary to calculate the cost function
c ---------------------------------------------------------------------
c
         assign = .true.
         call assign_params(assign,
     $        gasex_fact, sg_gasex_fact,
     $        D_fco2_corr,sg_D_fco2_corr,
     $        diff_fact,sg_diff_fact,D_dcdz,sg_D_dcdz,
     $        ddC13_dC_bml,sg_ddC13_dC_bml,
     $        const_lent,lent_const,sg_lent_const,
     $        dlent_dh,sg_dlent_dh,   
     $        D_dc13_org,sg_D_dc13_org, 
     $        u,sg_u,dc_dh,sg_dc_dh,ddC13_dh,sg_ddC13_dh,   
     $        p_init,sg_p,p_init)   
c
c ---------------------------------------------------------------------
c     loop over all monte carlo runs
c ---------------------------------------------------------------------
c
         do nr = 1,nrun
c
            if (mod(nr,1) .eq. 0) then
               write(unit,*) '--MAIN: running mc simulation # ',nr
            endif
c
c ---------------------------------------------------------------------
c     add a normal distributed random component to the parameters
c ---------------------------------------------------------------------
c
            gasex_fact_r = gasex_fact + sg_gasex_fact * 
     $           dble(gasdev(seed))
            D_fco2_corr_r = D_fco2_corr + sg_D_fco2_corr *
     $           dble(gasdev(seed))
            diff_fact_r = diff_fact + sg_diff_fact *
     $           dble(gasdev(seed))
            D_dcdz_r = D_dcdz + sg_D_dcdz *
     $           dble(gasdev(seed))
            ddC13_dC_bml_r = ddC13_dC_bml + sg_ddC13_dC_bml *
     $           dble(gasdev(seed))
            if (const_lent) then
               lent_const_r = lent_const + sg_lent_const *
     $              dble(gasdev(seed))
            else
               dlent_dh_r = dlent_dh + sg_dlent_dh *
     $              dble(gasdev(seed))
            endif
            D_dc13_org_r = D_dc13_org + sg_D_dc13_org *
     $           dble(gasdev(seed))
            u_r = u + sg_u * dble(gasdev(seed))

            dc_dh_r = dc_dh + sg_dc_dh *
     $           dble(gasdev(seed))
            ddc13_dh_r = ddc13_dh + sg_ddc13_dh *
     $           dble(gasdev(seed))




c
            assign = .true.
            call assign_params(assign,
     $           gasex_fact_r, sg_gasex_fact,
     $           D_fco2_corr_r,sg_D_fco2_corr,
     $           diff_fact_r,sg_diff_fact,
     $           D_dcdz_r,sg_D_dcdz,
     $           ddC13_dC_bml_r,sg_ddC13_dC_bml,
     $           const_lent,lent_const_r,sg_lent_const,
     $           dlent_dh_r,sg_dlent_dh,   
     $           D_dc13_org_r,sg_D_dc13_org, 
     $           u_r,sg_u,dc_dh_r,sg_dc_dh,ddC13_dh_r,sg_ddC13_dh,   
     $           p,sg_p,p_init)
c
c
            write(*,'(/A/)') ' monte carlo parameters : '
c     
            write(*,'(4(a9,1x)/)') 
     $           'parameter',
     $           'initial  ',
     $           '1-sg     ',
     $           'current  '
            
            do n = 1,nparams
               write(*,'(5x,i2,3x,f8.4,3x,f8.4,3x,f8.4)') 
     $              n,p_init(n),sg_p(n),p(n)
            enddo
c
c
c ---------------------------------------------------------------------
c     call the simulation 
c ---------------------------------------------------------------------
c
            call run_simulation(sol_scheme,
     $            tstart,tend,ts,nstep,
     $            dc13_coeff,
     $            pv_relship,D_fco2_corr_r,gasex_fact_r,ws_coeff,
     $            comp_fco2_online,
     $            ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $            dt_ent_re,const_lent,lent_const_r,
     $            dlent_dh_r,ent_fact,
     $            const_Kz,Kz_const,
     $            const_vertgrad,dCdz_const,D_dcdz_r,
     $            ddC13_dC_bml_r,diff_fact_r,
     $            incl_adv,dC_dh_r,ddC13_dh_r,u_r,u_var,
     $            D_dc13_org_r,
c     
     $            D_c,D_c13,fluxes,
     $            c_obs,dc13_obs,h_obs,
     $            temp_obs,sal_obs,salk_obs,
     $            fco2_o_obs,pco2_a_obs,fco2_a_obs,pco2_a_dry,
     $            dc13_a_obs,ws_obs,dcdz_obs,kz_obs,
     $            c_sim,
     $            int_rates,int_fluxes,nyears,
     $            chisq,rsq,closure,amplitude,amplitude_obs,
     $            c_min,c_max,c_min_obs,c_max_obs,costfn,
     $            avg_fract,dc13_prc,avg_kex,error)
c
            if (error) goto 9900
c
c ---------------------------------------------------------------------
c     keep certain information of each run
c ---------------------------------------------------------------------
c
c            write(*,*) '--MAIN: keep statistics...'
c
            costfn_min_old = costfn_min
c
            call keep_mc_statistics(
     $           nr,nstep,
     $           D_c,fluxes,c_sim,
     $           int_rates,int_fluxes,nyears,
     $           int_rates_std,int_fluxes_std,costfn_std,
     $           chisq,closure,amplitude,amplitude_obs,costfn,
     $           p,sg_p,p_init,
     $           mc_D_c,mc_fluxes,mc_c_sim,mc_weight,
     $           nstep2,rsmpl_step,
     $           var_int_rates,var_int_fluxes,
     $           wvar_int_rates,wvar_int_fluxes,
     $           sum_costfn,costfn_min,costfn_max)
c
c ---------------------------------------------------------------------
c     if a new minimum has been found then keep the minimum
c       parameters
c ---------------------------------------------------------------------
c
            if (costfn_min .ne. costfn_min_old) then
               do np = 1,nparams
                  p_opt(np) = p(np)
               enddo
            endif
c     
c ---------------------------------------------------------------------
c     end of the loop over all monte carlo runs
c ---------------------------------------------------------------------
c
         enddo
c
c ---------------------------------------------------------------------
c     calculate monte-carlo statistics
c ---------------------------------------------------------------------
c
         write(*,*) '--MAIN: calculating mc statistics...'
c
         call calc_mc_statistics(nrun,     
     $        int_rates_std,int_fluxes_std,nyears,costfn_std,
     $        var_int_rates,var_int_fluxes,
     $        wvar_int_rates,wvar_int_fluxes,
     $        sum_costfn,
     $        mc_D_c,mc_fluxes,mc_c_sim,mc_weight,nstep2,
     $        sg_int_rates,sg_int_fluxes,
     $        wsg_int_rates,wsg_int_fluxes,
     $        mean_costfn,
     $        mean_d_c,mean_fluxes,mean_c_sim,
     $        wmean_d_c,wmean_fluxes,wmean_c_sim,
     $        sg_d_c,sg_fluxes,sg_c_sim,
     $        wsg_d_c,wsg_fluxes,wsg_c_sim)
c
c ---------------------------------------------------------------------
c     write the parameters and integrated results of the
c       monte carlo run to screen and to a file
c ---------------------------------------------------------------------
c     
         do n = 1,2
c
            if (n .eq. 1) then
               unit = 6         ! output to screen
            elseif (n .eq. 2) then
               unit = 10
               filename = 'mc_results.'//exp_name
               open (unit,file = filename, form = 'formatted', 
     $              status = 'unknown')
            endif
                        
            call write_mc_output(unit,
     $           notes,exp_name,
     $           tstart,tend,ts,sol_scheme,
     $           dc13_coeff,
     $           pv_relship,D_fco2_corr,gasex_fact, ws_coeff,         
     $           ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $           dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $           const_Kz,Kz_const,
     $           const_vertgrad,dCdz_const,
     $           ddC13_dC_bml,diff_fact,
     $           incl_adv,dC_dh,ddC13_dh,u,
     $           D_dc13_org,
c          
     $           nrun,
     $           p_init,p_opt,sg_p,
c
     $           int_rates_std,int_fluxes_std,nyears,
     $           sg_int_rates,sg_int_fluxes,
     $           wsg_int_rates,wsg_int_fluxes,
     $           chisq_std,rsq_std,closure_std,
     $           costfn_std,mean_costfn,costfn_min,costfn_max,
     $           error)
c
            if (n .eq. 2) then
               close(unit)
            endif
            if (error) then
               write(*,*) '--MAIN: an error occured during writing.'
               goto 9900
            endif
         enddo
c
c ---------------------------------------------------------------------
c     generate NCAR plots
c ---------------------------------------------------------------------
c     
         if (plot_out) then
c
            write(*,*) '--MAIN: plotting monte carlo results...'
c
            call generate_mc_ncar_plots(nrun,
     $           tstart,tend,ts,nstep,
     $           D_C_std,fluxes_std,c_sim_std,c_obs,
     $           mean_d_c,mean_fluxes,mean_c_sim,
     $           wmean_d_c,wmean_fluxes,wmean_c_sim,
     $           sg_D_c,sg_fluxes,sg_c_sim,
     $           wsg_D_c,wsg_fluxes,wsg_c_sim,
     $           nstep2,rsmpl_step,
     $           exp_name,notes)
         endif        
c
c ---------------------------------------------------------------------
c     write the rates of changes to a file
c ---------------------------------------------------------------------
c     
         unit = 12
         filename = 'mc_rates.'//exp_name
         open (unit,file = filename, form = 'formatted', 
     $        status = 'unknown')
c
         call write_mc_rates_tofile(unit,tstart,ts,
     $        nstep2,rsmpl_step,
     $        D_c_std,c_sim_std,c_obs,
     $        mean_d_c,mean_c_sim,
     $        wmean_d_c,wmean_c_sim,
     $        sg_D_c,sg_c_sim,
     $        wsg_D_c,wsg_c_sim,
     $        error)
c
         close(unit)
         if (error) then
            write(*,*) '--MAIN: an error occured during writing.'
            goto 9900
         endif
c
c ---------------------------------------------------------------------
c     write the fluxes to a file
c ---------------------------------------------------------------------
c     
         unit = 14
         filename = 'mc_fluxes.'//exp_name
         open (unit,file = filename, form = 'formatted', 
     $        status = 'unknown')
c
         call write_mc_fluxes_tofile(unit,tstart,ts,
     $        nstep2,rsmpl_step,
     $        fluxes_std,
     $        mean_fluxes,wmean_fluxes,sg_fluxes,wsg_fluxes,
     $        error)
         close(unit)
         if (error) then
            write(*,*) '--MAIN: an error occured during writing.'
            goto 9900
         endif
c
c ---------------------------------------------------------------------
c     compute the annually smoothed rates and fluxes 
c         i.e. 365 day running averages
c ---------------------------------------------------------------------
c     
         call compute_mc_runavg(
     $     tstart,tend,ts,nstep, nstep2,
     $     rsmpl_step, D_c_std, mean_d_c,sg_D_c, fluxes_std,
     $     mean_fluxes,sg_fluxes, D_C_std_runavg, mean_d_c_runavg,
     $     sg_D_c_runavg, fluxes_std_runavg, mean_fluxes_runavg,
     $     sg_fluxes_runavg,error)
c
c ---------------------------------------------------------------------
c     write the flux anomalies to a file
c ---------------------------------------------------------------------
c     
         unit = 20
         filename = 'flux_mc_anomalies.'//exp_name
         open (unit,file = filename, form = 'formatted', 
     $        status = 'unknown')
c
         write(*,*) '--MAIN: saving running averages to file...'
c
         call write_mc_fluxanomalies_tofile(
     $        unit,tstart,ts,nstep,nstep2,rsmpl_step,
     $        fluxes_std_runavg, mean_fluxes_runavg,sg_fluxes_runavg,
     $        error)

         close(unit)
         if (error) then
            write(*,*) '--MAIN: an error occured during writing.'
            goto 9900
         endif
c
c =====================================================================
c     run the simulation in sensitivity mode
c =====================================================================
c
      elseif (sim_mode .eq. 'sn') then
c
         write(*,*) '--MAIN: determining sensitivities...'
c
c ---------------------------------------------------------------------
c     calculate the standard simulation 
c ---------------------------------------------------------------------
c
         write(*,*) '--MAIN: running the standard simulation...'
c
         call run_simulation(sol_scheme,
     $        tstart,tend,ts,nstep,
     $        dc13_coeff,
     $        pv_relship,D_fco2_corr,gasex_fact,ws_coeff,  
     $        comp_fco2_online,  
     $        ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $        dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $        const_Kz,Kz_const,
     $        const_vertgrad,dCdz_const,D_dcdz,
     $        ddC13_dC_bml,diff_fact,
     $        incl_adv,dC_dh,ddC13_dh,u,u_var,
     $        D_dc13_org,
     $        D_c_std,D_c13,fluxes_std,
     $        c_obs,dc13_obs,h_obs,
     $        temp_obs,sal_obs,salk_obs,
     $        fco2_o_obs,pco2_a_obs,fco2_a_obs,pco2_a_dry,
     $        dc13_a_obs,ws_obs,dcdz_obs,kz_obs,
     $        c_sim_std,
     $        int_rates_std,int_fluxes_std,nyears,
     $        chisq_std,rsq_std,closure_std,
     $        amplitude_std,amplitude_obs,
     $        c_min_std,c_max_std,c_min_obs,c_max_obs,costfn_std,
     $        avg_fract,dc13_prc,avg_kex,error)
c
c
         if (error) goto 9900
         write(*,*) '--MAIN: standard simulation completed'
c
c ---------------------------------------------------------------------
c     assign the individual parameters to the p_init arrays
c ---------------------------------------------------------------------
c
         assign = .true.
         call assign_params(assign,
     $        gasex_fact, sg_gasex_fact,
     $        D_fco2_corr,sg_D_fco2_corr,
     $        diff_fact,sg_diff_fact,D_dcdz,sg_D_dcdz,
     $        ddC13_dC_bml,sg_ddC13_dC_bml,
     $        const_lent,lent_const,sg_lent_const,
     $        dlent_dh,sg_dlent_dh,   
     $        D_dc13_org,sg_D_dc13_org, 
     $        u,sg_u,dc_dh,sg_dc_dh,ddC13_dh,sg_ddC13_dh,   
     $        p_init,sg_p,p_init) 
c
c ---------------------------------------------------------------------
c     loop over all parameters
c ---------------------------------------------------------------------
c
         do np = 1,nparams
c
c ---------------------------------------------------------------------
c     loop over p+sg_p and p-sg_p
c ---------------------------------------------------------------------
c
            do l = 1,2
c
c ---------------------------------------------------------------------
c     determine parameters
c ---------------------------------------------------------------------
c
               do np2 = 1,nparams
                  p(np2) = p_init(np2)
               enddo
               if (l .eq. 1) then
                  p(np) = p_init(np) - sg_p(np)
               elseif (l .eq. 2) then
                  p(np) = p_init(np) + sg_p(np)
               endif
c     
               gasex_fact   = p(1)
               D_fco2_corr  = p(2)
               diff_fact    = p(3)
               D_dcdz       = p(4)
               ddC13_dC_bml = p(5)
               if (const_lent) then
                  lent_const   = p(6)
               else
                  dlent_dh     = p(6)
               endif
               D_dc13_org   = p(7)
               u            = p(8)

               dc_dh        = p(9)
               ddC13_dh     = p(10)

c
c ---------------------------------------------------------------------
c     run the simulation
c ---------------------------------------------------------------------
c
c         write(*,*) '--MAIN: sensitivity run for param # ',np
c
               call run_simulation(sol_scheme,
     $              tstart,tend,ts,nstep,
     $              dc13_coeff,
     $              pv_relship,D_fco2_corr,gasex_fact,ws_coeff,   
     $              comp_fco2_online,    
     $              ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $              dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $              const_Kz,Kz_const,
     $              const_vertgrad,dCdz_const,D_dcdz,
     $              ddC13_dC_bml,diff_fact,
     $              incl_adv,dC_dh,ddC13_dh,u,u_var,
     $              D_dc13_org,
     $              D_c,D_c13,fluxes,
     $              c_obs,dc13_obs,h_obs,
     $              temp_obs,sal_obs,salk_obs,
     $              fco2_o_obs,pco2_a_obs,fco2_a_obs,pco2_a_dry,
     $              dc13_a_obs,ws_obs,dcdz_obs,kz_obs,
     $              c_sim,
     $              int_rates,int_fluxes,nyears,
     $              chisq,rsq,closure,
     $              amplitude,amplitude_obs,
     $              c_min,c_max,c_min_obs,c_max_obs,costfn,
     $              avg_fract,dc13_prc,avg_kex,error)

c
               if (error) goto 9900
c         write(*,*) '--MAIN: sensitivity run completed '
c     
c ---------------------------------------------------------------------
c     keep the int_rates and int_fluxes for determining the
c      sensitivity
c ---------------------------------------------------------------------
c
               do i = 1,nyears
                  do j = 1,nproc_rates
                     int_rates_lr(l,i,j) = int_rates(i,j)
                  enddo
                  do j = 1,nproc_fluxes
                     int_fluxes_lr(l,i,j) = int_fluxes(i,j)
                  enddo
               enddo
c
c ---------------------------------------------------------------------
c     end of loop over p + sg_p and p - sg_p
c ---------------------------------------------------------------------
c
            enddo
c
c ---------------------------------------------------------------------
c     determine sensitivities; this is done using centered differences
c ---------------------------------------------------------------------
c
            do i = 1,nyears
               do j = 1,nproc_rates
                  d_int_rates_dp(np,i,j) = 0.5d0/sg_p(np) *
     $                 (int_rates_lr(2,i,j) - int_rates_lr(1,i,j))
               enddo
               do j = 1,nproc_fluxes
                  d_int_fluxes_dp(np,i,j) = 0.5d0/sg_p(np) *
     $                 (int_fluxes_lr(2,i,j) - int_fluxes_lr(1,i,j))
               enddo
            enddo
c
c ---------------------------------------------------------------------
c     end of loop over parameters
c ---------------------------------------------------------------------
c
         enddo
c
c ---------------------------------------------------------------------
c     reassign the original parameters
c ---------------------------------------------------------------------
c
         gasex_fact  =  p_init(1) 
         D_fco2_corr =  p_init(2) 
         diff_fact   =  p_init(3) 
         D_dcdz      =  p_init(4) 
         ddC13_dC_bml=  p_init(5) 
         if (const_lent) then
            lent_const  =  p_init(6)
         else
            dlent_dh    =  p_init(6)
         endif
         D_dc13_org  =  p_init(7) 
         u           =  p_init(8) 

         dc_dh       =  p_init(9) 
         ddC13_dh    =  p_init(10)

c
c ---------------------------------------------------------------------
c     write results of sensitivity studies to screen and file
c ---------------------------------------------------------------------
c
         do n = 1,2
c
            if (n .eq. 1) then
               unit = 6         ! output to screen
            elseif (n .eq. 2) then
               unit = 10
               filename = 'sens_results.'//exp_name
               open (unit,file = filename, form = 'formatted', 
     $              status = 'unknown')
            endif
                        
            call write_sn_output(unit,
     $           notes,exp_name,
     $           tstart,tend,ts,
     $           dc13_coeff,
     $           pv_relship,D_fco2_corr,gasex_fact,ws_coeff,          
     $           ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $           dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $           const_Kz,Kz_const,
     $           const_vertgrad,dCdz_const,
     $           ddC13_dC_bml,diff_fact,
     $           incl_adv,dC_dh,ddC13_dh,u,
     $           D_dc13_org,
c          
     $           p_init,sg_p,
c
     $           int_rates_std,int_fluxes_std,nyears,
     $           d_int_rates_dp,d_int_fluxes_dp,
     $           chisq_std,rsq_std,closure_std,
     $           costfn_std,
     $           error)
c
            close(unit)
            if (error) then
               write(*,*) '--MAIN: an error occured during writing.'
               goto 9900
            endif
         enddo
c
c =====================================================================
c     end of sensitivity mode
c =====================================================================
c
      endif
c
c ---------------------------------------------------------------------
c     program terminated normally
c ---------------------------------------------------------------------
c
      write(*,9020)
 9020 format(/'program done. That was all, folks !'/
     $        '     see also gmeta file for output '//)
      stop
c
c ---------------------------------------------------------------------
c     program terminated abnormally
c ---------------------------------------------------------------------
c
 9900 write(*,9030)
 9030 format(/'program exited abnormally !')
      stop


      end
