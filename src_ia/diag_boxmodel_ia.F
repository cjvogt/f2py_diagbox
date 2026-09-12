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
# include "diagboxmod.h"
# include "simulation.h"
# include "gasex_params.h"
# include "diffent_params.h"
# include "adv_params.h"
# include "ncp_params.h"
# include "periods.h"
# include "processes.h"
# include "slabs.h"
# include "anomalies.h"
# include "runavg.h"

# include "intval.h"
# include "obs_values.h"
# include "statistics.h"
# include "averages.h"
c
# include "conj_grad.h"
# include "conj_grad_common.h"
# include "op_statistics.h"
c
# include "monte_carlo.h"
# include "mc_statistics.h"
# include "mc_plot_slabs.h"
c
# include "sensitivity.h"
c
# include "inverse_calc.h"
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
#ifdef with_horiz_grad
            dc_dh_r = dc_dh + sg_dc_dh *
     $           dble(gasdev(seed))
            ddc13_dh_r = ddc13_dh + sg_ddc13_dh *
     $           dble(gasdev(seed))
#else
            dc_dh_r = dc_dh
            ddc13_dh_r = ddc13_dh
#endif
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
#ifdef with_horiz_grad
               dc_dh        = p(9)
               ddC13_dh     = p(10)
#endif
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
#ifdef with_horiz_grad
         dc_dh       =  p_init(9) 
         ddC13_dh    =  p_init(10)
#endif
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
