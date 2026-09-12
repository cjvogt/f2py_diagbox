c =====================================================================
c
c     SUBROUTINE RUN_SIMULATION
c 
c     PURPOSE: This is the main routine for doing the diagnostic
c                seasonal calculations in the mixed layer
c
c     VARIABLES:
c 
c     IN : sol_scheme      : solution scheme to be used
c          tstart          : start time of integration
c          tend            : end time of integration
c          ts              : time step
c          nstep           : number of steps
c          dc13_coeff      : character specifying which dC13 coefficients
c                              are used                         
c          pv_relship      : character indicating the relationship
c                               that is used for calculating the
c                               piston-velocity
c                              <lm> : Liss and Merlivat
c                              <wa> : Wanninkhof
c          D_fco2_corr     : correction term for mixed layer fco2
c          gasex_fact      : air-sea exchange multiplication factor
c                              (for lm usually 1.7447)
c          ws_coeff        : character specifying wind speed coefficient
c          comp_foc2_online: logical indicating computation of pCO2
c          ent_scheme      : character indicating entrainment scheme
c                              to be used
c                              <in> : integrated vertical gradient scheme
c                              <mc> : mass conservation scheme
c                              <ep> : episodic event scheme
c          h_coeff         : character specifying which h coeff are taken
c          h_th            : depth at which thermocline values for 
c                                 mc entrainment are taken
c          c_th            : C concentration at h_th [umol kg-1]
c          dC13_th         : dC13 concentration at h_th [per mil]
c          dt_ent_re       : recurrence time for entrainment
c          const_lent      : logical indicating use of ent_length
c          lent_const      : const entrainment length scale (for ma scheme)
c          dlent_dh        : ratio of lent to mixed layer depth []
c          ent_fact        : entrainment factor
c          const_Kz        : logical indicating if constant Kz or
c                              seasonally varying Kz is used
c          Kz_const        : value of constant Kz 
c          const_vertgrad  : logical indicating if constant or variable
c                              vertical gradient is used
c          dCdz_const      : constant vertical C gradient 
c          D_dcdz          : correction term for dc/dz
c          ddC13_dC_bml    : ratio of ddC13/dz to dC/dz below mixed layer
c          diff_fact       : vertical diffusion mulitplication factor
c          incl_adv        : logical indicating if effect of advection
c                              is included
c          dC_dh           : horizontal gradient of C
c          ddC13_dh        : horizontal gradient of dC13
c          u               : horizontal velocity
c          D_dc13_org      : change in dC13 value of organic matter
c                          
c     OUT: D_c             : rates of change of the different processes
c          fluxes          : fluxes of the different processes
c          D_c13           : rates of change of C13 of the different proc 
c          u_var           : calculated velocity
c          c_obs           : observed seasonal cycle of C
c          dc13_obs        : observed seasonal cycle of dC13
c          h_obs           : observed seasonal cycle of h
c          temp_obs        : observed seasonal cycle of temp
c          sal_obs         : observed seasonal cycle of sal
c          salk_obs        : observed seasonal cycle of salk
c          fco2_o_obs      : observed seasonal cycle of oceanic fco2
c          pco2_a_obs      : observed seasonal cycle of atmospheric pco2
c          pco2_a_dry      : observed seasonal cycle of atmospheric pco2
c          fco2_a_obs      : observed seasonal cycle of atmospheric fco2
c          dc13_a_obs      : observed seasonal cycle of atmospheric dc13
c          ws_obs          : observed seasonal cycle of windspeed
c          dcdz_obs        : observed seasonal cycle of dcdz
c          kz_obs          : observed seasonal cycle of Kz
c
c          c_sim           : simulated seasonal cycle of c
c          int_rates       : integrated rates of change
c          int_fluxes      : integrated fluxes
c          nyears          : number of analysis years
c
c          chisq           : chi**2 between simulated and obs rates
c          rsq             : r**2 between simulated and obs rates
c          closure         : annual lack of closure
c          amplitude       : simulated amplitude
c          amplitude_obs   : observed amplitude
c          c_min           : value of minimum C
c          c_max           : value of maximum C
c          c_min_obs       : value of observed minimum C
c          c_max_obs       : value of observed maximum C
c          costfn          : cost function (without parameter terms)
c          avg_fract       : average fractionation factor
c          dc13_prc        : dc13 value of processes
c          avg_kex         : average kex
c
c     REVISION:
c
c     date     author  remarks
c   
c     31.07.96   ng    first implementation
c     08.08.96   ng    included matear entrainment, and more seasonal
c                         observations
c     14.08.96   ng    included more statistics
c      5.09.96   ng    included even more statistics
c                ng    included possibility to do backward simulation
c     15.10.99   ng    adapted for interannual variability run, ATTENTION
c                        time units have changed from d-1 to yr-1
c     07.06.00   ng    added salk_obs
c
c =====================================================================
c
      subroutine run_simulation(sol_scheme,
     $     tstart,tend,ts,nstep,
     $     dc13_coeff,
     $     pv_relship,D_fco2_corr,gasex_fact,ws_coeff,
     $     comp_fco2_online, 
     $     ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $     dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $     const_Kz,kz_const,
     $     const_vertgrad,dCdz_const,D_dcdz,
     $     ddC13_dC_bml,diff_fact,
     $     incl_adv,dC_dh,ddC13_dh,u,u_var,
     $     D_dc13_org,
     $     D_c,D_c13,fluxes,
     $     c_obs,dc13_obs,h_obs,
     $     temp_obs,sal_obs,salk_obs,
     $     fco2_o_obs,pco2_a_obs,fco2_a_obs,pco2_a_dry,
     $     dc13_a_obs,ws_obs,dcdz_obs,kz_obs,
     $     c_sim,
     $     int_rates,int_fluxes,nyears,
     $     chisq,rsq,closure,amplitude,amplitude_obs,
     $     c_min,c_max,c_min_obs,c_max_obs,costfn,
     $     avg_fract,dc13_prc,avg_kex,
     $     error)
c
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
      implicit none
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "periods.h"
# include "slabs.h"
c
# include "gasex_params.h"
# include "diffent_params.h"
# include "adv_params.h"
# include "ncp_params.h"
c
# include "intval.h"
# include "averages.h"
# include "input.h"
# include "solution.h"
# include "obs_values.h"
# include "c13_params.h"
# include "statistics.h"
c
      logical error
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision rates_to_fluxes,c13_to_delta
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision rho,sec_day,day_year
      parameter (rho = 1026.2d0,sec_day = 24.d0 * 3600.d0)
      parameter (day_year = 365.)
      integer k,i
      double precision rstar,dc13_i(5),D_dc13(5)
      double precision dc13_new,dc13_old,dc13
      double precision d_dc13_calc,D_c13_calc
      double precision dc13_tmp,r_prc
      double precision sum,sum_c13
      double precision ratio,denom,denom2
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
CC 
 
      write(*,*) '--RUN_SIMULATION: starting the simulation...'

      write(*,*) '     start time of integration : ',tstart
      write(*,*) '     end time of integration   : ',tend
      write(*,*) '     time step                 : ',ts
      write(*,*) '     # of steps                : ',nstep
      write(*,*) '     ent_scheme                : ',ent_scheme

c
c --------------------------------------------------------------------
c     initialize all integrated values and sums (set them to zero) 
c         especially important for Monte Carlo simulations
c --------------------------------------------------------------------
c       
      call set_sums_to_zero(
     $     sum_fract,sum_kex,
     $     sum_dCdz, sum_ddC13dz,
     $     sum_dc13_prc,nsum_dc13_prc)
c
c ---------------------------------------------------------------------
c     calculate initial values
c ---------------------------------------------------------------------
c
      time = tstart
c
      call init_val(time,dc13_coeff,h_coeff,C_old,C13_old,h_old,
     $     C_sim,error)
      if (error) then
         write(*,*) '--RUN_SIMULATION: an error occured in init_val'
         goto 9900
      endif
c
c ---------------------------------------------------------------------
c     loop over the simulation
c ---------------------------------------------------------------------
c 
      do k=1,nstep
c
c ---------------------------------------------------------------------
c     print out time every year
c ---------------------------------------------------------------------
c
         if (mod(k,400) .eq. 1) then
            write(*,*) '--RUN_SIMULATION: k,time : ',k,time
         endif
c
c ---------------------------------------------------------------------
c     get input values
c ---------------------------------------------------------------------
c
         call input(time,ts,dc13_coeff,h_coeff,ws_coeff,
     $        comp_fco2_online,
     $        c_new,c13_new,h_new,
     $        c,c13,h,
     $        temp,sal,ws,
     $        salk,fco2_o,pco2_a,dc13_a,dcdz_var,kz_var,error)
         if (error) then
            write(*,*) '--RUN_SIMULATION: an error occured in input.'
            goto 9900
         endif
c
c -----------------------------------------------------------------------
c     for debugging purposes
c -----------------------------------------------------------------------
c
#ifdef debugging
C         if (mod(k,90) .eq. 0) then
            write(*,*) '--RUN_SIMULATiON: time = ',time
            write(*,*) '   k             : ',k 
            write(*,*) '   temp,sal      : ',temp,sal
            write(*,*) '   dic,c13,dc13  : ',c_new,c13_new,
     $           c13_to_delta(c13_new,c_new)
            write(*,*) '   fco2_o,pco2_a : ',fco2_o,pco2_a
            write(*,*) '   dc13atm,ws    : ',dc13_a,ws
            write(*,*) '   kz,h          : ',kz_var,h
            write(*,*) '   hnew,hold     : ',h_new,h_old
C         endif
#endif
c     
c
c -----------------------------------------------------------------------
c     end of for debugging purposes
c -----------------------------------------------------------------------
c
c
c ---------------------------------------------------------------------
c     assign observed derivatives of c to arrays 
c ---------------------------------------------------------------------
c
         D_c(7,k) = c_new - c_old 
         fluxes(7,k) = rates_to_fluxes(D_c(7,k),ts,h)
c
c --------------------------------------------------------------------
c     calculate the observed rate of change of C13
c       this is !! NOT !! equal to C13_new - C13_old
c --------------------------------------------------------------------
c
         dC13_old = c13_to_delta(c13_old,c_old)
         dc13_new = c13_to_delta(c13_new,c_new)
         dc13 = c13_to_delta(c13,c)
c
         D_C13(7,k) = r_std / (1.d0 + r_std) * 
     $        (c * 1.d-3 * (dc13_new - dc13_old) +
     $        1.d-3 * dc13 * (c_new - c_old) +
     $        (c_new - c_old))
c     
c ---------------------------------------------------------------------
c     calculate contribution of air-sea exchange
c ---------------------------------------------------------------------
c
         call gasexchange(k,ts,
     $        pv_relship,D_fco2_corr,gasex_fact,
     $        c,c13,h,temp,sal,fco2_o,pco2_a,dc13_a,ws,
     $        sum_kex,fco2_a_obs,pco2_a_dry,
     $        D_c,D_c13,fluxes)
c
c ---------------------------------------------------------------------
c     calculate contribution of vertical diffusion
c ---------------------------------------------------------------------
c
         call vert_diff(k,ts,
     $        ent_scheme,
     $        h_th,c_th,dc13_th,
     $        const_kz,kz_const,kz_var,
     $        const_vertgrad,dcdz_const,dcdz_var,D_dcdz,
     $        ddc13_dc_bml,diff_fact,
     $        c,c13,h,
     $        D_c,D_c13,fluxes)

c ---------------------------------------------------------------------
c     calculate contribution of entrainment
c ---------------------------------------------------------------------
c
         if (ent_scheme .eq. 'in') then
c
            call intvertgrad_ent(k,time,ts,
     $           const_vertgrad,dcdz_const,dcdz_var,D_dcdz,
     $           ddc13_dc_bml,ent_fact,
     $           h_new,h_old,
     $           c,c13,h,
     $           D_c,D_c13,fluxes)
c
         elseif (ent_scheme .eq. 'mc') then
c
c            write(*,*) 'mc entrainment'
            call masscons_ent(k,ts,
     $           h_th,c_th,dc13_th,ent_fact,
     $           h_new,h_old,
     $           c,c13,h,
     $           D_c,D_c13,fluxes)
c
         elseif (ent_scheme .eq. 'ep') then
c
            call episevent_ent(k,time,ts,
     $           dc13_coeff,h_coeff,dt_ent_re,
     $           const_vertgrad,dcdz_const,dcdz_var,D_dcdz,
     $           ddc13_dc_bml,ent_fact,
     $           h_new,h_old,
     $           c,c13,h,
     $           D_c,D_c13,fluxes,error)
c
         elseif (ent_scheme .eq. 'ma') then
c
            call matear_ent(k,time,ts,
     $           const_vertgrad,dCdz_const,dCdz_var,D_dcdz,
     $           ddC13_dC_bml,const_lent,lent_const,dlent_dh,ent_fact,
     $           h_new,h_old,
     $           c,c13,h,
     $           D_C,D_c13,fluxes)
c
         endif
c
c ---------------------------------------------------------------------
c     calculate contribution of horizontal advection
c          part I
c ---------------------------------------------------------------------
c
         if (incl_adv) then
            call advection1(k,ts,sol_scheme,
     $           u,dc_dh,ddc13_dh,
     $           c,c13,h,
     $           dc13_dh,
     $           D_C,D_c13,fluxes)
         else
            D_c(4,k) = 0.0d0
            D_c13(4,k) = 0.0d0
            fluxes(4,k) = 0.0d0
         endif
c
c ---------------------------------------------------------------------
c     calculate contribution of Net community production
c ---------------------------------------------------------------------
c
         call net_comm_prod(k,ts,sol_scheme,
     $        D_dc13_org,
     $        c13_old,c13_new,c_old,c_new,
     $        c,c13,h,temp,sal,fco2_o,dc13_dh,dc_dh,
     $        D_C,D_c13,fluxes,dc13_org,sum_fract)
c
c ---------------------------------------------------------------------
c     calculate contribution of horizontal advection
c          part II (only for sol_scheme 'adv')
c ---------------------------------------------------------------------
c
         if (sol_scheme .eq. 'adv') then
            call advection2(k,ts,sol_scheme,
     $           dc_dh,dc13_dh,
     $           c,c13,h,
     $           u_var,D_C,D_c13,fluxes)
         endif
c     
c ---------------------------------------------------------------------
c     assign observed values to arrays
c       fco2_a_obs is already done in subroutine gasexchange
c       for fco2_o_obs it is necessary to subtract correction term
c ---------------------------------------------------------------------
c
         c_obs(k) = c_old
         dc13_obs(k) = c13_to_delta(c13_old,c_old)
         h_obs(k) = h_old
c
         temp_obs(k) = temp
         sal_obs(k) = sal
         ws_obs(k) = ws
         fco2_o_obs(k) = fco2_o - D_fco2_corr
         pco2_a_obs(k) = pco2_a
c         write(*,*) 'pco2_a,fco2_a : ',k,pco2_a_obs(k),fco2_a_obs(k)
         dc13_a_obs(k) = dc13_a
         salk_obs(k) = salk
         dCdz_obs(k) = dcdz_var
         kz_obs(k) = kz_var
c
c ---------------------------------------------------------------------
c     calculate the remaining components, depending on solution scheme
c        sol. scheme 'adv'
c ---------------------------------------------------------------------
c
         if (sol_scheme .eq. 'adv') then
            sum = 0.d0
            sum_c13 = 0.d0
            do i = 1,5
               sum = sum + D_C(i,k)
               sum_c13 = sum_c13 + D_C13(i,k)
            enddo
            D_C(6,k) = sum
            D_c13(6,k) = sum_c13
            fluxes(6,k) = rates_to_fluxes(D_C(6,k),ts,h)
            c_sim(k+1) = C_new
c
c ---------------------------------------------------------------------
c     standard solution scheme
c ---------------------------------------------------------------------
c
         elseif (sol_scheme .eq. 'std') then
            c_sim(k+1) = c_sim(k) + D_c(6,k)
         endif         
c
c ---------------------------------------------------------------------
c     calculate remaining fluxes : depth component, total calc and
c         total observed
c ---------------------------------------------------------------------
c
         fluxes(8,k) = c
     $        * (h_new-h_old)/(ts*sec_day*day_year)*rho*1.d-6 
         fluxes(9,k) = fluxes(6,k) + fluxes(8,k) ! total calc
         fluxes(10,k) = fluxes(7,k) + fluxes(8,k) ! total obs


c ---------------------------------------------------------------------
c     calculate dC13 values of the five processes
c ---------------------------------------------------------------------
CJML Fix floating point exceptions
CJML         do i = 1,4
CJML            if ((D_C(i,k) .ne. 0.d0) .and. 
CJML     $           (D_C13(i,k)/D_C(i,k) .ne. 1.d0)) then
CJML               nsum_dc13_prc(i) = nsum_dc13_prc(i) + 1
CJML               r_prc = (D_C13(i,k)/D_C(i,k)) / 
CJML     $              (1-(D_C13(i,k)/D_C(i,k)))
CJML               dc13_tmp = (r_prc/r_std - 1.d0) * 1.d3
CJML               sum_dc13_prc(i) = sum_dc13_prc(i) + dc13_tmp
CJML            endif
CJML         enddo
c
         do i = 1,4
            denom = D_C(i,k)
            if (abs(denom) > 10.d0*tiny(denom)) then
               ratio = D_C13(i,k) / denom
               if (abs(ratio-1.d0) > 10.d0*epsilon(1.d0)) then
                  nsum_dc13_prc(i) = nsum_dc13_prc(i) + 1
                  denom2 = 1.d0 - ratio
                  denom2 = sign(max(abs(denom2), 10.d0*tiny(1.d0)),
     $                           denom2)
                  r_prc = ratio / denom2
                  dc13_tmp = (r_prc/r_std - 1.d0) * 1.d3
                  sum_dc13_prc(i) = sum_dc13_prc(i) + dc13_tmp
               endif
            endif
         enddo
c
c ---------------------------------------------------------------------
c     for debugging purposes
c ---------------------------------------------------------------------
c
#ifdef debugging

C         if (mod(k,90) .eq. 0) then
            write(*,9000) time,k
 9000       format(/'--RUN_SIMULATION: for debugging '/
     $              '      time : ',f6.1/
     $              '      k    : ',i3/)
            write(*,9010) (D_C(i,k),i=1,4)
 9010       format('   D_C_ex, D_C_diff, D_C_ent, D_C_adv : ',
     $           4(f8.6,1x))
            sum = 0.d0
            do i = 1,4
               sum = sum + D_C(i,k)
            enddo
            write(*,9020) sum
 9020       format('   sum of these 4 processes           : ',f8.6)
c
            write(*,9100) (D_C13(i,k),i=1,4)
 9100       format('   D_C13_ex, D_C13_diff, D_C13_ent, D_C13_adv : ',
     $           4(f10.8,1x))
            sum_c13 = 0.d0
            do i = 1,4
               sum_c13 = sum_c13 + D_C13(i,k)
            enddo
            write(*,9120) sum_c13
 9120       format('   sum of these 4 processes           : ',f10.8)
c
            write(*,9140) D_C13(7,k)
 9140       format('   observed rate of change of C13     : ',f10.8)
c
            dc13_old = c13_to_delta(c13_old,c_old)
            dc13_new = c13_to_delta(c13_new,c_new)
            dc13 = c13_to_delta(c13,c)
            write(*,9160) C * 1.d-3*(dc13_new - dc13_old)
 9160       format('   C*ddc13/dt                         : ',f10.8)
            write(*,9180) 1d-3*(dc13)*D_C(7,k)
 9180       format('   dc13*dc/dt                         : ',f10.8)
            write(*,9185) D_C(7,k)
 9185       format('   dc/dt                              : ',f10.8)
c
            write(*,9260) D_C(6,k)
 9260       format('   calculated rate of change          : ',f8.6)
c
            write(*,9280) D_C(7,k)
 9280       format('   observed rate of change            : ',f8.6)
c
            write(*,*)
cc
c         write(*,*)
c         D_c13_calc = 0.d0
c         do i = 1,5
c            D_c13_calc = D_c13_calc + D_c13(i,k)
c         enddo
c         write(*,*) 'D_c13_calc, D_c13_obs : ',D_c13_calc,
c     $        D_c13(6,k)
cc     
c         write(*,*)
c         do i = 1,5
c            if (D_c(i,k) .ne. 0.d0) then
c               rstar = D_c13(i,k) / D_c(i,k)
c               dc13_i(i) = 1.d3 * (1.d0/r_std*(rstar/(1.d0-rstar)) 
c     $              - 1.d0)
c            else
c               dc13_i(i) = 0.d0
c            endif
c         enddo
c         write(*,*) 'dc13(ex),dc13(diff),dc13(ent),dc13(adv),',
c     $        'dc13(ncp) : ',(dc13_i(i),i=1,5)
cc
c         write(*,*)
c         dc13_old = c13_to_delta(c13_old,c_old)
c         dc13_new = c13_to_delta(c13_new,c_new)
cc
c         D_dc13_calc = 0.d0
c         do i = 1,5
c            D_dc13(i) = (D_C(i,k) * dc13_i(i) + c_old * dc13_old)/
c     $           (c_old + D_C(i,k)) - dc13_old
c            D_dc13_calc = d_dc13_calc + D_dc13(i)
c         enddo
c         write(*,*) 'D_dc13(ex),D_dc13(diff),D_dc13(ent),D_dc13(adv),',
c     $        'D_dc13(ncp) : ',(D_dc13(i),i=1,5)
c         write(*,*) 'D_dc13(calc),D_dc13(obs) : ',D_dc13_calc,
c     $        dc13_new - dc13_old
cc
cc
c         write(*,*) 
c         write(*,*) 'D_c(ex), D_c(diff),D_c(ent),D_c(adv),D_C(ncp) : ',
c     $        (D_c(i,k),i=1,5)
cc     
c         write(*,*)
c         write(*,*) 'c_new, c_old : ',c_new,c_old
c         write(*,*) 'c13_new,c13_old : ',c13_new,c13_old
c         write(*,*) 'dc13_new,dc13_old : ',dc13_new,dc13_old
c         write(*,*) 'D_C_obs : ',c_new - c_old
c         write(*,*) 'D_C_sim : ',D_C(5,k)
c         write(*,*) 'D_C13_obs : ',c13_new - c13_old
c         write(*,*) 'D_dc13_obs : ',dc13_new - dc13_old
c         write(*,*)
C         endif
#endif
c
c ---------------------------------------------------------------------
c     set new values to old values
c ---------------------------------------------------------------------
c
         c_old = c_new
         c13_old = c13_new
         h_old = h_new
c
c ---------------------------------------------------------------------
c     update time
c ---------------------------------------------------------------------
c
         time = time + ts
c
c ---------------------------------------------------------------------
c     end of loop over k
c ---------------------------------------------------------------------
c
      enddo
c
c ---------------------------------------------------------------------
c     calculate integrated values
c ---------------------------------------------------------------------
c
      call calc_int_val(D_c,D_c13,fluxes,tstart,tend,ts,nstep,
     $     int_fluxes,int_rates,nyears)
c
c ---------------------------------------------------------------------
c     calculate run statistic and averages
c ---------------------------------------------------------------------
c
      call calc_run_stat(D_C,c_sim,c_obs,tstart,tend,ts,nstep,
     $     sol_scheme,
     $     sum_dc13_prc,nsum_dc13_prc,sum_fract,sum_kex,
     $     dc13_prc,avg_fract,avg_kex,
     $     chisq,rsq,closure,amplitude,amplitude_obs,
     $     c_min,c_max,c_min_obs,c_max_obs,
     $     costfn)
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c     
      return
c
c --------------------------------------------------------------------
c     error handling
c --------------------------------------------------------------------
c        
 9900 write(*,9920) 
 9920 format(/'--RUN_SIMULATION: an error occured.')
      error = .true.
      return
      end

