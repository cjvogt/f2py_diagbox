c =====================================================================
c
c     DOUBLE PRECISION FUNCTION FUNC
c
c     PURPOSE: this is the routine which calculates the value of
c                the cost function which is to minimize.
c              communication with the main program (for all the
c              the run parameters except those given through p
c              is done through common blocks, defined in
c              conj_grad_common.h
c     
c     VARIABLES:
c 
c     IN:   p            : input array of parameters
c     OUT:  func         : value of function at x
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     13.08.96   ng     first implementation
c 
c =====================================================================
c
      double precision function func(p)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "conj_grad.h"
# include "sg_params.h"
# include "simulation.h"
# include "gasex_params.h"
# include "diffent_params.h"
# include "adv_params.h"
# include "ncp_params.h"
# include "conj_grad_common.h"
c
c --------------------------------------------------------------------
c     variables to communicate with run_simulation
c --------------------------------------------------------------------
c 
#include "diagboxmod.h"
#include "processes.h"
#include "periods.h"
#include "intval.h"
#include "slabs.h"
#include "obs_values.h"
#include "statistics.h"
#include "averages.h"
      logical error
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c 
      integer np
c
c --------------------------------------------------------------------
c     external functions
c --------------------------------------------------------------------
c 
      double precision cost_function
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     info
c --------------------------------------------------------------------
c 
c      write(*,*) '--FUNC: calling run_simulation '
c      write(*,*) '    parameter values : '
c      do np = 1,nparams
c         write(*,*) '    np,p(np) : ',np,p(np)
c      enddo
c      write(*,*) 'incl_adv : ',incl_adv
c      write(*,*) 'const_lent : ',const_lent
c
c --------------------------------------------------------------------
c     assign conjugate gradient parameter values to regular parameters
c --------------------------------------------------------------------
c 
      assign = .false.
      call assign_params(assign,
     $     gasex_fact, sg_gasex_fact,
     $     D_fco2_corr,sg_D_fco2_corr,
     $     diff_fact,sg_diff_fact,D_dcdz,sg_D_dcdz,
     $     ddC13_dC_bml,sg_ddC13_dC_bml,
     $     const_lent,lent_const,sg_lent_const,
     $     dlent_dh,sg_dlent_dh,   
     $     D_dc13_org,sg_D_dc13_org, 
     $     u,sg_u,dc_dh,sg_dc_dh,ddC13_dh,sg_ddC13_dh,   
     $     p,sg_p,p_init) 
c
c --------------------------------------------------------------------
c     run the simulation
c --------------------------------------------------------------------
c 
c      write(*,*) '--FUNC: calling simulation...'
c
         call run_simulation(
     $     tstart,tend,ts,nstep,
     $     dc13_coeff,
     $     pv_relship,D_fco2_corr,gasex_fact,ws_coeff,      
     $     ent_scheme,h_coeff,
     $     h_th,c_th,dc13_th,
     $     dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $     const_Kz,Kz_const,
     $     const_vertgrad,dCdz_const,D_dcdz,
     $     ddC13_dC_bml,diff_fact,
     $     incl_adv,dC_dh,ddC13_dh,u,
     $     D_dc13_org,
c
     $     D_c,D_c13,fluxes,
     $     c_obs,dc13_obs,h_obs,
     $     temp_obs,sal_obs,
     $     fco2_o_obs,pco2_a_obs,fco2_a_obs,
     $     dc13_a_obs,ws_obs,dcdz_obs,kz_obs,
     $     c_sim,
     $     int_rates,int_fluxes,
     $     chisq,rsq,closure,amplitude,amplitude_obs,
     $     c_min,c_max,c_min_obs,c_max_obs,costfn,
     $     avg_fract,dc13_prc,avg_kex,
     $     error)
c
c --------------------------------------------------------------------
c     calculate cost function : part of parameter deviation
c -------------------------------------------------------------------- 
c
         par_sum = 0.d0
         do np = 1,nparams
            par_sum = par_sum + p(np)**2
c            par_sum = par_sum + ((p(np) - p_init(np))/sg_p(np))**2
         enddo
c
c --------------------------------------------------------------------
c     calculate total cost function
c --------------------------------------------------------------------
c 
         func = cost_function(chisq,closure,par_sum,nparams,
     $        amplitude,amplitude_obs,
     $        c_min,c_max,c_min_obs,c_max_obs,
     $        chisq_term,clos_term,par_term,ampl_term)
c
c
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c     
      return

      end











