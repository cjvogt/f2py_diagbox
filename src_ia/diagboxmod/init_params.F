c =====================================================================
c
c     SUBROUTINE INIT_PARAMS
c 
c     PURPOSE: this routine initializes the values of the various
c                parameters
c
c     VARIABLES:
c 
c     IN : 
c
c     OUT: 
c 
c     REVISION:
c
c     date     author  remarks
c   
c     29.08.01   ng    first implementation, used to be an include
c                        file with data statements
c
c =====================================================================
c
      subroutine init_params(incl_adv,plot_out,
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
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "diagboxmod.h"
# include "simulation.h"
# include "adv_params.h"
# include "gasex_params.h"
# include "ncp_params.h"
# include "diffent_params.h"
c     
c --------------------------------------------------------------------
c     statement
c --------------------------------------------------------------------
c     
      write(*,*) '--INIT_PARAMS: initializing parameters...'
c     
c --------------------------------------------------------------------
c     general parameters
c --------------------------------------------------------------------
c     
      incl_adv = .false. 
      plot_out = .false. 
c     
c --------------------------------------------------------------------
c     simulation control parameters
c --------------------------------------------------------------------
c     
      tstart = 1983.0d0
      tend = 2001.d0
      ts = 0.0025d0 
c     
c --------------------------------------------------------------------
c     air-sea gasexchange parameters
c --------------------------------------------------------------------
c     
      pv_relship = 'wa'
c
      D_fco2_corr = 0.0d0 
c      sg_D_fco2_corr = 10.d0 
      sg_D_fco2_corr = 1.d0 
c
      gasex_fact = 1.0d0 
      sg_gasex_fact = 0.1d0
c
      ws_coeff = 'ih'
c     
c --------------------------------------------------------------------
c     diffusion - entrainment parameters
c --------------------------------------------------------------------
c     h_coeff = 'a'

      ent_scheme = 'ma'
      const_kz = .false.
c
      kz_const = 2.d-5
      sg_kz_const = 1.d-5
c
      const_vertgrad = .false.
      dCdz_const = 0.45d0
c
      D_dcdz = 0.0d0
      sg_D_dcdz = 0.1d0
c
      h_th    = 200.d0
      c_th    = 2066.31d0
      dc13_th = 1.28d0
c
      const_lent = .true.
      lent_const = 12.d0
c      sg_lent_const = 6.d0
      sg_lent_const = 3.d0
c
      dlent_dh = 0.11
      sg_dlent_dh = 0.4
c
      ddC13_dC_bml = -0.0048d0 
      sg_ddC13_dC_bml = 0.0008d0
c
      ent_fact = 1.0d0
      diff_fact = 1.0d0
      sg_diff_fact = 0.5d0
      sg_ent_fact = 0.5d0
c     
c --------------------------------------------------------------------
c     advection parameters
c --------------------------------------------------------------------
c     
      dC_dh = 1.1d-5
      sg_dC_dh = 0.3d-5
c
      ddC13_dh = -1.0d-7
      sg_ddC13_dh = 1.0d-7
c
      u = 0.05d0
      sg_u = 0.05d0
c     
c --------------------------------------------------------------------
c     net community production parameters
c --------------------------------------------------------------------
c     
      D_dc13_org = 0.0d0
      sg_D_dc13_org = 3.0d0
c     
c --------------------------------------------------------------------
c     return to calling routine
c --------------------------------------------------------------------
c     
      return
      end
