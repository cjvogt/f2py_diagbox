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
