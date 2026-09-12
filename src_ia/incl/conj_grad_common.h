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
      
