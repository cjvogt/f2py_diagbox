c===================== include file init_params.h =========================
c
c initialize values of parameters. for definition of the various
c    parameters see the individual *.h files
c
c REVISIONS:
c
c date     author    remarks
c
c 2.08.96  ng        first implementation
c 13.8.96  ng        changed l_* to sg_*
c
c========================================================================
c     
c --------------------------------------------------------------------
c     general parameters
c --------------------------------------------------------------------
c     
      data incl_adv / .false. /
      data plot_out / .false. /
c     
c --------------------------------------------------------------------
c     simulation control parameters
c --------------------------------------------------------------------
c     
      data tstart,tend /83.0d0, 99.d0/
      data ts / 0.0025d0 /
c     
c --------------------------------------------------------------------
c     air-sea gasexchange parameters
c --------------------------------------------------------------------
c     
      data pv_relship /'wa'/
c
      data D_fco2_corr / 0.0d0 /
      data sg_D_fco2_corr / 10.d0 /
c
      data gasex_fact / 1.0d0 /
      data sg_gasex_fact /0.3d0/
c
      data ws_coeff /'ih'/
c     
c --------------------------------------------------------------------
c     diffusion - entrainment parameters
c --------------------------------------------------------------------
c     data h_coeff /'a'/

      data ent_scheme /'ma'/
      data const_kz /.false./
c
      data kz_const /2.d-5/
      data sg_kz_const /1.d-5/
c
      data const_vertgrad /.false./
      data dCdz_const /0.45d0/
c
      data D_dcdz /0.0d0/
      data sg_D_dcdz /0.1d0/
c
      data h_th    /200.d0/
      data c_th    /2066.31d0/
      data dc13_th /1.28d0/
c
      data const_lent /.true./
      data lent_const /12.d0/
      data sg_lent_const /6.d0/
c
      data dlent_dh /0.11/
      data sg_dlent_dh /0.4/
c
      data ddC13_dC_bml /-0.0048d0 /
      data sg_ddC13_dC_bml /0.0008d0/
c
      data ent_fact /1.0d0/
      data diff_fact /1.0d0/
      data sg_diff_fact /0.5d0/
      data sg_ent_fact /0.5d0/
c     
c --------------------------------------------------------------------
c     advection parameters
c --------------------------------------------------------------------
c     
      data dC_dh /1.1d-5/
      data sg_dC_dh /0.3d-5/
c
      data ddC13_dh /-1.0d-7/
      data sg_ddC13_dh /1.0d-7/
c
      data u /0.05d0/
      data sg_u /0.05d0/
c     
c --------------------------------------------------------------------
c     net community production parameters
c --------------------------------------------------------------------
c     
      data D_dc13_org / 0.0d0 /
      data sg_D_dc13_org /3.0d0/

