c===================== include file sg_params.h =========================
c
c variables defining special uncertainties
c
c REVISIONS:
c
c date     author    remarks
c
c 13.8.96  ng        first implementation
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c
c sg_c          dp  s   1-sigma uncertainty of C [umol kg-1]
c sg_dcdt       sp  s   1-sigma uncertainty of dcdt [umol kg-1 day-1]
c                          (estimated from figure)
c
c===========================================================================
c
      double precision sg_c,sg_dcdt
      parameter(sg_c = 2.d0, sg_dcdt = 0.02d0)
      
