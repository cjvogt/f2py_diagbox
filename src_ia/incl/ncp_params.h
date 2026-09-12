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
