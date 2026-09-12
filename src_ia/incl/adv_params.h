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
