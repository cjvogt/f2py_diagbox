c===================== include file diffent_params.h =========================
c
c parameter of vertical diffusion and entrainment
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 13.8.96  ng        changed l_* to sg_*
c 09.9.96  ng        included const_lent
c
c variables:
c
c name            type   description
c-----------------------------------------------------------------------------
c                 
c ent_scheme      c  s   character indicating entrainment scheme to be used
c                            <in> : integrated vertical gradient scheme (DEF)
c                            <mc> : mass conservation scheme
c                            <ep> : episodic event scheme
c                            <ma> : matear scheme
c                        the <in> scheme is the new scheme developed by ng
c                                      for station BATS
c                        the <mc> scheme is the scheme of Thomas Stocker
c                                      with the assumption that the 
c                                      concentrations at a specified depth
c                                      remain constant
c                        the <ep> scheme is the original scheme developed by
c                                      C.D. Keeling used in the station 'S'
c                                      study
c                        the <ma> scheme is the scheme used by Matear (1995)
c                                      in his model (uses an entrainment
c                                      length scale)
c h_coeff         c  s   character specifying which coefficients for
c                            mixed layer depth coefficients should be taken
c                            <b> : bats          <s> : station S  
c                            <c> : bats CTD      <a> : bats avg CTD 
c                 
c dt_ent_re      dp  s   recurrence period of entrainment (only important
c                            for the <ep> scheme
c const_Kz        l  s   logical indicating whether const Kz is used or
c                            seasonally varying Kz
c kz_const       dp  s   value of constant kz [m2 s-1]
c const_vertgrad  l  s   logical indicating whether const vertical gradients
c                            are used or seasonally varying
c dCdz_const     dp  s   value of constant dCdz [umol kg-1 m-1]
c D_dcdz         dp  s   correction term for dCdz [umol kg-1 m-1]
c dd13C_dC_bml   dp  s   ratio of ddC13/dz and dC/dz below mixed layer
c                           [per mil kg umol-1]
c diff_fact      dp  s   multiplication factor for vertical diffusion
c ent_fact       dp  s   multiplication factor for entrainment
c h_th           dp  s   depth at which thermocline values for mc entrainment
c                           scheme are taken
c c_th           dp  s   C concentration at h_th [umol kg-1]
c dC13_th        dp  s   dC13 concentration at h_th [per mil]
c 
c const_lent     dp  l   logical indicating use of a constant ent_length
c lent_const     dp  s   constant entrainment length scale 
c                                (for Matear entrainment) [m]
c dlent_dh       dp  s   slope of ent length vs h [no units]
c
c sg_ddC13_dC_bml dp  s   sigma uncertainty of dd13C_dC_bml (for OP and MC)
c sg_diff_fact    dp  s   sigma uncertainty of diff_fact (for OP and MC)
c sg_ent_fact             sigma uncertainty of ent_fact (for OP and MC)
c sg_kz_const     dp  s   sigma uncertainty of kz (for OP and MC)
c sg_lent_const   dp  s   sigma uncertainty of ent_length (for OP and MC)
c sg_D_dcdz       dp  s   sigma uncertainty of dcdz
c                
c=============================================================================
c
      character ent_scheme*2,h_coeff*1
      logical const_kz,const_vertgrad,const_lent
      double precision dt_ent_re
      double precision kz_const,dCdz_const,D_dcdz
      double precision ddC13_dC_bml,diff_fact,ent_fact
      double precision h_th,C_th,dC13_th
      double precision lent_const,dlent_dh
c
      double precision sg_ddC13_dC_bml
      double precision sg_ent_fact,sg_diff_fact
      double precision sg_kz_const,sg_lent_const,sg_D_dcdz
      double precision sg_dlent_dh
      
