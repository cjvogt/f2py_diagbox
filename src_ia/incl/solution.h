c================== include file solution.h =========================
c
c variables for sensitivity analysis
c
c REVISIONS:
c
c date     author    remarks
c
c 21.8.96  ng        first implementation
c 18.10.99 ng        adapted for interannual variability runs
c
c variables:
c
c name            type   description
c---------------------------------------------------------------------------
c 
c u_var           dp a   inverse calculated horizontal velocity m s-1
c 
c===========================================================================
c
      double precision u_var(nstepmax)
