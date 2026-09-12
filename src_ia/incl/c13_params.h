cc=================== include file c13_params.h =============================
c
c parameters for C13 and dC13
c
c REVISIONS:
c
c date     author    remarks
c
c 2.08.96  ng        first implementation
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c
c r_Std          dp  s  C13/C12 ratio of PDB standard
c rstar_Std      dp  s  C13/(C12+C13) ratio of PDB standard
c
c---------------------------------------------------------------------------
c
      double precision r_std,rstar_std
      parameter (r_std = 0.0112372d0, rstar_std = 0.011112328d0)



