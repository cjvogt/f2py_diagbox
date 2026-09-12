c======================= include file averages.h ==============================
c
c variables for calculating annual averages
c
c REVISIONS:
c
c date     author    remarks
c
c 31.07.96  ng        first implementation
c 10.09.96  ng        included dC13 values for different processes
c
c variables:
c
c name           type   description
c-----------------------------------------------------------------------------
c 
c sum_fract     dp  s   sum for calculating average fractionation factor
c sum_kex       dp  s   sum for calculating average kex
c sum_dCdz      dp  s   sum for calculating average vertical DIC gradient
c sum_ddC13dz   dp  s   sum for calculating average vertical dC13 gradient
c sum_dC13_prc  dp  a   sum for calculating average dC13 values of the proc
c nsum_dc13_prc  i  a   number of values used for summing sum_dc13_prc
c
c avg_fract     dp  s   average fractionation factor [per mil]
c dc13_prc      dp  a   average dc13 value of the processes [per mil]
c avg_kex       dp  s   average kex value [mol m-2 s-1 uatm-1]
c
c-------------------------------------------------------------------------
c
      double precision sum_fract,sum_kex,sum_dCdz,sum_ddC13dz
      double precision avg_fract,avg_kex
      double precision sum_dC13_prc(4)
      integer nsum_dc13_prc(4)
      double precision dc13_prc(4)



