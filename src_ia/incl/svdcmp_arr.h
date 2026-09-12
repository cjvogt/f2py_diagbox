c===================== include file svdcmp_arr.h =========================
c
c variables for singular value decomposition
c
c REVISIONS:
c
c date     author    remarks
c
c 15.8.96  ng        first implementation
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c        
c w            dp  a   diagonal vector returned from singular value 
c                          decomposition routine svdcmp 
c v            dp  a   matrix v returned from svdcmp
c
c===========================================================================
c
      double precision w(nparams), v(nparams)
c
      
