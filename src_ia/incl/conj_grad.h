c===================== include file conj_grad.h =========================
c
c variables for finding minimum using conjugate gradient technique
c
c REVISIONS:
c
c date     author    remarks
c
c 12.8.96  ng        first implementation
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c     
c assign         l  s   logical for subroutine assign_parameters
c                        true : assign values to parameters
c                        false: assign parameters to values   
c p             dp  a   array of parameters
c p_opt         dp  a   array of optimum parameters
c sg_p          dp  a   sigma of parameters
c p_init        dp  a   array of initial parameters
c nparams        i  s   dimension of parameters
c npmax          i  s   maximum dimension of p
c iter           i  s   number of iterations
c fret          dp  s   value of function at minimum
c hessian       dp  a   hessian matrix
c covm          dp  a   inverse of hessian matrix, covariance matrix
c
c sg_c          dp  s   1-sigma uncertainty of C [umol kg-1]
c sg_dcdt       sp  s   1-sigma uncertainty of dcdt [umol kg-1 day-1]
c                          (estimated from figure)
c
c===========================================================================
c
      logical assign
      double precision ftol,fret
      integer nruns
      integer nparams,iter
#ifndef wo_optimization_ent
#  ifdef with_horiz_grad
      parameter(nparams = 10)
#  else
      parameter(nparams = 8)
#  endif
#else
#  ifdef with_horiz_grad
      parameter(nparams = 9)
#  else
      parameter(nparams = 7)
#  endif
#endif
      double precision p(nparams),p_opt(nparams),p_init(nparams)
      double precision sg_p(nparams)
      double precision hessian(nparams,nparams),covm(nparams,nparams)
c
      
