c =====================================================================
c
c     SUBROUTINE ASSIGN_PARAMS
c 
c     PURPOSE: assigns and deassigns parameters for optimization,
c                monte carlo analysis and sensitivity tests.
c              action is determined by the value of the logical
c              variable assign.
c
c     VARIABLES:
c
c     IN : assign        : logical variable specifying action
c                             <TRUE>  assign variables to paramets p
c                             <FALSE> assign parameters to variables
c                        
c     IN/OUT: 
c          gasex_fact      : value of gasexchange factor
c          sg_gasex_fact   : 1-sg uncertainty of gasexchange factor
c          D_fco2_corr     : value of correction term for fCO2 ocean
c          sg_D_fco2_corr  : 1-sg uncertainty of correction term for fCO2 ocean
c          diff_fact       : value of diffusion factor
c          sg_diff_fact    : 1-sg uncertainty of diffusion factor
c          D_dcdz          : value of correction term for dcdz
c          sg_D_dcdz       : 1-sg uncertainty of correction term for dcdz
c          ddC13_dC_bml    : value of ratio of vertical dC13 and C gradients
c          sg_ddC13_dC_bml : 1-sg uncertainty of ratio of vertical dC13 and 
c                                C gradients
c          const_lent      : logical specifyin if lent is constant or not
c          lent_const      : value of constant entrainment length
c          sg_lent_const   : 1-sg uncertainty of constant entrainment length
c          dlent_dh        : value of slope of entrainment length with h
c          sg_dlent_dh     : 1-sg uncertainty of slope of entrainment 
c                                length with h
c          D_dc13_org      : value of correction term for dC13org
c          sg_D_dc13_org   : 1-sg uncertainty  of correction term for dC13org
c          u               : value of horizontal adv velocity
c          sg_u            : 1-sg uncertainty of horizontal adv velocity
c          dc_dh           : value of horizontal C gradient
c          sg_dc_dh        : 1-sg uncertainty of horizontal C gradient
c          ddC13_dh        : value of horizontal dC13 gradient
c          sg_ddC13_dh     : 1-sg uncertainty of horizontal dC13 gradient
c         
c          p               : value of parameters
c          sg_p            : 1-sg uncertainty of parameters
c          p_init          : initial value of parameters
c          
c
c     REVISION:
c
c     date     author  remarks
c   
c     11.09.96   ng    first implementation
c
c =====================================================================
c
      subroutine assign_params(
     $     assign,
     $     gasex_fact, sg_gasex_fact,D_fco2_corr,sg_D_fco2_corr,
     $     diff_fact,sg_diff_fact,D_dcdz,sg_D_dcdz,
     $     ddC13_dC_bml,sg_ddC13_dC_bml,
     $     const_lent,lent_const,sg_lent_const,dlent_dh,sg_dlent_dh,   
     $     D_dc13_org,sg_D_dc13_org, 
     $     u,sg_u,dc_dh,sg_dc_dh,ddC13_dh,sg_ddC13_dh,   
     $     p,sg_p,p_init)          
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "diagboxmod.h"
# include "gasex_params.h"
# include "diffent_params.h"
# include "adv_params.h"
# include "ncp_params.h"
# include "conj_grad.h"
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
c      
c =====================================================================
c     begin of executable code
c =====================================================================
c
      if (assign) then
c
         write(*,*) '--ASSIGN_PARAMS: assigning parameters...'
c
c --------------------------------------------------------------------
c     assign variables to parameters
c       normalize them with
c       the formula p' = (p - p_init)/sigma
c --------------------------------------------------------------------
c     
#ifndef wo_optimization_ent
         p(1) = gasex_fact
         p(2) = D_fco2_corr
         p(3) = diff_fact
         p(4) = D_dcdz  
         p(5) = ddC13_dC_bml
         if (const_lent) then
            p(6) = lent_const
         else
            p(6) = dlent_dh
         endif
         p(7) = D_dc13_org
         p(8) = u
#  ifdef with_horiz_grad
         p(9) = dc_dh
         p(10) = ddC13_dh
#  endif
         sg_p(1) = sg_gasex_fact
         sg_p(2) = sg_D_fco2_corr
         sg_p(3) = sg_diff_fact
         sg_p(4) = sg_D_dcdz
         sg_p(5) = sg_ddC13_dC_bml
         if (const_lent) then
            sg_p(6) = sg_lent_const
         else
            sg_p(6) = sg_dlent_dh
         endif
         sg_p(7) = sg_D_dc13_org
         sg_p(8) = sg_u
#  ifdef with_horiz_grad
         sg_p(9) = sg_dc_dh
         sg_p(10) = sg_ddC13_dh
#  endif
c
c
#else
         p(1) = gasex_fact
         p(2) = D_fco2_corr
         p(3) = diff_fact
         p(4) = D_dcdz  
         p(5) = ddC13_dC_bml
         p(6) = D_dc13_org
         p(7) = u
#  ifdef with_horiz_grad
         p(8) = dc_dh
         p(9) = ddC13_dh
#  endif
         sg_p(1) = sg_gasex_fact
         sg_p(2) = sg_D_fco2_corr
         sg_p(3) = sg_diff_fact
         sg_p(4) = sg_D_dcdz
         sg_p(5) = sg_ddC13_dC_bml
         sg_p(6) = sg_D_dc13_org
         sg_p(7) = sg_u
#  ifdef with_horiz_grad
         sg_p(8) = sg_dc_dh
         sg_p(9) = sg_ddC13_dh
#  endif
#endif
c
c --------------------------------------------------------------------
c     assign parameters to variables
c --------------------------------------------------------------------
c     
      else
         write(*,*) '--ASSIGN_PARAMS: deassigning parameters...'
#ifndef wo_optimization_ent
         gasex_fact   =  sg_p(1)*p(1) + p_init(1)
         D_fco2_corr  =  sg_p(2)*p(2) + p_init(2)
         diff_fact    =  sg_p(3)*p(3) + p_init(3)
         D_dcdz       =  sg_p(4)*p(4) + p_init(4)
         ddC13_dC_bml =  sg_p(5)*p(5) + p_init(5)
         if (const_lent) then
            lent_const   =  sg_p(6)*p(6) + p_init(6)
         else
            dlent_dh     =  sg_p(6)*p(6) + p_init(6)
         endif
         D_dc13_org   =  sg_p(7)*p(7) + p_init(7)
         u            =  sg_p(8)*p(8) + p_init(8)
#  ifdef with_horiz_grad
         dC_dh        =  sg_p(9)*p(9) + p_init(9)
         ddC13_dh     =  sg_p(10)*p(10) + p_init(10)
#  endif
#else 
         gasex_fact   =  sg_p(1)*p(1) + p_init(1)
         D_fco2_corr  =  sg_p(2)*p(2) + p_init(2)
         diff_fact    =  sg_p(3)*p(3) + p_init(3)
         D_dcdz       =  sg_p(4)*p(4) + p_init(4)
         ddC13_dC_bml =  sg_p(5)*p(5) + p_init(5)
         D_dc13_org   =  sg_p(6)*p(6) + p_init(6)
         u            =  sg_p(7)*p(7) + p_init(7)
#  ifdef with_horiz_grad
         dC_dh        =  sg_p(8)*p(8) + p_init(8)
         ddC13_dh     =  sg_p(9)*p(9) + p_init(9)
#  endif
#endif
      endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

