c =====================================================================
c
c     SUBROUTINE NET_COMM_PROD
c 
c     PURPOSE: this routine calculates the contribution of the
c                  net community production. This contribution
c                  is calculated from the difference between
c                  the observed rate of change of C13 and the
c                  sum of the rate of changes from the other
c                  processes (air-sea exchange, vertical diffusion,
c                             entrainment and advection)
c              the calculation scheme has been drastically changed from
c                  the original version.
c                  
c
c     VARIABLES:
c 
c     IN : k              : time loop counter
c          ts             : time step
c          D_dc13_org     : correction term for dC13 of organic matter
c          C13_old        : obs C13 in the mixed layer at timestep t
c          C13_new        : obs C13 in the mixed layer at timestep t+ts
c          c13            : obs C13 in the mixed layer at timestep t+0.5*ts
c          c_old          : obs C in the mixed layer at timestep t
c          c_new          : obs C in the mixed layer at timestep t+ts
c          c              : obs C in the mixed layer at timestep t+0.5*ts
c          h              : obs mixed layer depth
c          temp           : obs mixed layer temp
c          sal            : obs mixed layer sal
c          fco2_o         : obs mixed layer fco2
c
c     OUT: D_C            : modified array holding the rates of change for C
c          D_C13          : modified array holding the rates of change for C13
c          fluxes         : modified array holding the fluxes
c          dc13_org       : dc13 value of organic matter
c          sum_fract      : the sum of the fractionation factor for avg
c 
c     REVISION:
c
c     date     author  remarks
c   
c      2.08.96   ng    first implementation
c      5.09.96   ng    included possibility to do backwards integration
c      9.09.96   ng    calculate ncp using the more exact equation given
c                         in the paper
c     17.10.99   ng    included sol_scheme, now solving for u, rather
c                          than for dC/dt, adapting for interannual var.
c                          runs
c
c =====================================================================
c
      subroutine net_comm_prod(k,ts,sol_scheme,
     $     D_dc13_org,
     $     c13_old,c13_new,c_old,c_new,
     $     c,c13,h,temp,sal,fco2_o,dc13_dh,dc_dh,
     $     D_C,D_c13,fluxes,dc13_org,sum_fract)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
# include "input.h"
# include "ncp_params.h"
# include "c13_params.h"
# include "adv_params.h"
# include "averages.h"
      integer k
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision rates_to_fluxes,solco2kg ! functions
      double precision c13_to_delta,delta_to_ratio ! function
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision sol
      integer i
      double precision dc13_old,dc13_new,dc13
      double precision fract_tmp
      double precision term1,term2,term3,sum_c,sum_c13
      double precision r_new,r_old,r_org,r_oc
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     calculate the solubility of CO2 using Weiss, 1974
c       mol C kg-1 atm-1, this is directly numerically
c       equivalent ot umol C kg-1 uatm-1
c --------------------------------------------------------------------
c        
      sol = solco2kg(temp,sal)  ! units [mol C kg-1 atm-1 = umol kg-1 uatm-1]
c
c --------------------------------------------------------------------
c     calculate the dC13 of organic matter following
c       the relationship given by Rau et al, 1989
c     include also directly the correction
c --------------------------------------------------------------------
c   
      dc13_org = -0.8d0 * (sol * fCO2_o) - 12.6d0 ! units [per mil]
c
      dc13_org = 1.d-3 * (dc13_org + D_dc13_org) ! no units
      r_org = delta_to_ratio(dc13_org*1.d3)
c
c      write(*,*) k,dc13_org
c
c --------------------------------------------------------------------
c     calculate dc13 at time t+0.5*ts, t and t+ts
c --------------------------------------------------------------------   
c
      dc13 = c13_to_delta(c13,c)
      dC13_old = c13_to_delta(c13_old,c_old)
      dc13_new = c13_to_delta(c13_new,c_new)
c
c --------------------------------------------------------------------
c     calculate average fractionation factor
c --------------------------------------------------------------------   
c
      fract_tmp = 1.d3 * dc13_org - dc13
      sum_fract = sum_fract + fract_tmp
c
c =====================================================================
c     this is the part as calculated in the old model
c     standard solution scheme
c =====================================================================
c
      if (sol_scheme .eq. 'std') then
c
c ---------------------------------------------------------------------
c     calculate the sums of all processes except net community production
c ---------------------------------------------------------------------
c
        sum_c = 0.d0 
        sum_c13 = 0.d0
        do i = 1,4
           sum_c = sum_c + D_C(i,k)
           sum_c13 = sum_c13 + D_C13(i,k)
        enddo
c
c ---------------------------------------------------------------------
c     calculate the three different terms; accurate version
c ---------------------------------------------------------------------
c
        r_oc = delta_to_ratio(dc13)
        r_new = delta_to_ratio(dc13_new)
        r_old = delta_to_ratio(dc13_old)
c     
        term1 = (r_org*(1.d0+r_oc)/(r_org-r_oc)) * sum_c
        term2 = (1.d0+r_org)/((1.d0+r_oc)*(r_org-r_oc)) * c *
     $        (r_new - r_old)
        term3 = - ((1.d0+r_oc)*(1+r_org))/(r_org-r_oc) * sum_c13


        D_C(6,k) = term1 + term2 + term3
        fluxes(6,k) = rates_to_fluxes(D_C(6,k),ts,h)
c
c ---------------------------------------------------------------------
c     calculate the rate of change of net community production
c      and the flux
c ---------------------------------------------------------------------
c
        D_C(5,k) = D_C(6,k) - sum_c
        fluxes(5,k) = rates_to_fluxes(D_C(5,k),ts,h)

c =====================================================================
c     this is solution scheme for adv
c =====================================================================
c
      elseif (sol_scheme .eq. 'adv') then         
c      
c --------------------------------------------------------------------
c     sum up contributions by gasex, diff and entrainment
c --------------------------------------------------------------------
c         
         sum_c = 0.d0
         sum_c13 = 0.d0
          do i = 1,3
            sum_c = sum_c + D_C(i,k)
            sum_c13 = sum_c13 + D_C13(i,k)
         enddo
c      
c --------------------------------------------------------------------
c     solve for net community production
c --------------------------------------------------------------------
c          
         D_C(5,k) = (dc13_dh*(D_C(7,k)-sum_c) 
     $         - dc_dh*(D_c13(7,k)-sum_c13)) /
     $         (dc13_dh - dc_dh*(r_org/(1+r_org)))
         
         fluxes(5,k) = rates_to_fluxes(D_C(5,k),ts,h)
c
         D_C13(5,k) = r_org / (1.d0 + r_org) * D_C(5,k)
      endif
c      
c --------------------------------------------------------------------
c     for debugging purposes
c --------------------------------------------------------------------
c
#ifdef debugging
c
      if (mod(k,90) .eq. 0) then
         write(*,*)
         write(*,*) '--NET_COMM_PROD: '
         write(*,*) '  k : ',k
         write(*,*) '  dc13_org                  : ',dc13_org
         write(*,*) '  r_org                     : ',r_org
         write(*,*) '  solubility                : ',sol
         write(*,*) '  fCO2_o                    : ',fco2_o
         if (sol_scheme .eq. 'std') then
            write(*,*) '  term1                     : ',term1
            write(*,*) '  term2                     : ',term2
            write(*,*) '  term3                     : ',term3
            write(*,*) '  dc_dt_sim                 : ',D_C(6,k)
         elseif (sol_scheme .eq. 'adv') then
            write(*,*) '  D_C(7,k)               : ',D_C(7,k)
            write(*,*) '  sum_c                  : ',sum_c
            write(*,*) '  dc_dh                  : ',dc_dh
            write(*,*) '  D_C13(7,k)             : ',d_c13(7,k)
            write(*,*) '  sum_c13                : ',sum_c13
            write(*,*) '  dc13_dh                : ',dc13_dh
            write(*,*) '  denominator            : ',
     $           dc13_dh - dc_dh*(r_org/(1+r_org))
            write(*,*) '  denominator second part: ',
     $           dc_dh*(r_org/(1+r_org))
            write(*,*) '  enuminator             : ',
     $           (dc13_dh*(D_C(7,k)-sum_c) 
     $         - dc_dh*(D_c13(7,k)-sum_c13))
            write(*,*) '  enuminator first part  : ',
     $           dc13_dh*(D_C(7,k)-sum_c)
            write(*,*) '  enuminator second part  : ',
     $           dc_dh*(D_c13(7,k)-sum_c13)

         endif
            
         write(*,*) '  dc_dt_ncp                 : ',D_C(5,k)
         write(*,*)
         write(*,*) '  D_C13_ncp                 : ',D_C13(5,k)
      endif
#endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

