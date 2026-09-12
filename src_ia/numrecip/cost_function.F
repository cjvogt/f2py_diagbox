c =====================================================================
c
c     DOUBLE PRECISION FUNCTION COST_FUNCTION
c
c     PURPOSE: this is the routine which calculates the value of
c                the cost function from the different terms
c     
c     VARIABLES:
c 
c     IN:   chisq         : value of chi**2
c           closure       : value of annual lack of closure
c           par_sum       : value of sum of parameter deviations
c           np            : number of parameters
c           amplitude     : value of simulated amplitude
c           amplitude_obs : value of observed amplitude
c           c_min         : value of minimum C
c           c_max         : value of maximum C
c           c_min_obs     : value of observed minimum C
c           c_max_obs     : value of observed maximum C
c     OUT:  cost_function : value of cost function
c           chisq_term    : value of chi**2 term
c           clos_term     : value of closure term
c           par_term      : value of parameter term
c           ampl_term     : value of amplitude term
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     14.08.96   ng     first implementation
c      5.09.96   ng     changed definition of cost function, introduced
c                        c_min,c_max,c_min_obs,c_max_obs
c     10.06.97   ng    included preprocessor option reviewer_costfn,
c                        (calculates chisq as deviations from C rather
c                         than D_C)
cc 
c =====================================================================
c
      double precision function cost_function(
     $     chisq,closure,par_sum,np,
     $     amplitude,amplitude_obs,
     $     c_min,c_max,c_min_obs,c_max_obs,
     $     chisq_term,clos_term,par_term,ampl_term)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none

#include "statistics.h"
#include "sg_params.h"
      integer np
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c 
      double precision scale_fact
#ifdef reviewer_costfn
      parameter (scale_fact = 1.d0)
#else
      parameter (scale_fact = 1.d-3)
#endif
c
      double precision chisq_weight,clos_weight
      double precision par_weight,ampl_weight
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     determine weight of different terms
c --------------------------------------------------------------------
c 
#ifdef reviewer_costfn
      chisq_weight = 1.d0/365.d0
      clos_weight = 0.d0
      par_weight = 1.d0
      ampl_weight = 0.0
#else
      chisq_weight = 1.d0/365.d0
      clos_weight = 1.d0
      par_weight = 1.d0
      ampl_weight = 1.0
#endif
c
c --------------------------------------------------------------------
c     calculate different terms
c --------------------------------------------------------------------
c 
#ifdef reviewer_costfn
      chisq_term = scale_fact * chisq_weight * chisq/((3.07d0)**2)
#else
      chisq_term = scale_fact * chisq_weight * chisq/((sg_dcdt)**2)
#endif
      clos_term = scale_fact * clos_weight * (closure/sg_c)**2
      par_term = scale_fact * par_weight * par_sum
c      ampl_term = scale_fact * ampl_weight * 
c     $     ((amplitude - amplitude_obs)/sg_c)**2
      ampl_term = scale_fact * ampl_weight *
     $     ((c_min - c_min_obs)**2 + (c_max - c_max_obs)**2)/(sg_c)**2
c
c --------------------------------------------------------------------
c     calculate total cost function
c --------------------------------------------------------------------
c 
      cost_function = chisq_term + clos_term + par_term + ampl_term
c
c         write(*,*) '--COST_FUNCTION: value of cost : ',cost_function
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c     
      return

      end











