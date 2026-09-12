c =====================================================================
c
c     SUBROUTINE CALC_RUN_STAT
c 
c     PURPOSE: this routine calculates the chi**2 and other statitical
c                variables that describe the behaviour of the model
c                run
c
c     VARIABLES:
c
c     IN : D_c           : rates of change of C
c          c_sim         : simulated seasonal cycle of C
c          c_obs         : observed seasonal cycle of C
c          tstart        : start time
c          tend          : end time
c          ts            : time step
c          nstep         : number of steps
c          sol_scheme    : solution scheme
c          sum_dc13_prc  : sum of dc13 values of processes
c          nsum_dc13_prc : number of values used for summing sum_dc13_prc
c          sum_fract     : sum for calc of avg fractionation factor
c          sum_kex       : sum for calc of avg kex
c                        
c     OUT: dc13_prc      : dc13 values of processes
c          chisq         : chi**2 of rate of change
c          rsq           : r**2 of rate of change
c          closure       : difference between end and beginning of year
c          amplitude     : calculated amplitude between min and max
c          amplitude_obs : observed amplitude between min and max
c          c_min         : value of minimum C
c          c_max         : value of maximum C
c          c_min_obs     : value of observed minimum C
c          c_max_obs     : value of observed maximum C
c          costfn        : value of cost function
c
c     REVISION:
c
c     date     author  remarks
c   
c      3.08.96   ng    first implementation
c     14.08.96   ng    included costfn, amplitude, amplitude_obs
c      5.09.96   ng    adapted to new version of cost function
c     16.12.96   ng    included avg_kex
c     10.06.97   ng    included preprocessor option reviewer_costfn,
c                        (calculates chisq as deviations from C rather
c                         than D_C)
c     17.10.99   ng    adapting for interannual variability runs
c
c
c =====================================================================
c
      subroutine calc_run_stat(D_C,c_sim,c_obs,tstart,tend,ts,nstep,
     $     sol_scheme,
     $     sum_dc13_prc,nsum_dc13_prc,sum_fract,sum_kex,
     $     dc13_prc,avg_fract,avg_kex,
     $     chisq,rsq,closure,amplitude,amplitude_obs,
     $     c_min,c_max,c_min_obs,c_max_obs,
     $     costfn)
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
# include "statistics.h"
# include "sg_params.h"
# include "averages.h"
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision sum,var,mean
      double precision obs_min,obs_max,sim_min,sim_max
      integer k,i
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision cost_function
c      
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--CALC_RUN_STAT : calculating run statistics...'
c
c --------------------------------------------------------------------
c     calculate mean dc13_prc
c --------------------------------------------------------------------
c        
      do i = 1,4
         if (nsum_dc13_prc(i) .ne. 0) then
            dc13_prc(i) = sum_dc13_prc(i) / 
     $           dble(float(nsum_dc13_prc(i)))
         else
            dc13_prc(i) = -99.99
         endif
      enddo
c
c --------------------------------------------------------------------
c     calculate mean fract_fact
c --------------------------------------------------------------------
c        
      avg_fract = sum_fract / dble(float(nstep))
c
c --------------------------------------------------------------------
c     calculate mean kex
c --------------------------------------------------------------------
c        
      avg_kex = sum_kex / dble(float(nstep))
c
c --------------------------------------------------------------------
c     do detailed analysis for fit in case of standard solution scheme
c --------------------------------------------------------------------
c  
      if (sol_scheme .eq. 'std') then
c
c --------------------------------------------------------------------
c     initialize start values
c --------------------------------------------------------------------
c        
         chisq = 0.d0
         sum = 0.d0
         var = 0.d0
c
c --------------------------------------------------------------------
c     loop over all timesteps, calculate chi**2 and sum for mean
c --------------------------------------------------------------------
c
         do k = 1,nstep
#ifdef reviewer_costfn
            chisq = chisq + (C_sim(k)-C_obs(k))**2
            sum = sum + C_obs(k)
#else
            chisq = chisq + (D_C(6,k)-D_C(7,k))**2
            sum = sum + D_C(7,k)
#endif
         enddo
c
c --------------------------------------------------------------------
c     calculate mean D_C(7,k): should be zero over 365 days
c --------------------------------------------------------------------
c        
         mean = sum/dble(float(nstep))
c
c --------------------------------------------------------------------
c     calculate variance
c --------------------------------------------------------------------
c        
        do k = 1,nstep
#ifdef reviewer_costfn
            var = var + ((D_C(6,k) - mean))**2
#else
            var = var + ((C_sim(k) - mean))**2
#endif
        enddo
c
c --------------------------------------------------------------------
c     calculate r2; r2 = (1 - chi2/var), where
c        var = sum(yi-y_mean)**2 ; variance
c --------------------------------------------------------------------
c
        if (var .ne. 0.0) then
           rsq = 1. - chisq/var
        else
           rsq = 1.
        endif
c
c --------------------------------------------------------------------
c     calculate annual closure [umol/kg]
c --------------------------------------------------------------------
c 
        if ((tend - tstart) .ge. 365.d0) then  ! doesn't work with ia
           k = int(365.d0/ts)
           closure = C_sim(k) - C_sim(1)
        else
           closure = -99.99
        endif
c
c --------------------------------------------------------------------
c     calculate amplitude
c --------------------------------------------------------------------
c 
         sim_max = -9.99e30
         sim_min =  9.99e30
         obs_max = -9.99e30
         obs_min =  9.99e30
c     
         do k = 1,nstep
            if (c_sim(k) .gt. sim_max) then
               sim_max = c_sim(k)
            endif
            if (c_sim(k) .lt. sim_min) then
               sim_min = c_sim(k)
            endif
            if (c_obs(k) .gt. obs_max) then
               obs_max = c_obs(k)
            endif
            if (c_obs(k) .lt. obs_min) then
               obs_min = c_obs(k)
            endif
         enddo
c     
         amplitude = sim_max - sim_min
         amplitude_obs = obs_max - obs_min
         c_min = sim_min
         c_max = sim_max
         c_min_obs = obs_min
         c_max_obs = obs_max
c
c      write(*,*) 'amplitude,amplitude_obs : ',amplitude,amplitude_obs
c
c --------------------------------------------------------------------
c     calculate cost function (without parameter terms)
c --------------------------------------------------------------------
c 
         par_sum = 0.d0
c     
         costfn = cost_function(chisq,closure,par_sum,0,
     $         amplitude,amplitude_obs,
     $         c_min,c_max,c_min_obs,c_max_obs,
     $         chisq_term,clos_term,par_term,ampl_term)
c     
c --------------------------------------------------------------------
c     end of detailed analysis for 'std' solution scheme
c --------------------------------------------------------------------
c     
      endif
c     
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c   
      return
      end

