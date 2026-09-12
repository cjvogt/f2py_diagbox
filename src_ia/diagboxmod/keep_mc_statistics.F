c =====================================================================
c
c     SUBROUTINE KEEP_MC_STATISTICS
c 
c     PURPOSE: this routine updates the statistics of the monte
c                  carlo runs (calculates variance)
c
c     VARIABLES:
c
c     IN : nr              : actual monte carlo run
c          nstep           : number of steps
c          D_c             : array of rates of change
c          fluxes          : array of fluxes
c          c_sim           : array of simulated C
c          int_rates       : array of integrated rates
c          int_fluxes      : array of integrated fluxes
c          nyears          : number of analysis years
c          int_rates_std   : array of integrated rates of standard run
c          int_fluxes_std  : array of integrated fluxes of standard run
c          costfn_std      : standard cost function
c          var_int_rates   : previous variance of int_rates
c          var_int_fluxes  : previous variance of int_fluxes
c          wvar_int_rates  : previous weighted variance of int_rates
c          wvar_int_fluxes : previous weighted variance of int_fluxes
c          chisq           : chi**2 of run
c          closure         : lack of closure of run
c          amplitude       : simulated amplitude
c          amplitude_obs   : observed amplitude    
c          costfn          : costfunction of this run
c          p               : parameters
c          p_init          : initial parameters
c          sg_p            : sigma of parameters   
c          sum_costfn      : previous sum of cost function
c          costfn_min      : previous minimum of cost function
c          costfn_max      : previous maximum of cost function
c                          
c     OUT: mc_D_c          : modified array holding resampled history of D_C
c          mc_fluxes       : modified array holding resampled history of fluxes
c          mc_c_sim        : modified array holding resampled history of c_sim
c          mc_weight       : weight of each run
c          nstep2          : number of samples in mc arrays
c          var_int_rates   : new variance of int_rates
c          var_int_fluxes  : new variance of int_fluxes
c          var_inv_costfn  : new variance of inverse of cost function(weights)
c          wvar_int_rates  : new weighted variance of int_rates
c          wvar_int_fluxes : new weighted variance of int_fluxes
c          sum_costfn      : nev sum of costfn
c          costfn_min      : new minimum of cost function
c          costfn_max      : new maximum of cost function
c
c     REVISION:
c
c     date     author  remarks
c   
c     16.08.96   ng    first implementation
c      2.09.96   ng    changed definition of weight and calc of costfunction
c     18.10.99   ng    adapted for interannual variability runs
c
c =====================================================================
c
      subroutine keep_mc_statistics(
     $     nr,nstep,
     $     D_c,fluxes,c_sim,
     $     int_rates,int_fluxes,nyears,
     $     int_rates_std,int_fluxes_std,costfn_std,
     $     chisq,closure,amplitude,amplitude_obs,costfn,
     $     p,sg_p,p_init,
     $     mc_D_c,mc_fluxes,mc_c_sim,mc_weight,
     $     nstep2,rsmpl_step,
     $     var_int_rates,var_int_fluxes,
     $     wvar_int_rates,wvar_int_fluxes,
     $     sum_costfn,costfn_min,costfn_max)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "simulation.h"
# include "processes.h"
# include "periods.h"
# include "intval.h"
# include "slabs.h"
# include "conj_grad.h"
# include "statistics.h"
# include "monte_carlo.h"
# include "mc_statistics.h"
# include "mc_plot_slabs.h"
c
      integer nr
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer i,j,np,k
      logical first
      data first /.true./
      double precision weight
      integer nsmpl
c
c ---------------------------------------------------------------------
c     save first between calls
c ---------------------------------------------------------------------
c
      save first
c
c ---------------------------------------------------------------------
c     external function
c ---------------------------------------------------------------------
c
      double precision cost_function
c      
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--KEEP_MC_STATISTICS: summing run statistics...'
c
c --------------------------------------------------------------------
c     initialize start values at the first time this subroutine is
c       called
c --------------------------------------------------------------------
c       
      if (first) then
         do i = 1,nyearmax
            do j = 1,nproc_rates
               var_int_rates(i,j) = 0.d0
            enddo
            do j = 1,nproc_fluxes
               var_int_fluxes(i,j) = 0.0d0
            enddo
         enddo
c
         sum_costfn = 0.d0
         costfn_min = 9.99d30
         costfn_max = -9.99d30
c
         first = .false.
c
         rsmpl_step = 40
         if ((nstep / rsmpl_step) .gt. nstepmax2) then
            rsmpl_step = nstep / nstepmax2 + 1
         endif
c
      endif
c
c --------------------------------------------------------------------
c     calculate cost function : changed 2.9.96 ; only cost function
c       without parameter term is calculated
c -------------------------------------------------------------------- 
c
cc
cc --------------------------------------------------------------------
cc     calculate cost function : part of parameter deviation
cc -------------------------------------------------------------------- 
cc
c      par_sum = 0.d0
c      do np = 1,nparams
cc     par_sum = par_sum + p(np)**2
c         par_sum = par_sum + ((p(np) - p_init(np))/sg_p(np))**2
c      enddo
cc     
cc --------------------------------------------------------------------
cc     calculate total cost function
cc --------------------------------------------------------------------
cc 
c      costfn = cost_function(chisq,closure,par_sum,nparams,
c     $     amplitude,amplitude_obs,
c     $     c_min,c_max,c_min_obs,c_max_obs,
c     $     chisq_term,clos_term,par_term,ampl_term)
cc     
      sum_costfn = sum_costfn + costfn
c
c --------------------------------------------------------------------
c     keep minimum and maximum of costfunction
c --------------------------------------------------------------------
c       
      if (costfn .lt. costfn_min) then
         costfn_min = costfn
      endif
      if (costfn .gt. costfn_max) then
         costfn_max = costfn
      endif
c
c --------------------------------------------------------------------
c     add the new values to the old variances
c     for weighted variances : weight each new value
c     use quotient of cost function to cost_function of standard 
c     run as weight
c --------------------------------------------------------------------
c
c      weight = costfn_std/costfn
c      weight = dexp((-costfn/costfn_std + 1.d0))
c      weight = dexp((-costfn/costfn_std + 1.d0)*0.116d0)
#ifdef reviewer_costfn
      weight = dexp(-((costfn-costfn_std)**2)/(4.d0**2.d0))
#else
      weight = dexp(-((costfn-costfn_std)**2)/(0.04d0**2.d0))
#endif
c
      do i = 1,nyears
         do j = 1,nproc_rates
            var_int_rates(i,j) = var_int_rates(i,j) + 
     $           (int_rates(i,j)-int_rates_std(i,j))**2
            wvar_int_rates(i,j) = wvar_int_rates(i,j) +
     $           ((int_rates(i,j)-int_rates_std(i,j))*weight)**2
         enddo
         do j = 1,nproc_fluxes
            var_int_fluxes(i,j) = var_int_fluxes(i,j) + 
     $           (int_fluxes(i,j)-int_fluxes_std(i,j))**2
            wvar_int_fluxes(i,j) = wvar_int_fluxes(i,j) +
     $           ((int_fluxes(i,j)-int_fluxes_std(i,j))*weight)**2
         enddo
      enddo
c
c --------------------------------------------------------------------
c     resample D_C,fluxes and C_sim
c --------------------------------------------------------------------
c 
      nsmpl = 0
c
      do k = 1,nstep
         if (mod(k,rsmpl_step) .eq. 1) then
            nsmpl = nsmpl + 1
            do i = 1,nproc_rates
               mc_d_c(nr,i,nsmpl) = d_c(i,k)
            enddo
            do i = 1,nproc_fluxes
               mc_fluxes(nr,i,nsmpl) = fluxes(i,k)
            enddo
            mc_c_sim(nr,nsmpl) = c_sim(k)
         endif
      enddo
c
      nstep2 = nsmpl
      mc_weight(nr) = weight
c
c --------------------------------------------------------------------
c     for debugging
c --------------------------------------------------------------------
c 
c      write(*,*) '--KEEP_MC_STATISTICS : mc run # ',nr
c      write(*,*) '   c_sim(k=1,10) : ',(c_sim(k),k=1,10)
c      write(*,*) '   mc_c_sim(nr,k=1,10) : ',(mc_c_sim(nr,k),k=1,10)
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

