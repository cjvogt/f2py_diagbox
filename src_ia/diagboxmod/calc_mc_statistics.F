c =====================================================================
c
c     SUBROUTINE CALC_MC_STATISTICS
c 
c     PURPOSE: this routine calculates the final statistics of the
c                monte carlo simulations
c
c     VARIABLES:
c
c     IN : nrun            : number of Monte Carlo runs
c          int_rates_std   : integrated rates of standard run
c          int_fluxes_std  : integrated fluxes of standard run
c          nyears          : number of analysis years
c          costfn_std      : costfn of standard run
c          var_int_rates   : variance of int_rates
c          var_int_fluxes  : variance of int_fluxes
c          wvar_int_rates  : weighted variance of int_rates
c          wvar_int_fluxes : weighted variance of int_fluxes
c          sum_costfn      : sum of cost function
c          mc_d_c          : resampled d_c results of mc runs
c          mc_fluxes       : resampled fluxes results of mc runs
c          mc_c_sim        : resampled c_sim results of mc runs
c          mc_weight       : weight of mc run i
c          nstep2          : number of samples in mc* arrays
c                          
c     OUT: sg_int_rates    : 1-sigma uncertainty of int_rates
c          sg_int_fluxes   : 1-sigma uncertainty of int_fluxes
c          wsg_int_rates   : weighted 1-sigma uncertainty of int_rates
c          wsg_int_fluxes  : weighted 1-sigma uncertainty of int_fluxes
c          mean_costfn     : mean costfunction
c          mean_d_c        : mean d_c of mc runs
c          mean_fluxes     : mean fluxes of mc runs
c          mean_c_sim      : mean c_sim of mc runs
c          wmean_d_c       : weighted mean d_c of mc runs
c          wmean_fluxes    : weighted mean fluxes of mc runs
c          wmean_c_sim     : weighted mean c_sim of mc runs
c          sg_d_c          : 1-sg uncertainty of d_c
c          sg_fluxes       : 1-sg uncertainty of fluxes
c          sg_c_sim        : 1-sg uncertainty of c_sim
c          wsg_d_c         : weighted 1-sg uncertainty of d_c
c          wsg_fluxes      : weighted 1-sg uncertainty of fluxes
c          wsg_c_sim       : weighted 1-sg uncertainty of c_sim
c
c     REVISION:
c
c     date     author  remarks
c   
c     16.08.96   ng    first implementation
c     20.08.96   ng    included weighted and unweighted means and sigmas
c     18.10.99   ng    adapted for interannual variability runs
c
c =====================================================================
c
      subroutine calc_mc_statistics(nrun,     
     $     int_rates_std,int_fluxes_std,nyears,costfn_std,
     $     var_int_rates,var_int_fluxes,
     $     wvar_int_rates,wvar_int_fluxes,
     $     sum_costfn,
     $     mc_D_c,mc_fluxes,mc_c_sim,mc_weight,nstep2,
     $     sg_int_rates,sg_int_fluxes,
     $     wsg_int_rates,wsg_int_fluxes,
     $     mean_costfn,
     $     mean_d_c,mean_fluxes,mean_c_sim,
     $     wmean_d_c,wmean_fluxes,wmean_c_sim,
     $     sg_d_c,sg_fluxes,sg_c_sim,
     $     wsg_d_c,wsg_fluxes,wsg_c_sim)
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
# include "statistics.h"
# include "monte_carlo.h"
# include "mc_statistics.h"
# include "mc_plot_slabs.h"

c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      integer i,j,k,nr
      double precision sum,wsum,ws
c      
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*) '--CALC_MC_STATISTICS: calculating monte carlo',
c     $     ' statistics...'
c
c --------------------------------------------------------------------
c     calculate standard deviation, first calculate variance
c --------------------------------------------------------------------
c
      if (nrun .gt. 1) then
         do i = 1,nyears
            do j = 1,nproc_rates
               var_int_rates(i,j) = var_int_rates(i,j)/
     $              float(nrun-1)
               sg_int_rates(i,j) = dsqrt(var_int_rates(i,j))
               wvar_int_rates(i,j) = wvar_int_rates(i,j)/
     $              float(nrun-1)
               wsg_int_rates(i,j) = dsqrt(wvar_int_rates(i,j))
            enddo
            do j = 1,nproc_fluxes
               var_int_fluxes(i,j) = var_int_fluxes(i,j)/
     $              float(nrun-1)
               sg_int_fluxes(i,j) = dsqrt(var_int_fluxes(i,j))
               wvar_int_fluxes(i,j) = wvar_int_fluxes(i,j)/
     $              float(nrun-1)
               wsg_int_fluxes(i,j) = dsqrt(wvar_int_fluxes(i,j))
            enddo
         enddo
      else
         do i = 1,nyears
            do j = 1,nproc_rates
               sg_int_rates(i,j)   = -99.99
               wsg_int_rates(i,j)  = -99.99
            enddo
            do j = 1,nproc_fluxes
               sg_int_fluxes(i,j)  = -99.99
               wsg_int_fluxes(i,j) = -99.99
            enddo
         enddo
      endif
c
c ---------------------------------------------------------------------
c     for debugging
c ---------------------------------------------------------------------
c
      do i=1,nyears
         write(*,*) (sg_int_fluxes(i,j),j=1,nproc_fluxes)
      enddo
c
c --------------------------------------------------------------------
c     calculate mean cost function
c --------------------------------------------------------------------
c      
      mean_costfn = sum_costfn / float(nrun)
c
c ---------------------------------------------------------------------
c     calculate plot averages of runs
c ---------------------------------------------------------------------
c
c      write(*,*) 'rates...'
      do i=1,nproc_rates
         do k = 1,nstep2
            sum = 0.d0
            wsum = 0.d0
            ws = 0.d0
            do nr = 1,nrun
               sum = sum + mc_d_c(nr,i,k)
               wsum = wsum + mc_d_c(nr,i,k)*mc_weight(nr)
               ws = ws + mc_weight(nr)
            enddo
c
            mean_d_c(i,k) = sum/float(nrun)
            wmean_d_c(i,k) = wsum/ws
         enddo
      enddo
c
c      write(*,*) 'fluxes...'
      do i=1,nproc_fluxes
         do k = 1,nstep2
            sum = 0.d0
            wsum = 0.d0
            ws = 0.d0
            do nr = 1,nrun
               sum = sum + mc_fluxes(nr,i,k)
               wsum = wsum + mc_fluxes(nr,i,k)*mc_weight(nr)
               ws = ws + mc_weight(nr)
            enddo
c
            mean_fluxes(i,k) = sum/float(nrun)
            wmean_fluxes(i,k) = wsum/ws
         enddo
      enddo
c      
c      write(*,*) 'c sim...'
      do k = 1,nstep2
         sum = 0.d0
         wsum = 0.d0
         ws = 0.d0
         do nr = 1,nrun
            sum = sum + mc_c_sim(nr,k)
            wsum = wsum + mc_c_sim(nr,k)*mc_weight(nr)
            ws = ws + mc_weight(nr)
         enddo
c
         mean_c_sim(k) = sum/float(nrun)
         wmean_c_sim(k) = wsum/ws
      enddo
c
c --------------------------------------------------------------------
c     for debugging
c --------------------------------------------------------------------
c 
c      do nr = 1,nrun
c         write(*,*) '--CALC_MC_STATISTICS : mc run # ',nr
c         write(*,*) '   mc_c_sim(nr,k=1,10) : ',
c     $        (mc_c_sim(nr,k),k=1,10)
c      enddo
c      write(*,*) 'mean_c_sim(k=1,10) : ',(mean_c_sim(k),k=1,10)
c      write(*,*) 'wmean_c_sim(k=1,10) : ',(wmean_c_sim(k),k=1,10)
c
c
c ---------------------------------------------------------------------
c     calculate standard deviations
c ---------------------------------------------------------------------
c
c      write(*,*) 'rates...'
      do i=1,nproc_rates
         do k = 1,nstep2
            sum = 0.d0
            wsum = 0.d0
            ws = 0.d0
            do nr = 1,nrun
               sum = sum +
     $              (mc_d_c(nr,i,k)-mean_d_c(i,k))**2
               wsum = wsum + 
     $              ((mc_d_c(nr,i,k)-mean_d_c(i,k))*mc_weight(nr))**2
            enddo
            sg_d_c(i,k) = dsqrt(sum/float(nrun-1))
            wsg_d_c(i,k) = dsqrt(wsum/float(nrun-1))
         enddo
      enddo
c
c      write(*,*) 'fluxes...'
      do i=1,nproc_fluxes
         do k = 1,nstep2
            sum = 0.d0
            wsum = 0.d0
            ws = 0.d0
            do nr = 1,nrun
               sum = sum +
     $              (mc_fluxes(nr,i,k)-mean_fluxes(i,k))**2
               wsum = wsum + 
     $              ((mc_fluxes(nr,i,k)-mean_fluxes(i,k))
     $              *mc_weight(nr))**2
            enddo
            sg_fluxes(i,k) = dsqrt(sum/float(nrun-1))
            wsg_fluxes(i,k) = dsqrt(wsum/float(nrun-1))
         enddo
      enddo
c
c      write(*,*) 'c_sim...'
      do k = 1,nstep2
         sum = 0.d0
         wsum = 0.d0
         ws = 0.d0
         do nr = 1,nrun
            sum = sum +
     $           (mc_c_sim(nr,k)-mean_c_sim(k))**2
            wsum = wsum + 
     $           ((mc_c_sim(nr,k)-mean_c_sim(k))*mc_weight(nr))**2
         enddo
         sg_c_sim(k) = dsqrt(sum/float(nrun-1))
         wsg_c_sim(k) = dsqrt(wsum/float(nrun-1))
      enddo
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

