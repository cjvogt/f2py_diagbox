c =====================================================================
c
c     SUBROUTINE compute_runavg
c
c     PURPOSE: computes the one year running averages of 
c                               the rates and fluxes
c     
c
c     VARIABLES:
c 
c     IN:  tstart          : start time of integration
c          tend            : end time of integration
c          ts              : time step
c          nstep           : number of steps
c          D_C             : rates of change
c          fluxes          : fluxes
c
c
c
c     OUT: D_C_runavg       : average timeseries of dC/dt
c          fluxes_runavg    : average timeseries of fluxes
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     29.08.01   ng     first implementation
c     05.05.02   ng     added running average of anomalies
c 
c =====================================================================
c
      subroutine compute_runavg(tstart,tend,ts,nstep,
     $        D_C,fluxes,
     $        D_C_anom,fluxes_anom,
     $        D_C_runavg,fluxes_runavg,
     $        D_C_a_runavg,fluxes_a_runavg)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
# include "runavg.h"
# include "anomalies.h"
c
      logical error
      integer unit
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      integer k,i,k2,kk
      integer nstep_per_year,nstep_per_halfyr,nstep_start,nstep_stop
      integer nsum(nstepmax)
      double precision sum1,sum2
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c --------------------------------------------------------------------
c     determine the number of timesteps per year
c --------------------------------------------------------------------
c     
      nstep_per_year = int(1.d0/ts)
c
      write(*,*) '--COMPUTE_RUNAVG: number of ts per year : ',
     $     nstep_per_year
c
c --------------------------------------------------------------------
c     determine the first and last steps
c --------------------------------------------------------------------
c     
      nstep_per_halfyr = int(nstep_per_year/2)
      nstep_start = nstep_per_halfyr+1
      nstep_stop = nstep - nstep_per_halfyr
      if (nstep_stop .le. nstep_start) then
         write(*,*) '--COMPUTE_RUNAVG: running average not possible,',
     $        ' because record is not long enough '
         write(*,*) '     nstep_start, nstep_stop : ',nstep_start,
     $        nstep_stop
         error = .true.
         goto 9900
      endif
c
      write(*,*) '--COMPUTE_RUNAVG: nstep_start,nstep_stop : ',
     $     nstep_start,nstep_stop
c
c --------------------------------------------------------------------
c     pad at the beginning
c --------------------------------------------------------------------
c     
      
      do k = 1,nstep_start-1
         do i = 1,nproc_rates
            D_C_runavg(i,k) = 0.0d0
            D_C_a_runavg(i,k) = 0.0d0
         enddo
         do i = 1,nproc_fluxes            
            fluxes_runavg(i,k) = 0.0d0
            fluxes_a_runavg(i,k) = 0.0d0
         enddo
      enddo

c         do i = 1,nproc_rates
c            sum = 0.d0
c            do kk = 1,2*k-1
c               sum = sum + D_C(i,kk)
c            enddo
c            D_C_runavg(i,k) = sum/(2*k-1)
c         enddo
c         do i = 1, nproc_fluxes
c            sum = 0.0d0
c            do kk = 1,2*k-1
c               sum = sum + fluxes(i,kk)
c            enddo
c            fluxes_runavg(i,k) = sum/(2*k-1)
c         enddo
c      enddo
             
c
c --------------------------------------------------------------------
c     do the central (the real) part
c --------------------------------------------------------------------
c     
      do k = nstep_start,nstep_stop
         do i = 1,nproc_rates
            sum1 = 0.d0
            sum2 = 0.d0
             do kk = k-nstep_per_halfyr,k+nstep_per_halfyr
               sum1 = sum1 + D_C(i,kk)
               sum2 = sum2 + D_C_anom(i,kk)
            enddo
            D_C_runavg(i,k) = sum1/(2*nstep_per_halfyr+1)
            D_C_a_runavg(i,k) = sum2/(2*nstep_per_halfyr+1)
         enddo
         do i= 1, nproc_fluxes
            sum1 = 0.0d0
            sum2 = 0.0d0
            do kk = k-nstep_per_halfyr,k+nstep_per_halfyr
               sum1 = sum1 + fluxes(i,kk)
               sum2 = sum2 + fluxes_anom(i,kk)
            enddo
            fluxes_runavg(i,k) = sum1/(2*nstep_per_halfyr+1)
            fluxes_a_runavg(i,k) = sum2/(2*nstep_per_halfyr+1)
         enddo
      enddo
c
c --------------------------------------------------------------------
c     pad at the end
c --------------------------------------------------------------------
c     
      
      do k = nstep_stop+1,nstep
         do i = 1,nproc_rates
            D_C_runavg(i,k) = 0.0d0
            D_C_a_runavg(i,k) = 0.0d0
         enddo
         do i = 1,nproc_fluxes            
            fluxes_runavg(i,k) = 0.0d0
            fluxes_a_runavg(i,k) = 0.0d0
         enddo
      enddo

c         do i = 1,nproc_rates
c            sum = 0.d0
c            do kk = nstep-2*(nstep-k),nstep
c               sum = sum + D_C(i,kk)
c            enddo
c            D_C_runavg(i,k) = sum/(2*(nstep-k)+1)
c         enddo
c         do i = 1, nproc_fluxes
c            sum = 0.0d0
c            do kk = 1,2*k-1
c               sum = sum + fluxes(i,kk)
c            enddo
c            fluxes_runavg(i,k) = sum/(2*(nstep-k)+1)
c         enddo
c      enddo
             
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c
      write(*,*) '--COMPUTE_RUNAVG: terminated.'
c    
      return
c
c --------------------------------------------------------------------
c     error handling
c --------------------------------------------------------------------
c     
 9900 write(*,*) '--COMPUTE_RUNAVG: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











