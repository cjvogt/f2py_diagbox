c =====================================================================
c
c     SUBROUTINE WRITE_OBS_TOFILE
c
c     PURPOSE: writes the observations to a file specified by
c                the unit number unit
c     
c
c     VARIABLES:
c 
c     IN:  unit            : unit number to write to
c          tstart          : start time of integration
c          ts              : time step
c          nstep           : number of steps
c          c_obs           : array holding the observed seasonal cycle of C
c          dc13_obs        : array holding the observed seasonal cycle of dC13
c          h_obs           : array holding the observed seasonal cycle of h
c          temp_obs        : observed seasonal cycle of temp
c          sal_obs         : observed seasonal cycle of sal
c          salk_obs        : observed seasonal cycle of salk
c          fco2_o_obs      : observed seasonal cycle of oceanic fco2
c          pco2_a_obs      : observed seasonal cycle of atmospheric xco2
c          fco2_a_obs      : observed seasonal cycle of atmospheric fco2
c          pco2_a_dry      : observed atmospheric pCO2
c          dc13_a_obs      : observed seasonal cycle of atmospheric dc13
c          ws_obs          : observed seasonal cycle of windspeed
c          dcdz_obs        : observed seasonal cycle of dcdz
c          kz_obs          : observed seasonal cycle of Kz
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     05.08.96   ng     first implementation
c     08.08.96   ng     included temp_obs,sal_obs etc
c     07.06.00   ng     adopted for interannual var run
c                         added salk_obs
c     26.10.02   ng     added atmospheric pCO2
c 
c =====================================================================
c
      subroutine write_obs_tofile(unit,tstart,ts,nstep,
     $     c_obs,dc13_obs,h_obs,
     $     temp_obs,sal_obs,
     $     salk_obs,
     $     fco2_o_obs,pco2_a_obs,fco2_a_obs,pco2_a_dry,
     $     dc13_a_obs,
     $     ws_obs,dcdz_obs,kz_obs,
     $     error)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
# include "obs_values.h"
c
      logical error
      integer unit
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      integer k
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      error = .false.
      time = tstart
c
c --------------------------------------------------------------------
c     loop over all time steps
c --------------------------------------------------------------------
c     
      do k = 1,nstep
c --------------------------------------------------------------------
c     write observations to file
c --------------------------------------------------------------------
c
         write(unit,9000) 
     $        time,
     $        c_obs(k),dc13_obs(k),h_obs(k),
     $        temp_obs(k),sal_obs(k),salk_obs(k),
     $        fco2_o_obs(k),pco2_a_obs(k),fco2_a_obs(k),pco2_a_dry(k),
     $        dc13_a_obs(k),
     $        ws_obs(k),dcdz_obs(k),kz_obs(k)
c
 9000    format(f10.3,1x,f6.1,1x,f7.3,1x,f5.1,1x,
     $        f5.2,1x,f6.3,1x,f6.1,1x,
     $        4(f5.1,1x),f7.3,1x,
     $        f5.2,1x,f5.3,1x,e10.4)
c
c --------------------------------------------------------------------
c     increase time
c --------------------------------------------------------------------
c
         time = time + ts
c
c --------------------------------------------------------------------
c     end of loop over all simulated values
c --------------------------------------------------------------------
c
      enddo
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c     
      return
c
c --------------------------------------------------------------------
c     error handling
c --------------------------------------------------------------------
c     
 9900 write(unit,*) '--WRITE_OBS_TOFILE: ',
     $     'an error occurred while writing...'
      error = .true.
c
      return

      end











