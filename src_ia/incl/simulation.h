c======================= include file simulation.h ==========================
c
c variables for controlling the simulation
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 15.10.99 ng        adapted for interannual variability runs
c
c variables:
c
c name          type   description
c-----------------------------------------------------------------------------
c 
c tstart        dp s   start time of simulation  [decimal year]
c tend          dp s   end time of simulation    [decimal year]
c ts            dp s   time step of simulation   [decimal year]
c nstepmax      dp s   maximum number of steps  
c nstep         i  s   number of steps
c time          dp s   actual time               [decimal year] 
c
c-------------------------------------------------------------------------
c
      double precision tstart,tend,ts,time
      integer nstep,nstepmax
c
      parameter (nstepmax = 10000)



