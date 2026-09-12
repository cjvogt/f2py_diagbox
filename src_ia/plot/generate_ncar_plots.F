c =====================================================================
c
c     SUBROUTINE GENERATE_NCAR_PLOTS
c 
c     PURPOSE: this is the routine for creating ncar plots of the
c                 results
c
c     VARIABLES:
c 
c     IN : tstart          : start time of integration
c          tend            : end time of integration
c          ts              : time step
c          nstep           : number of steps
c          D_c             : rates of change of the different processes
c          fluxes          : fluxes of the different processes
c          D_c_avg         : mean seasonal cycle
c          D_C_anom        : anomalies of the mean seas cycle
c          fluxes_avg      : mean seasonal fluxes
c          fluxes_anom     : anomalies of the mean seas fluxes
c          D_C_runavg      : running average D_C
c          fluxes_runavg   : running average fluxes
c          D_C_a_runavg    : running average D_C anomalies
c          fluxes_a_runavg : running average fluxes anomalies
c          c_sim           : simulated seasonal cycle of c
c          c_obs           : observed seasonal cycle of C
c          exp_name        : experiment name
c          notes           : notes
c     OUT: error           : logical indicating that an error occured
c
c     REVISION:
c
c     date     author  remarks
c   
c     05.08.96   ng    first implementation
c     05.09.96   ng    included possibility to run backwards
c     17.10.99   ng    adapted for interannual variability run, ATTENTION
c                        time units have changed from d-1 to yr-1
c     19.10.99   ng    added plotting of u_var
c     29.08.01   ng    added plotting of anomalies
c                            plotting of running anomalies
c     05.05.02   ng    added plotting of seasonal plots and running
c                            mean of anomalies
c
c =====================================================================
c
      subroutine generate_ncar_plots(tstart,tend,ts,nstep,
     $     D_C,fluxes,     
     $     D_C_avg,D_C_anom,fluxes_avg,fluxes_anom,
     $     D_C_runavg,fluxes_runavg,
     $     D_C_a_runavg,fluxes_a_runavg,
     $     c_sim,c_obs,u_var,exp_name,notes)
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
      implicit none
c
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
# include "solution.h"
# include "anomalies.h"
# include "runavg.h"

# include "plot_sets.h"
# include "ncar.h"
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      character datum*17,datestring*25
c
      real vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt
c
      integer ls,i,k
c
      real xposl(1),xposr(1),ypost(1),yposb(1)
      data xposl(1),xposr(1),yposb(1),ypost(1) /0.15,0.95,0.2,0.8/
c
      real xpos,ypos
c
      double precision sec_day,g_mol,day_year
      parameter (sec_day = 24.d0 * 3600.d0, g_mol = 12.01d0)
      parameter (day_year = 365.d0)
c
      integer nstep_per_year,nstep_per_halfyr
      integer nstep_start,nstep_stop
c
c =================== begin of executable code ========================
c
c
c --------------------------------------------------------------------
c     determine the number of timesteps per year
c --------------------------------------------------------------------
c     
      nstep_per_year = int(1.d0/ts)
      nstep_per_halfyr = int(nstep_per_year/2)
      nstep_start = nstep_per_halfyr+1
      nstep_stop = nstep - nstep_per_halfyr
c
      write(*,*) '--GENERATE_NCAR_PLOTS: '
      write(*,*) '    nstep_per_year, nstep_per_halfyr : ',
     $     nstep_per_year,nstep_per_halfyr
      write(*,*) '    nstep_start, nstep_stop : ',
     $     nstep_start,nstep_stop
c
c ---------------------------------------------------------------------
c     Open GKS and turn clipping off
c ---------------------------------------------------------------------
c
      call opngks
      call gsclip (0)
c
c ---------------------------------------------------------------------
c     sets up a color table (only background and first foreground color)
c ---------------------------------------------------------------------
c
      call bkgroundcolor
c
c---------------------------------------------------------------------
c     timeseries plots should be generated
c---------------------------------------------------------------------
c
      plottype = 'ts'
      xautoscale = .true.
      tickformat(1) = '(f7.0)'
c
c ====================================================================
c     plot the seasonal mean rates
c         plot as a function of month
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      plottype = 'sa'
      xautoscale = .false.
      xplot(1) = 0
      xplot(2) = 365.
      nticks(1,1) = 12
      nticks(1,2) = 3
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nstep_per_year = int(1.d0/ts)
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = nstep_per_year
      enddo
c
      time = (tstart - int(tstart))*365

      do k = 1,nstep_per_year
         do i = 1,5
            xvalplot(i,k) = sngl(time)
            yvalplot(i,k) = sngl(D_C_avg(i,k)/(ts*365.d0))
         enddo
         time = time + ts*365
      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_rates_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting seasonal mean rates'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: adding titles etc...'
c
      title = 'Interannual Diagn. Box Model (IDBM): Seasonal Mean Cycle'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'dC/dt [$G$l$R$mol kg$S$-1$N$ d$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the seasonal mean fluxes
c         plot as a function of month
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      plottype = 'sa'
      xautoscale = .false.
      xplot(1) = 0
      xplot(2) = 365.
      nticks(1,1) = 12
      nticks(1,2) = 3
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nstep_per_year = int(1.d0/ts)
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = nstep_per_year
      enddo
c
      time = (tstart - int(tstart))*365

      do k = 1,nstep_per_year
         do i = 1,5
            xvalplot(i,k) = sngl(time)
            yvalplot(i,k) = sngl(fluxes_avg(i,k) 
     $           * sec_day * day_year)
         enddo
         time = time + ts*365
      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_fluxes_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting data seasonal mean',
     $     ' fluxes...'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: adding titles etc...'
c
      title = 'Interannual Diagn. Box Model (IDBM): Seasonal Mean Cycle'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'fluxes [mol m$S$-2$N$ y$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the interannual rates
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      xautoscale = .true.
      yautoscale = .true.
      auto_style = .true.
      plottype = 'ts'
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = nstep
      enddo
c
      time = tstart

      do k = 1,nstep
         do i = 1,5
            xvalplot(i,k) = sngl(time)
            yvalplot(i,k) = sngl(D_C(i,k)/(ts*365.0d0))
         enddo
         time = time + ts
      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_rates_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting interannual rates...'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: adding titles etc...'
c
      title = 'Interannual Diagn. Box Model (IDBM)'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'dC/dt [$G$l$R$mol kg$S$-1$N$ d$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the fluxes
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = nstep
      enddo
c
      time = tstart
      do k = 1,nstep
         do i = 1,5
            xvalplot(i,k) = sngl(time)
            yvalplot(i,k) = sngl(fluxes(i,k) 
     $           * sec_day * day_year)
         enddo
         time = time + ts

      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_fluxes_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting ',
     $     'interannual fluxes...'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      title = 'Interannual Diagn. Box Model (IDBM)'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'fluxes [mol m$S$-2$N$ y$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the velocity
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nsets = 1
      do i = 1,nsets
         nvalplot(i) = nstep
      enddo
c
      time = tstart

      do k = 1,nstep
         do i = 1,5
            xvalplot(i,k) = sngl(time)
            yvalplot(i,k) = sngl(u_var(k))
         enddo
         time = time + ts
      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      legendlabel(1) = 'horizontal velocity'
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
c      write(*,*) '--GENERATE_NCAR_PLOTS: plotting velocity...'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
c      write(*,*) '--GENERATE_NCAR_PLOTS: adding titles etc...'
c
      title = 'Interannual Diagn. Box Model (IDBM)'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'u [m s$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the rates running averages
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c        don't plot the first half and last half
c---------------------------------------------------------------------
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = 0
      enddo
c
      time = tstart + nstep_per_halfyr*ts

      do k = nstep_start,nstep_stop
         do i = 1,5         
            nvalplot(i) = nvalplot(i) + 1
            xvalplot(i,k-nstep_start+1) = sngl(time)
            yvalplot(i,k-nstep_start+1) = 
     $           sngl(D_C_runavg(i,k)/(ts*365.d0))
         enddo
         time = time + ts
      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_rates_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting ',
     $     'running avg rates...'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: adding titles etc...'
c
      title = 'Interannual Diagn. Box Model (IDBM): Running avg.'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'dC/dt [$G$l$R$mol kg$S$-1$N$ d$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the flux running averages
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = 0
      enddo
c
      time = tstart+nstep_per_halfyr*ts
      do k = nstep_start,nstep_stop
         do i = 1,5
            nvalplot(i) = nvalplot(i) + 1
            xvalplot(i,k-nstep_start+1) = sngl(time)
            yvalplot(i,k-nstep_start+1) = sngl(fluxes_runavg(i,k) 
     $           * sec_day * day_year)
         enddo
         time = time + ts

      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_fluxes_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting ',
     $     'running avg fluxes...'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      title = 'Interannual Diagn. Box Model (IDBM): Running Avg.'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'fluxes [mol m$S$-2$N$ y$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the rates anomalies
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = nstep
      enddo
c
      time = tstart

      do k = 1,nstep
         do i = 1,5
            xvalplot(i,k) = sngl(time)
            yvalplot(i,k) = sngl(D_C_anom(i,k)/(ts*365.d0))
         enddo
         time = time + ts
      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_rates_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting data...'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: adding titles etc...'
c
      title = 'Interannual Diagn. Box Model (IDBM): Anomalies'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'dC/dt [$G$l$R$mol kg$S$-1$N$ d$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the flux anomalies
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = nstep
      enddo
c
      time = tstart
      do k = 1,nstep
         do i = 1,5
            xvalplot(i,k) = sngl(time)
            yvalplot(i,k) = sngl(fluxes_anom(i,k) 
     $           * sec_day * day_year)
         enddo
         time = time + ts

      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_fluxes_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      title = 'Interannual Diagn. Box Model (IDBM): Anomalies'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'fluxes [mol m$S$-2$N$ y$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the running averages of the rate anomalies
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c        don't plot the first half and last half
c---------------------------------------------------------------------
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = 0
      enddo
c
      time = tstart + nstep_per_halfyr*ts

      do k = nstep_start,nstep_stop
         do i = 1,5         
            nvalplot(i) = nvalplot(i) + 1
            xvalplot(i,k-nstep_start+1) = sngl(time)
            yvalplot(i,k-nstep_start+1) = 
     $           sngl(D_C_a_runavg(i,k)/(ts*365.d0))
         enddo
         time = time + ts
      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_rates_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting ',
     $     'running avg rate anomalies...'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: adding titles etc...'
c
      title = 'Interannual Diagn. Box Model (IDBM): Running avg.'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'Anomaly dC/dt [$G$l$R$mol kg$S$-1$N$ d$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the running averages of the flux anomalies
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nsets = 5
      do i = 1,nsets
         nvalplot(i) = 0
      enddo
c
      time = tstart+nstep_per_halfyr*ts
      do k = nstep_start,nstep_stop
         do i = 1,5
            nvalplot(i) = nvalplot(i) + 1
            xvalplot(i,k-nstep_start+1) = sngl(time)
            yvalplot(i,k-nstep_start+1) = sngl(fluxes_a_runavg(i,k) 
     $           * sec_day * day_year)
         enddo
         time = time + ts

      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      do i = 1,5
         legendlabel(i) = proc_fluxes_name(i)
      enddo
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting ',
     $     'running avg flux anomalies...'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      title = 'Interannual Diagn. Box Model (IDBM): Running Avg.'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'Anomaly fluxes [mol m$S$-2$N$ y$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the observed and simulated rates
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c      same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nsets = 2
      do i = 1,nsets
         nvalplot(i) = nstep
      enddo
c
      time = tstart

      do k = 1,nstep
         xvalplot(1,k) = sngl(time)
         yvalplot(1,k) = sngl(D_C(6,k)/(ts*365.d0))
         xvalplot(2,k) = sngl(time)
         yvalplot(2,k) = sngl(D_C(7,k)/(ts*365.d0))
         time = time + ts
      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.05
      legendpos(2) = 0.95
      legendlabel(1) = proc_rates_name(6)
      legendlabel(2) = proc_rates_name(7)
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.2)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting ',
     $     'observed and simulated rates..'

      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      title = 'Interannual Diagn. Box Model (IDBM)'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'dC/dt [$G$l$R$mol kg$S$-1$N$ d$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     advance frame
c---------------------------------------------------------------------
c
      call frame
c
c ====================================================================
c     plot the observed and simulated seasonal cycle of C
c ====================================================================
c
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c       same for auto_style
c---------------------------------------------------------------------
c
      yautoscale = .true.
      auto_style = .true.
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
      nsets = 2
      do i = 1,nsets
         nvalplot(i) = nstep
      enddo
c
      time = tstart

      do k = 1,nstep
         xvalplot(1,k) = sngl(time)
         yvalplot(1,k) = sngl(C_sim(k))
         xvalplot(2,k) = sngl(time)
         yvalplot(2,k) = sngl(C_obs(k))
       time = time + ts

      enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
      legend = .true.
      legendpos(1) = 0.75
      legendpos(2) = 0.95
      legendlabel(1) = 'simulated'
      legendlabel(2) = 'observed'
c
c---------------------------------------------------------------------
c     set other plot parameters
c---------------------------------------------------------------------
c
      colored = .true.
      tickformat(2) = '(f7.1)'
c
c---------------------------------------------------------------------
c     define viewport
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(xposl(1),xposr(1),yposb(1),ypost(1),
     $     wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
      write(*,*) '--GENERATE_NCAR_PLOTS: plotting ',
     $     'observed and simulated concentrations..'
c
      call plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
      title = 'Interannual Diagn. Box Model (IDBM)'
      subtitle = 'Experiment : '//exp_name
      axislabel(2) = 'C [$G$l$R$mol kg$S$-1$N$]'
c
      call pcsetc('FC','$')     ! sets function code character
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
c
      xpos = xposl(1)+(xposr(1)-xposl(1))/2.0
      ypos = ypost(1)+(1.0 - ypost(1))/1.5
      call plchhq(xpos,ypos, title, 0.017, 0.,0.)
      xpos = xposl(1)+0.25*(xposr(1)-xposl(1))
      ypos = ypost(1)+(1.0 - ypost(1))/2.5
      call plchhq(xpos,ypos, subtitle, 0.012, 0.,-1.)
      ypos = ypost(1)+(1.0 - ypost(1))/4.0
      call plchhq(xpos,ypos, notes, 0.012, 0.,-1.)
c
c      ypos = yposb(1)-(0.6*yposb(1))
c      call plchhq(xpos,ypos, axislabel(1), 0.012, 0., -1.)
c
      xpos = 0.0 + 0.2*(xposl(1))
      ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
      call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
      call pcsetc('FC','$')     ! sets function code character
      call getdate(datum)
      datestring = '(ng) '//datum
      call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     close frame and close GKS
c---------------------------------------------------------------------
c
      call frame
      call clsgks
c
c---------------------------------------------------------------------
c     return to calling routine
c---------------------------------------------------------------------
c
      return
c
c
c---------------------------------------------------------------------
c     error handling
c---------------------------------------------------------------------
c
 9900 write(*,9910)
 9910 format(//'An error occured.')
      return

      end
