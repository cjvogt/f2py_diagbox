c =====================================================================
c
c     SUBROUTINE GENERATE_MC_NCAR_PLOTS
c 
c     PURPOSE: this is the routine for creating ncar plots of the
c                 results of the monte carlo runs. 
c
c     VARIABLES:
c 
c     IN : nrun            : number of total monte carlo runs
c          tstart          : start time of integration
c          tend            : end time of integration
c          ts              : time step
c          nstep           : number of steps
c          D_c_std         : rates of change of the different processes
c          fluxes_std      : fluxes of the different processes
c          c_sim_std       : simulated seasonal cycle of c
c          c_obs           : observed seasonal cycle of C
c          exp_name        : experiment name
c          notes           : notes
c
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
c          nstep2          : number of samples in sg* arrays
c          rsmpl_step      : resampling time step of sg* arrays
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISION:
c
c     date     author  remarks
c   
c     19.08.96   ng    first implementation based on generate_ncar_plots
c     17.10.99   ng    adapted for interannual variability run, ATTENTION
c                        time units have changed from d-1 to yr-1
c
c =====================================================================
c
      subroutine generate_mc_ncar_plots(nrun,
     $     tstart,tend,ts,nstep,
     $     D_C_std,fluxes_std,c_sim_std,c_obs,
     $     mean_d_c,mean_fluxes,mean_c_sim,
     $     wmean_d_c,wmean_fluxes,wmean_c_sim,
     $     sg_D_c,sg_fluxes,sg_c_sim,
     $     wsg_D_c,wsg_fluxes,wsg_c_sim,
     $     nstep2,rsmpl_step,
     $     exp_name,notes)
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
# include "monte_carlo.h"
# include "mc_plot_slabs.h"

# include "plot_sets.h"
# include "ncar.h"
c
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      character datum*17,datestring*25
c
      real vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt
c
      integer ls,i,k,k2,nr,l
c
      real xposl(1),xposr(1),ypost(1),yposb(1)
      data xposl(1),xposr(1),yposb(1),ypost(1) /0.15,0.95,0.2,0.8/
c
      real xpos,ypos
c
      double precision sec_day,day_year,g_mol
      parameter (sec_day = 24.d0 * 3600.d0, g_mol = 12.01d0)
      parameter (day_year = 365.)
c
c
c
c =================== begin of executable code ========================
c
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
c     determine whether seasonal or a timeseries plots
c        should be generated
c---------------------------------------------------------------------
c
      plottype = 'ts'
      xautoscale = .true.
      tickformat(1) = '(f7.0)'
c
c ====================================================================
c     plot the rates with weighted and unweighted uncertainties
c ====================================================================
c
c      write(*,*) '--GENERATE_MC_NCAR_PLOTS: plotting rates...'
c     
c---------------------------------------------------------------------
c     set xautoscale to false for the moment
c---------------------------------------------------------------------
c
      yautoscale = .true.
c     
c---------------------------------------------------------------------
c     loop over weighted and unweighted uncertainties
c---------------------------------------------------------------------
c
      do l=1,2
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
         nsets = 15
         do i = 1,nsets
            nvalplot(i) = nstep2
         enddo
c     
         time = tstart
         do k = 1,nstep2
            do i = 1,5
               xvalplot(i,k) = sngl(time)
               if (l .eq. 1) then
                  yvalplot(i,k) = sngl(mean_D_C(i,k)/ts)
               elseif (l .eq. 2) then
                  yvalplot(i,k) = sngl(wmean_D_C(i,k)/ts)
               endif
            enddo
            do i = 6,10
               xvalplot(i,k) = sngl(time)
               if (l .eq. 1) then
                  yvalplot(i,k) = 
     $                 sngl((mean_D_C(i-5,k)+sg_d_c(i-5,k))/ts)
               elseif (l .eq. 2) then
                  yvalplot(i,k) = 
     $                 sngl((wmean_D_C(i-5,k)+wsg_d_c(i-5,k))/ts)
               endif                  
            enddo
            do i = 11,15
               xvalplot(i,k) = sngl(time)
               if (l .eq. 1) then
                  yvalplot(i,k) = 
     $                 sngl((mean_D_C(i-10,k)-sg_d_c(i-10,k))/ts)
               elseif (l .eq. 2) then
                  yvalplot(i,k) = 
     $                 sngl((wmean_D_C(i-10,k)-wsg_d_c(i-10,k))/ts)
               endif
            enddo
            time = time + float(rsmpl_step)*ts
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
         do i = 6,15
            legendlabel(i) = ''
         enddo
c
c---------------------------------------------------------------------
c     determine the line and color style
c---------------------------------------------------------------------
c
         auto_style = .false.
         line_dashpat(1) = 1
         line_dashpat(2) = 2
         line_dashpat(3) = 3
         line_dashpat(4) = 4
         line_dashpat(5) = 1
         do i=6,15
            line_dashpat(i) = 2
         enddo
         line_color(1) = 2
         line_color(2) = 3
         line_color(3) = 4
         line_color(4) = 5
         line_color(5) = 6
         do i=6,15
            line_color(i) = line_color(mod(i-1,5)+1)
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
     $        wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
         call plot_results(plottype,nsets,
     $        xvalplot,yvalplot,nvalplot,
     $        xautoscale,yautoscale,xplot,yplot,
     $        nticks,tickformat,
     $        legend,legendpos,legendlabel,
     $        auto_style,line_dashpat,line_color,
     $        colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
         title = 'Interannual Diagnostic Box Model (IDBM)'

         if (l .eq. 1) then
            subtitle = 'unweighted uncertainties; Exp. : '//exp_name
         elseif (l .eq. 2) then
            subtitle = 'weighted uncertainties; Exp. : '//exp_name
         endif
            
         axislabel(2) = 'dC/dt [$G$l$R$mol kg$S$-1$N$ d$S$-1$N$]'
c     
         call pcsetc('FC','$')  ! sets function code character
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
c     ypos = yposb(1)-(0.6*yposb(1))
c     call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c     
         xpos = 0.0 + 0.2*(xposl(1))
         ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
         call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c     
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
         call pcsetc('FC','$')  ! sets function code character
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
c---------------------------------------------------------------------
c     end of loop
c---------------------------------------------------------------------
c
      enddo
c
c ====================================================================
c     plot the fluxes
c ====================================================================
c
c      write(*,*) '--GENERATE_MC_NCAR_PLOTS: plotting fluxes...'
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c---------------------------------------------------------------------
c
      yautoscale = .true.
c     
c---------------------------------------------------------------------
c     loop over weighted and unweighted uncertainties
c---------------------------------------------------------------------
c
      do l=1,2
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
         nsets = 15
         do i = 1,nsets
            nvalplot(i) = nstep2
         enddo
c     
         time = tstart
         do k = 1,nstep2
            do i = 1,5
               xvalplot(i,k) = sngl(time)
               if (l .eq. 1) then
                  yvalplot(i,k) = sngl(mean_fluxes(i,k) 
     $                 * sec_day * day_year * g_mol)
               elseif (l .eq. 2) then
                  yvalplot(i,k) = sngl(wmean_fluxes(i,k) 
     $                 * sec_day * day_year * g_mol)  
               endif
            enddo
            do i = 6,10
               xvalplot(i,k) = sngl(time)
               if (l .eq. 1) then
                  yvalplot(i,k) = 
     $                 sngl((mean_fluxes(i-5,k)+sg_fluxes(i-5,k))
     $                 *sec_day* day_year * g_mol)
               elseif (l .eq. 2) then
                  yvalplot(i,k) = 
     $                 sngl((wmean_fluxes(i-5,k)+wsg_fluxes(i-5,k))
     $                 *sec_day*day_year*g_mol)
               endif
            enddo
            do i = 11,15
               xvalplot(i,k) = sngl(time)
               if (l .eq. 1) then
                  yvalplot(i,k) = 
     $                 sngl((mean_fluxes(i-10,k)-sg_fluxes(i-10,k))
     $                 *sec_day*day_year*g_mol)
               elseif (l .eq. 2) then
                  yvalplot(i,k) = 
     $                 sngl((wmean_fluxes(i-10,k)-wsg_fluxes(i-10,k))
     $                 *sec_day*day_year*g_mol)
               endif
            enddo
            time = time + float(rsmpl_step)*ts
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
         do i = 6,15
            legendlabel(i) = ''
         enddo
c
c---------------------------------------------------------------------
c     determine the line and color style
c---------------------------------------------------------------------
c
         auto_style = .false.
         line_dashpat(1) = 1
         line_dashpat(2) = 2
         line_dashpat(3) = 3
         line_dashpat(4) = 4
         line_dashpat(5) = 1
         do i=6,15
            line_dashpat(i) = 2
         enddo
         line_color(1) = 2
         line_color(2) = 3
         line_color(3) = 4
         line_color(4) = 5
         line_color(5) = 6
         do i=6,15
            line_color(i) = line_color(mod(i-1,5)+1)
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
     $        wdl,wdr,wdb,wdt,ls)
c     
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
         call plot_results(plottype,nsets,
     $        xvalplot,yvalplot,nvalplot,
     $        xautoscale,yautoscale,xplot,yplot,
     $        nticks,tickformat,
     $        legend,legendpos,legendlabel,
     $        auto_style,line_dashpat,line_color,
     $        colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
         title = 'Seasonal Diagnostic Box Model (SDBM)'
c
         if (l .eq. 1) then
            subtitle = 'unweighted uncertainties; Exp. : '//exp_name
         elseif (l .eq. 2) then
            subtitle = 'weighted uncertainties; Exp. : '//exp_name
         endif
c
         axislabel(2) = 'fluxes [gC m$S$-2$N$ d$S$-1$N$]'
c     
         call pcsetc('FC','$')  ! sets function code character
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
c     ypos = yposb(1)-(0.6*yposb(1))
c     call plchhq(xpos,ypos, axislabel(1), 0.012, 0., 0.)
c     
         xpos = 0.0 + 0.2*(xposl(1))
         ypos = yposb(1) + (ypost(1)-yposb(1))/2.0
         call plchhq(xpos,ypos,axislabel(2), 0.012, 90., 0.)
c
c---------------------------------------------------------------------
c     draw date
c---------------------------------------------------------------------
c
         call pcsetc('FC','$')  ! sets function code character
         call getdate(datum)
         datestring = '(ng) '//datum
         call plchhq(.1,.03,datestring,0.007,0.,-1.)
c
c---------------------------------------------------------------------
c     call frame
c---------------------------------------------------------------------
c
         call frame
c
c---------------------------------------------------------------------
c     end of loop
c---------------------------------------------------------------------
c
      enddo
c
c ====================================================================
c     plot the observed and simulated rates
c ====================================================================
c
c      write(*,*) '--GENERATE_MC_NCAR_PLOTS: plotting obs rates...'
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c---------------------------------------------------------------------
c
      yautoscale = .true.
c     
c---------------------------------------------------------------------
c     loop over weighted and unweighted uncertainties
c---------------------------------------------------------------------
c
      do l=1,2
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
         nsets = 4
         do i = 1,nsets
            nvalplot(i) = nstep2
         enddo
c     
         time = tstart
         do k = 1,nstep2
            xvalplot(1,k) = sngl(time)
            xvalplot(2,k) = sngl(time)
            xvalplot(3,k) = sngl(time)
            xvalplot(4,k) = sngl(time)
            yvalplot(2,k) = sngl(mean_D_C(7,k)/ts)
            if (l .eq. 1) then
               yvalplot(1,k) = sngl(mean_D_C(6,k)/ts)
               yvalplot(3,k) = sngl((mean_d_c(6,k)+sg_d_c(6,k))/ts)
               yvalplot(4,k) = sngl((mean_d_c(6,k)-sg_d_c(6,k))/ts)
            elseif (l .eq. 2) then
               yvalplot(1,k) = sngl(wmean_D_C(6,k)/ts)
               yvalplot(3,k) = sngl((wmean_d_c(6,k)+wsg_d_c(6,k))/ts)
               yvalplot(4,k) = sngl((wmean_d_c(6,k)-wsg_d_c(6,k))/ts)
            endif
            time = time + float(rsmpl_step) * ts
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
         do i = 3,4
            legendlabel(i) = ''
         enddo
c
c---------------------------------------------------------------------
c     determine the line and color style
c---------------------------------------------------------------------
c
         auto_style = .false.
         line_dashpat(1) = 1
         line_dashpat(2) = 1
         line_dashpat(3) = 2
         line_dashpat(4) = 2
         line_color(1) = 2
         line_color(2) = 3
         line_color(3) = 2
         line_color(4) = 2
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
     $        wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
         call plot_results(plottype,nsets,
     $        xvalplot,yvalplot,nvalplot,
     $        xautoscale,yautoscale,xplot,yplot,
     $        nticks,tickformat,
     $        legend,legendpos,legendlabel,
     $        auto_style,line_dashpat,line_color,
     $        colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
         title = 'Seasonal Diagnostic Box Model (SDBM)'
c
         if (l .eq. 1) then
            subtitle = 'unweighted uncertainties; Exp. : '//exp_name
         elseif (l .eq. 2) then
            subtitle = 'weighted uncertainties; Exp. : '//exp_name
         endif
c
         axislabel(2) = 'dC/dt [$G$l$R$mol kg$S$-1$N$ d$S$-1$N$]'
c     
         call pcsetc('FC','$')  ! sets function code character
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
         call pcsetc('FC','$')  ! sets function code character
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
c---------------------------------------------------------------------
c     end of loop
c---------------------------------------------------------------------
c
      enddo
c
c ====================================================================
c     plot the observed and simulated seasonal cycle of C
c ====================================================================
c
c      write(*,*) '--GENERATE_MC_NCAR_PLOTS: plotting C and C_sim...'
c
c---------------------------------------------------------------------
c     set yautoscale to true for the moment
c---------------------------------------------------------------------
c
      yautoscale = .true.
c     
c---------------------------------------------------------------------
c     loop over weighted and unweighted uncertainties
c---------------------------------------------------------------------
c
      do l=1,2
c
c---------------------------------------------------------------------
c     prepare the sets for plotting of the rates of change
c        convert to single precision for plotting
c---------------------------------------------------------------------
c
         nvalplot(1) = nstep
         nvalplot(2) = nstep
         nvalplot(3) = nstep2
         nvalplot(4) = nstep2
         nvalplot(5) = nstep2
c
         nsets = 5
c     
         time = tstart
         do k = 1,nstep
            xvalplot(1,k) = sngl(time)
            yvalplot(1,k) = sngl(c_sim_std(k))
            xvalplot(2,k) = sngl(time)
            yvalplot(2,k) = sngl(C_obs(k))
            time = time + ts
         enddo
c     
         time = tstart
         do k = 1,nstep2
            k2 = (k-1)*rsmpl_step + 1

            xvalplot(3,k) = sngl(time)
            xvalplot(4,k) = sngl(time)
            xvalplot(5,k) = sngl(time)
            if (l .eq. 1) then
               yvalplot(4,k) = sngl(c_sim_std(k2)+sg_c_sim(nstep2))
               yvalplot(5,k) = sngl(c_sim_std(k2)-sg_c_sim(nstep2))
               yvalplot(3,k) = sngl(mean_c_sim(k))
            elseif (l .eq. 2) then
               yvalplot(3,k) = sngl(wmean_c_sim(k))
               yvalplot(4,k) = sngl(c_sim_std(k2)+wsg_c_sim(nstep2))
               yvalplot(5,k) = sngl(c_sim_std(k2)-wsg_c_sim(nstep2))
            endif
              
            time = time + float(rsmpl_step)*ts
         enddo
c
c---------------------------------------------------------------------
c     determine the legend
c---------------------------------------------------------------------
c
         legend = .true.
         legendpos(1) = 0.70
         legendpos(2) = 0.95
         legendlabel(1) = 'simulated std'
         legendlabel(2) = 'observed'
         legendlabel(3) = 'simulated mean'
         do i = 4,5
            legendlabel(i) = ''
         enddo
         
c
c---------------------------------------------------------------------
c     determine the line and color style
c---------------------------------------------------------------------
c
         auto_style = .false.
         line_dashpat(1) = 1
         line_dashpat(2) = 4
         line_dashpat(3) = 1
         line_dashpat(4) = 2
         line_dashpat(5) = 2
         line_color(1) = 3
         line_color(2) = 2
         line_color(3) = 4
         line_color(4) = 4
         line_color(5) = 4
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
     $        wdl,wdr,wdb,wdt,ls)
c
c---------------------------------------------------------------------
c     call routine which does the plotting 
c---------------------------------------------------------------------
c
         call plot_results(plottype,nsets,
     $        xvalplot,yvalplot,nvalplot,
     $        xautoscale,yautoscale,xplot,yplot,
     $        nticks,tickformat,
     $        legend,legendpos,legendlabel,
     $        auto_style,line_dashpat,line_color,
     $        colored)
c
c---------------------------------------------------------------------
c     draw title and axis labels
c---------------------------------------------------------------------
c
         title = 'Seasonal Diagnostic Box Model (SDBM)'
c
         if (l .eq. 1) then
            subtitle = 'unweighted uncertainties; Exp. : '//exp_name
         elseif (l .eq. 2) then
            subtitle = 'weighted uncertainties; Exp. : '//exp_name
         endif
c
         axislabel(2) = 'C [$G$l$R$mol kg$S$-1$N$]'
c     
         call pcsetc('FC','$')  ! sets function code character
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
         call pcsetc('FC','$')  ! sets function code character
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
c---------------------------------------------------------------------
c     end of loop
c---------------------------------------------------------------------
c     
      enddo
c
c---------------------------------------------------------------------
c     close frame and close GKS if last run
c---------------------------------------------------------------------
c
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
