c =====================================================================
c
c     SUBROUTINE plot_results
c
c     PURPOSE: plots results of diagnostic box model. 
c
c     VARIABLES:
c 
c     IN:  plottype       : type of plot 'sa' or 'ts'
c          nsets          : number of total sets
c          xvalplot       : two dim array holding the x-values to be plotted 
c          yvalplot       : two dim array holding the y-values to be plotted
c          nvalplot       : number of xy-values to be plotted
c          xautoscale     : logical indicating that the x-axis is scaled auto-
c                             matically
c          yautoscale     : dito for y-axis
c          xplot          : extrema of x-axis
c          yplot          : extrema of y-axis
c          nticks         : two dim array holding the number of tick marks
c          tickformat     : character holding the x-tickformat
c          legend         : logical indicating use of legend
c          legendpos      : one dim array holding legend position
c          legendlabel    : one dim array holding legend labels
c          auto_style     : logical indicating if user or automated scheme
c                             is used for line colors and patterns
c          line_dashpat   : array holding dashpatterns
c          line_color     : array holding color patterns
c          colored        : logical specifying if lines are colored
c                         
c     OUT: none
c
c     CHANGES:
c
c     date     author   remarks
c
c     5.08.96   ng      adapted from plot_timeseries.F
c 
c =====================================================================
c
      subroutine plot_results(plottype,nsets,
     $     xvalplot,yvalplot,nvalplot,
     $     xautoscale,yautoscale,xplot,yplot,
     $     nticks,tickformat,
     $     legend,legendpos,legendlabel,
     $     auto_style,line_dashpat,line_color,
     $     colored)

c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c
      implicit none
#include"plot_sets.h"
#include"ncar.h"
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      integer ns,j,i
c
      integer nvalplot1
      real traplot1(nsetmax*nvalplotmax)
      real xvalplot1(nvalplotmax)
      real yvalplot1(nvalplotmax)
      real axismin,axismax
      integer nmajticks,nminticks
c
      real vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt
c
      real xmin,xmax,ymin,ymax,xmintmp,xmaxtmp,ymintmp,ymaxtmp
      real xpos,ypos,mjtickdistx,mjtickdisty
      real xpos2(2),ypos2(2)
      integer xposi
      real ypost,yposb,hobarlength
      real legposuc(2)
c
      integer ls,ic,il
c
      integer dashpat,colindx
c
      character xlab*16, ylab*19, monthlab(12)*3
      character str*40
c
      data monthlab /'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
     $               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c      write(*,*)
c      write(*,*)
c      write(*,*) '--PLOT_RESULTS: '
c      write(*,*) '     number of sets to plot : ',nsets
c
c---------------------------------------------------------------------
c     if autoscale, find xplot, yplot and nticks
c---------------------------------------------------------------------
c
      if (xautoscale) then
c
         nvalplot1 = 0
c
         do ns = 1,nsets
            do j = 1,nvalplot(ns)
               nvalplot1 = nvalplot1 + 1
               traplot1(nvalplot1) = xvalplot(ns,j)
            enddo
         enddo
c
         call calc_grid(traplot1,nvalplot1,
     $        axismin,axismax,nmajticks,nminticks)
c
         xplot(1) = axismin
         xplot(2) = axismax
         nticks(1,1) = nmajticks
         nticks(1,2) = nminticks
c
      endif
c
      if (yautoscale) then
c
         nvalplot1 = 0
c
         do ns = 1,nsets
            do j = 1,nvalplot(ns)
               nvalplot1 = nvalplot1 + 1
               traplot1(nvalplot1) = yvalplot(ns,j)
            enddo
         enddo
c
         call calc_grid(traplot1,nvalplot1,
     $        axismin,axismax,nmajticks,nminticks)
c
         yplot(1) = axismin
         yplot(2) = axismax
         nticks(2,1) = nmajticks
         nticks(2,2) = nminticks
c
      endif
c      
c =====================================================================
c     end of preparations, start with NCAR specific part
c =====================================================================
c
c      write(*,*) '--PLOT_RESULTS: starting to draw '
c
c      write(*,*) '   xplot : ',xplot
c      write(*,*) '   yplot : ',yplot
c
c---------------------------------------------------------------------
c     set up a color table (only colors 2 and higher) if color
c---------------------------------------------------------------------
c
      if (colored) then
         call simplecolor
      endif
c
c---------------------------------------------------------------------
c     get viewport coordinates, afterwards set user coordinates
c---------------------------------------------------------------------
c
      call getset(vpl,vpr,vpb,vpt,wdl,wdr,wdb,wdt,ls)
      call set(vpl,vpr,vpb,vpt,xplot(1),xplot(2),yplot(1),yplot(2),ls)
c
c---------------------------------------------------------------------
c     set tick mark length and direction (outside ticks)
c       major = -12, minor = -8
c---------------------------------------------------------------------
c
      call tick4(-12,-8,-12,-8)
c
c---------------------------------------------------------------------
c     set labeling type to 1, plchhq is used, makes nicer labels
c---------------------------------------------------------------------
c
      call gaseti('LTY - labeling type',1)
c
c --------------------------------------------------------------------
c     make grid using gridal, since ixlb and iylb are zero, labmod
c       does not have to be called
c --------------------------------------------------------------------
c
      call gridal(nticks(1,1),nticks(1,2),nticks(2,1),nticks(2,2),
     $     0,0,5,0.,0.)
c
c---------------------------------------------------------------------
c     force NCAR to use high quality labels for lines
c---------------------------------------------------------------------
c
      call gstxfp(1,2)
c
c---------------------------------------------------------------------
c     draw the the lines
c---------------------------------------------------------------------
c
c      write(*,*) '--PLOT_RESULTS: drawing the lines '
c
CJML      do 4000, ns=1,nsets
C
c
C         do 3100, j=1,nvalplot(ns)
C            xvalplot1(j) = xvalplot(ns,j)
C            yvalplot1(j) = yvalplot(ns,j)
C 3100    continue
c
C         il = 1                 ! connects with lines
c
C         if (colored) then
C            if (auto_style) then
C               colindx = ns+1
C            else
C               colindx = line_color(ns)
C            endif
C         else
C            colindx = 1
C         endif
c
C         if (auto_style) then
C            dashpat = mod(ns-1,5)+1
C         else
C            dashpat = line_dashpat(ns)
C         endif
c
C         call gsplci(colindx)   ! select color of lines
C         call gsln(dashpat)     ! set dashpat
Cc         str = '$$$$$$$$$$$'//legendlabel(ns)//'$$$$$$$'
Cc         call dashdc(str,15,20) ! adds legendlabel to line
C         call curved(xvalplot1,yvalplot1,nvalplot(ns))
C 4000 continue
      do 4000, ns = 1, nsets
c        gather only finite points for this set
         nvalplot1 = 0
         do 3100, j = 1, nvalplot(ns)
c           filter NaNs/Infs (NaN .ne. NaN)
            if (xvalplot(ns,j) .eq. xvalplot(ns,j) .and.
     $          yvalplot(ns,j) .eq. yvalplot(ns,j)) then
               nvalplot1 = nvalplot1 + 1
               xvalplot1(nvalplot1) = xvalplot(ns,j)
               yvalplot1(nvalplot1) = yvalplot(ns,j)
            endif
 3100    continue

c        only draw if we have at least 2 valid points
         if (nvalplot1 .ge. 2) then
            il = 1

            if (colored) then
               if (auto_style) then
                  colindx = ns + 1
               else
                  colindx = line_color(ns)
               endif
            else
               colindx = 1
            endif

            if (auto_style) then
               dashpat = mod(ns-1,5)+1
            else
               dashpat = line_dashpat(ns)
            endif

            call gsplci(colindx)
            call gsln(dashpat)
            call curved(xvalplot1, yvalplot1, nvalplot1)
         endif
 4000 continue

c
c---------------------------------------------------------------------
c     draw the tick labels
c---------------------------------------------------------------------
c
c      write(*,*) '--PLOT_RESULTS: drawing the tick labels '
      call gsln(1)              ! select dash pattern solid
      call gstxci(1)            ! select black color for points
      call gsplci(1)            ! select black color for lines
c
c---------------------------------------------------------------------
c     draw the x-tick labels for timeseries plot
c---------------------------------------------------------------------
c
c      write(*,*) '--PLOT_RESULTS: drawing x-tick labels...'
c      write(*,*) '  format = ',tickformat(1)
      ypos = yplot(1) - 0.05*(yplot(2)-yplot(1))
      mjtickdistx = (xplot(2)-xplot(1))/max(1,nticks(1,1))
c
      if (plottype .eq. 'ts') then
         do i=1,(nticks(1,1))
c            xpos = xplot(1) + (i-0.5)*mjtickdistx
            xpos = xplot(1) + (i-1)*mjtickdistx
c            xposi = int(xpos)
            write(xlab,tickformat(1)) xpos
            call plchhq(xpos,ypos,xlab,0.011,0.,0.)
         enddo
c
c---------------------------------------------------------------------
c     draw the x-tick labels for seasonally aggregated plot
c---------------------------------------------------------------------
c
      elseif (plottype .eq. 'sa') then
         do 5100, i=1,nticks(1,1)
            xpos = xplot(1) + (i-0.5)*mjtickdistx
            call plchhq(xpos,ypos,monthlab(i),0.011,0.,0.)
 5100    continue
      endif
c
c---------------------------------------------------------------------
c     draw the y-tick labels
c---------------------------------------------------------------------
c
c      write(*,*) '--PLOT_RESULTS: drawing y-tick labels...'
c      write(*,*) '  format = ',tickformat(2)
      xpos = xplot(1) - 0.05*(xplot(2)-xplot(1))
      mjtickdisty = (yplot(2)-yplot(1))/max(1,nticks(2,1))
                  
      do 5300, i=1,(nticks(2,1)+1)
         ypos = yplot(1) + (i-1)*mjtickdisty
         write(ylab,tickformat(2)) ypos
         call plchhq(xpos,ypos,ylab,0.011,0.,0.)
 5300 continue
     
c
c---------------------------------------------------------------------
c     draw the legend if desired
c---------------------------------------------------------------------
c
      if (legend) then
c         write(*,*) '--PLOT_RESULTS: drawing legend'
         call pcsetc('FC','$')  ! sets function code character
         legposuc(1) = xplot(1)+legendpos(1)*(xplot(2)-xplot(1))
         legposuc(2) = yplot(1)+legendpos(2)*(yplot(2)-yplot(1))
c
         do ns = 1,nsets
            if (legendlabel(ns) .ne. '') then
               xpos2(1) = legposuc(1)
               ypos2(1) = legposuc(2)-(ns-1)*0.05*(yplot(2)-yplot(1))
               xpos2(2) = legposuc(1)+0.05*(xplot(2)-xplot(1))
               ypos2(2) = ypos2(1)
c     
               if (colored) then ! set color
                  if (auto_style) then
                     colindx = ns + 1
                  else
                     colindx = line_color(ns)
                  endif
                  call gsplci (colindx) ! select color of lines
               endif
c     
               if (auto_style) then
                  dashpat = mod(ns-1,5)+1
               else
                  dashpat = line_dashpat(ns)
               endif
               call gsln(dashpat) ! set dashpatc
               call curved(xpos2,ypos2,2)
               if (colored) then ! set color back to black
                  call gstxci(1)
                  call gsplci (1)
               endif
               call gsln(1)
               xpos = xpos2(2) + 0.05*(xplot(2)-xplot(1))
               ypos = ypos2(1)
               call plchhq(xpos,ypos,legendlabel(ns),0.009,0.,-1.)
            endif
         enddo
         
      else
      endif
      
c
c --------------------------------------------------------------------
c     return to calling program
c --------------------------------------------------------------------
c
c      write(*,*)'PLOT_RESULTS: successfully completed '
      return
      end
