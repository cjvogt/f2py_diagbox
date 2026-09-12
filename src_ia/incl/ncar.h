c===================== include file ncar.h ===========================
c
c variables for NCAR plotting
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c
c variables:
c
c name          type   description
c-----------------------------------------------------------------------------
c
c plottype      c  s   character specifying plottype 
c                       'sa' : seasonal plot
c                       'ts' : timeseries plot
c
c xautoscale    l  s   logical specifying the use of xplot given by user
c yautoscale    l  s   logical specifying the use of yplot given by user
c yplot         r  s   extrema of y-axis to be plotted
c xplot         r  s   extrema of x-axis to be plotted
c
c colored       l  s   logical specifying if color are used
c
c nticks        i  a   two dimensional array holding the number of tick
c                        marks nticks(axis,type)
c                        dir:   1: x-axis; 2: y-axis
c                        type:  1: major ; 2: minor
c
c tickformat    c  a   array of character holding the tick formats
c
c legend        l  s   logical specifying if a legend is drawn
c legendpos     r  a   array holding the coordinates of the legend
c                        (fractional position within the viewport)
c legendlabel   c  a   array holding the legendlabels
c
c axislabel     c  a   array holding the x- and y-axis labels
c title         c  s   title of plot
c subtitle      c  s   subtitle of plot
c auto_style    l  s   logical indicating if user or automated style for
c                        color and dash pattern of lines is used
c line_color    i  a   array holding colors of lines
c line_dashpat  i  a   array holding dashpatterns
c
c=====================================================================
c
      character plottype*2
      character tickformat(2)*7
      character axislabel(2)*10
      real yplot(2),xplot(2)
c
      logical xautoscale,yautoscale,colored
c
      logical legend
      character legendlabel(nsetmax)*30
      real legendpos(2)
c
      integer nticks(2,2)
c
      character title*60,subtitle*50
c
      logical auto_style
      integer line_color(nsetmax),line_dashpat(nsetmax)
