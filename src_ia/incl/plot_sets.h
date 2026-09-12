c===================== include file plotsets.h ===========================
c
c variables containing the values for plotting
c
c REVISIONS:
c
c date     author    remarks
c
c 05.08.96  ng       first implementation
c
c variables:
c
c name           type   description
c---------------------------------------------------------------------
c
c xvalplot       r  a   two dim array holding the x-values to be 
c                        plotted (xvalplot(nsets,ntraplot))  
c yvalplot       r  a   two dim array holding the y-values to be  
c                        plotted (yvalplot(nsets,ntraplot))
c nvalplotmax    i  s   maximum number of points to plot per set
c nvalplot       i  a   one dim array holding the sizes of the sets 
c                        to be plotted
c nsetmax        i  s   maximum number of sets
c nsets          i  s   number of sets  
c
c=====================================================================
c
      integer nvalplotmax,nsetmax,nsets
      parameter(nvalplotmax = 8000,nsetmax = 15)
c
      real xvalplot(nsetmax,nvalplotmax)
      real yvalplot(nsetmax,nvalplotmax)
c
      integer nvalplot(nsetmax)

