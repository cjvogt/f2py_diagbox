c =====================================================================
c
c     SUBROUTINE bkgroundcolor  
c
c     PURPOSE: sets up a the background color and the first foreground
c               color. other colors must be set by a different routine
c
c     VARIABLES:
c     
c     IN: none
c =====================================================================
c
      subroutine bkgroundcolor
c
c ---------------------------------------------------------------------
c     background color white
c ---------------------------------------------------------------------
c
      call gscr(1,0,1.,1.,1.)
c
c ---------------------------------------------------------------------
c     foreground color black
c ---------------------------------------------------------------------
c
      call gscr(1,  1, 0.0,  0.0,  0.0)
c
c ---------------------------------------------------------------------
c     foreground color grey
c ---------------------------------------------------------------------
c
      call gscr(1, 2, 0.7, 0.7, 0.7)

      return
      end
