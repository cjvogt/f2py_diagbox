c =====================================================================
c
c     SUBROUTINE simplecolor
c
c     PURPOSE: sets up a simple color table starting from index 2.
c              background and first foreground color must be set
c              by a different routine
c
c     VARIABLES:
c 
c     IN:      none 
c     OUT:     none
c
c     CHANGES:
c
c     date     author   remarks
c
c     22.12.94   ng     first implementation
c     19.08.96   ng     special color table for diag_box_model
c 
c =====================================================================
c
      subroutine simplecolor

c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c red
      call gscr(1,  2, 1.0, 0.0, 0.0)
c blue
      call gscr(1,  3, 0.0, 0.0, 1.0)
c darkviolet
      call gscr(1,  4, 0.6, 0.0, 0.8)
c purple
      call gscr(1,  5, 1.0, 0.0, 1.0)
c green
      call gscr(1,  6, 0.0, 1.0, 0.0)
c yellow
      call gscr(1,  7, 1.0, 1.0, 0.0)
c greenyellow
      call gscr(1,  8, 0.7, 1.0, 0.2)
c chartreuse
      call gscr(1,  9, 0.5, 1.0, 0.0)
c celeste
      call gscr(1, 10, 0.2, 1.0, 0.5)
c orangered
      call gscr(1, 11, 1.0, 0.0, 0.2)
c aqua
      call gscr(1, 12, 0.0, 0.9, 1.0)
c deepskyblue
      call gscr(1, 13, 0.0, 0.75, 1.0)
c royalblue
      call gscr(1, 14, 0.25, 0.45, 0.95)
c slateblue
      call gscr(1, 15, 0.4, 0.35, 0.8)
c gold
      call gscr(1, 16, 1.0, 0.85, 0.0)
c lavender
      call gscr(1, 17, 0.8, 0.8, 1.0)
c orchid
      call gscr(1, 18, 0.85, 0.45, 0.8)
c orange
      call gscr(1, 19, 1.0, 0.65, 0.0)
c ??
c
c --------------------------------------------------------------------
c     return to calling program
c --------------------------------------------------------------------
c
      return
      end
