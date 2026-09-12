c =====================================================================
c
c     SUBROUTINE getdate
c
c     PURPOSE: returns time and date in special format (without seconds)
c              taken from johannes sander.
c     
c     VARIABLES:
c 
c     IN:  
c
c     OUT: datum    : character*17 holding the date
c 
c =====================================================================
c
      subroutine getdate(datum)

      implicit none
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c
      integer i
      character *24 datum1
      character *17 datum
      character *1 datum2(24)
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
      call fdate(datum1)        ! implicit routine
      read(datum1,'(24a)')(datum2(i),i=1,24)
      write(datum,'(17a)')(datum2(i),i=5,16),(datum2(i),i=20,24)
      return
      end
