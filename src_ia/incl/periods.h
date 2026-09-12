c======================= include file periods.h ==============================
c
c variables for the periods that are investigated 
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c 17.10.99 ng        adapted for interannual runs, periods are replaced with
c                       years
c
c variables:
c
c name          type   description
c-----------------------------------------------------------------------------
c 
c nyears        i  s   number of years to analyze
c year_start    dp a   start of analysis year
c 
c  
c------------------------------------------------------------------------
c
      integer nyearmax,nyears
      parameter(nyearmax = 20)
c
      double precision year_start
c
      data year_start / 0.0 /



      

