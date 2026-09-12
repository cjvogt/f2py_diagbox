c======================= include file processes.h ============================
c
c variables for the processes that are investigated 
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
c nproc_fluxes      i  s   number of processes for the fluxes
c nproc_rates       i  s   number of processes for the rates of change
c proc_fluxes_name  c  a   array holding the names of the flux processes
c proc_rates_name   c  a   array holding the names of the rates processes
c  
c------------------------------------------------------------------------
c
      integer nproc_fluxes,nproc_rates
      parameter(nproc_fluxes = 10, nproc_rates = 7)
c
      character proc_fluxes_name(nproc_fluxes)*4
      character proc_rates_name(nproc_rates)*4
c
      data proc_rates_name  /'ex  ','diff', 'ent ', 'adv ', 'ncp ',
     $                       'sim ','obs '/
      data proc_fluxes_name /'ex  ','diff', 'ent ', 'adv ', 'ncp ',
     $                       'sim ','obs ', 'zcmp', 'tcal', 'tobs'/



      

