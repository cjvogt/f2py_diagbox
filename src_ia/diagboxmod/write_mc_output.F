c =====================================================================
c
c     SUBROUTINE WRITE_MC_OUTPUT
c
c     PURPOSE: writes the parameters and the integrated values of
c               the standard run to file or screen, depending
c               on value of unit
c     
c
c     VARIABLES:
c 
c     IN:  unit            : unit number to write to
c          notes           : additional notes for run
c          tstart          : start time of integration
c          tend            : end time of integration
c          ts              : time step
c          nstep           : number of steps
c          sol_scheme      : solution scheme
c          dc13_coeff      : character specifying which dC13 coefficients
c                              are used       
c                            <o1> : original data with 1st order fit   
c                            <o2> : original data with 2nd order fit   
c                            <a2> : artficial data with 2nd order fit   
c                            <n1> : adjusted data with 1st order fit
c                            <n2> : adjusted data with 2nd order fit   
c                            <m2> : adj. manip. data with 2nd o. fit  
c          pv_relship      : character indicating the relationship
c                               that is used for calculating the
c                               piston-velocity
c                              <lm> : Liss and Merlivat
c                              <wa> : Wanninkhof
c          D_fco2_corr     : correction term for mixed layer fco2
c          gasex_fact      : air-sea exchange multiplication factor
c                              (for lm usually 1.7447)
c          ws_coeff        : character specifying wind speed coefficient
c                              <ih> : Isemer and Hasse, 1985
c                              <ae> : AEROCE
c          ent_scheme      : character indicating entrainment scheme
c                              to be used
c                              <in> : integrated vertical gradient scheme
c                              <mc> : mass conservation scheme
c                              <ep> : episodic event scheme
c          h_coeff         : character specifying which h coeff are taken
c          h_th            : depth at which thermocline values for 
c                                 mc entrainment are taken
c          c_th            : C concentration at h_th [umol kg-1]
c          dC13_th         : dC13 concentration at h_th [per mil]
c          dt_ent_re       : recurrence time for entrainment
c          const_lent      : logical indicating use of lent_const
c          lent_const      : entrainment length scale for matear ent. [m]
c          dlent_dh        : ratio of lent to mixed layer depth []
c          ent_fact        : entrainment multiplication factor
c          const_Kz        : logical indicating if constant Kz or
c                              seasonally varying Kz is used
c          Kz_const        : value of constant Kz 
c          const_vertgrad  : logical indicating if constant or variable
c                              vertical gradient is used
c          dCdz_const      : constant vertical C gradient 
c          ddC13_dC_bml    : ratio of ddC13/dz to dC/dz below mixed layer
c          diff_fact       : vertical diffusion mulitplication factor
c          incl_adv        : logical indicating if effect of advection
c                              is included
c          dC_dh           : horizontal gradient of C
c          ddC13_dh        : horizontal gradient of dC13
c          u               : horizontal velocity
c          D_dc13_org      : change in dC13 value of organic matter
c                          
c          nrun            : number of runs
c          p_init          : inital parameter
c          p_opt           : optimum parameters
c          sg_p            : uncertainty of parameters
c          int_rates_std   : integrated rates of change of standard run
c          int_fluxes_std  : integrated fluxes of standard run
c          nyears          : number of analysis years
c          sg_int_rates    : 1-sigma uncertainty of int_rates
c          wsg_int_rates   : weighted 1-sigma uncertainty of int_rates
c          sg_int_fluxes   : 1-sigma uncertainty of int_fluxes
c          wsg_int_fluxes  : weighted 1-sigma uncertainty of int_fluxes
c          chisq_std       : chi**2 between simulated and obs rates of std run
c          rsq_std         : r**2 between simulated and obs rates of std run
c          closure_std     : annual lack of closure of std run
c          costfn_std      : value of cost function of std run
c          mean_costfn     : mean cost function of all runs
c          costfn_min      : minimum of cost function
c          costfn_max      : maximum of cost function
c
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     19.08.96   ng     first implementation
c     18.10.99   ng     adapted for interannual variability
c     20.10.99   ng     added  some improvement
c 
c =====================================================================
c
      subroutine write_mc_output(unit,
     $     notes,exp_name,
     $     tstart,tend,ts,sol_scheme,
     $     dc13_coeff,
     $     pv_relship,D_fco2_corr,gasex_fact,ws_coeff,    
     $     ent_scheme,h_coeff,h_th,c_th,dc13_th,
     $     dt_ent_re,const_lent,lent_const,dlent_dh,ent_fact,
     $     const_Kz,Kz_const,
     $     const_vertgrad,dCdz_const,
     $     ddC13_dC_bml,diff_fact,
     $     incl_adv,dC_dh,ddC13_dh,u,
     $     D_dc13_org,
c
     $     nrun,
     $     p_init,p_opt,sg_p,
c
     $     int_rates_std,int_fluxes_std,nyears,
     $     sg_int_rates,sg_int_fluxes,
     $     wsg_int_rates,wsg_int_fluxes,
     $     chisq_std,rsq_std,closure_std,
     $     costfn_std,mean_costfn,costfn_min,costfn_max,
     $     error)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "periods.h"
# include "slabs.h"
c
# include "gasex_params.h"
# include "diffent_params.h"
# include "adv_params.h"
# include "ncp_params.h"
c
# include "intval.h"
# include "conj_grad.h"
# include "averages.h"
# include "monte_carlo.h"
# include "mc_statistics.h"

      logical error
      integer unit
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      integer i,ny,colons,i1,i2,ii,jmax
      integer n
      character date*17
      character*15 tmp(nyearmax)
      character*30 tmp2(nproc_fluxes)
      real avg_int_fluxes(nproc_fluxes),
     $     avg_sg_int_fluxes(nproc_fluxes)
      real avg_int_rates(nproc_rates),
     $     avg_sg_int_rates(nproc_rates)
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      error = .false.
c
c ---------------------------------------------------------------------
c     write general information 
c ---------------------------------------------------------------------
c
      if (unit .ne. 6) then
         write(unit,9000) 
 9000    format(//'------------------------------------------',
     $            '---------',/
     $            ' Results of monte carlo simulations of the ',/
     $            '    diagnostic box model',/
     $            '    of the seasonal carbon cycle at ',/
#ifdef BATS
     $            '    the U.S. JGOFS station BATS ',//
#elif BERM
     $            '    the Hydrostation S ',//
#else
     $            '    BATS/Station S combined ',//
#endif
     $            '    written by Nicolas Gruber, 1996-2001 ',/
     $            '    this is version 4.2 ',//
     $            '    based on the program model.F developed',/
     $            '    during 1989 to 1992',/
     $            '------------------------------------------',
     $            '--------'/)
c
      endif
c
c ---------------------------------------------------------------------
c     get date
c ---------------------------------------------------------------------
c
      call getdate(date)
c
c ---------------------------------------------------------------------
c     write general specification of this run
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A)')' general specifications for this model run :'
      write(unit,'(2A)')    '   date               : ',date
      write(unit,'(2A)')    '   name of experiment : ',exp_name 
      write(unit,'(A,F7.1)')'   start time         : ',tstart
      write(unit,'(A,F7.1)')'   end time           : ',tend
      write(unit,'(A,F7.5)')'   time step (years)  : ',ts
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' specifications for location :'
#ifdef BATS
      write(unit,'(A)')     '   U.S. JGOFS station BATS '
#elif BERM
      write(unit,'(A)')     '   Hydrostation S '
#else
      write(unit,'(A)')     '   BERM/BATS combined'
      
#endif
c
c ---------------------------------------------------------------------
c     write specification for solution scheme
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)')' specifications for solution scheme :'
c
      if (sol_scheme .eq. 'std') then
         write(unit,'(2A)') '   Standard (original) solution scheme ',
     $        'used, i.e. dC/dt is solved for'
      elseif (sol_scheme .eq. 'adv') then
         write(unit,'(2A)') '   Advection solution scheme used, ',
     $        'i.e. u is solved for'
      endif
c
c ---------------------------------------------------------------------
c     write specification for air-sea gasexchange
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)')' specifications for air-sea gasexchange :'
c
      if (pv_relship .eq. 'wa') then
         write(unit,'(2A)') '   Wanninkhof (1992) quadratic ',
     $        'piston velocity relationship used'
      elseif (pv_relship .eq. 'lm') then
         write(unit,'(2A)') '   Liss and Merlivat (1987) linear ',
     $        'piston velocity relationship used'
      endif
c         
      write(unit,'(A,F8.2)')'   correction term for fCO2 oc [uatm]) :', 
     $     D_fCO2_corr
      write(unit,'(2A,F7.4)')'   multiplication factor for air-sea ', 
     $     'gasexchange : ',gasex_fact
c
      if (ws_coeff .eq. 'ih') then
         write(unit,'(2A)')  '   wind-speed coefficients of Isemer',
     $        ' and Hasse'
      elseif (ws_coeff .eq. 'ae') then
         write(unit,'(A)')  '   wind-speed coefficients of AEROCE'
      endif
         
c
c ---------------------------------------------------------------------
c     write specification for vertical processes
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' specifications for vertical diffusion',
     $     ' and entrainment :'
c
c
      if (ent_scheme .eq.  'mc') then
c
         write(unit,'(A)') 
     $    '     mass conservation entrainment scheme '
         write(unit,'(A,F6.1)') 
     $        '     depth where thermocline values are taken : ',h_th
         write(unit,'(A,F6.1)')
     $        '     C concentration at this depth : ',c_th
         write(unit,'(A,F6.2)')
     $        '     dC13 value at this depth      : ',dc13_th
c
      elseif (ent_scheme .eq. 'ep') then
         write(unit,'(A)')
     $        '     episodic event entrainment scheme '
         write(unit,'(A,F6.1)')   
     $        '     entrainment recurrence period : ', dt_ent_re
c     
      elseif (ent_scheme .eq. 'in') then
         write(unit,'(A)')
     $        '     integrated vertical gradient scheme '
      elseif (ent_scheme .eq. 'ma') then
         write(unit,'(A)')
     $        '     Matear (1995) entrainment scheme '
         if (const_lent) then
            write(unit,'(A,F6.1)')
     $           '     constant entrainment length : ',lent_const
         else
            write(unit,'(A/,2A,F6.4)')
     $           '     variable entrainment length ',
     $           '     ratio of entrainment length to ',
     $           'mixed layer depth : ',dlent_dh
         endif
      endif
c
      if (h_coeff .eq. 's') then
         write(unit,'(A)') 
     $        '     h coefficients from station S (1981-89)'
      elseif (h_coeff .eq. 'b') then
         write(unit,'(A)')
     $        '     h coefficients from BATS (1991-94)'
      elseif (h_coeff .eq. 'c') then
         write(unit,'(A)')
     $        '     h coefficients from BATS CTD (1989-93)'
      elseif (h_coeff .eq. 'a') then
         write(unit,'(A)')
     $        '     h coefficients from avg BATS CTD (1989-93)'
      endif
c
      write(unit,'(A,F5.2)')
     $    '     multiplication factor for entrainment : ',ent_fact
c
      if (const_kz) then
         write(unit,'(A,e10.3)') 
     $    '     seasonally constant Kz : ',kz_const
      else
         write(unit,'(A)')
     $    '     seasonally varying Kz '
      endif
c
      if (ent_scheme .ne. 'mc') then
c
         if (const_vertgrad) then
            write(unit,'(A,f7.3)') 
     $           '     seasonally constant dC/dz [umol kg-1 m-1 : ',
     $           dCdz_const
         else
            write(unit,'(A)')
     $           '     seasonally varying dC/dz'
         endif
c
         write(unit,'(A,/A,f10.5)')
     $        '     ratio of vertical dC13 and C gradients ',
     $        '     below the mixed layer [per mil kg umol-1] : ',
     $        ddC13_dC_bml
c
      endif
c
      write(unit,'(A,F5.2)')
     $    '     multiplication factor for vertical diffusion : ',
     $     diff_fact
c
c ---------------------------------------------------------------------
c     write specification for advection
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' specifications for horizontal advection :'
c
      if (incl_adv) then
         if (sol_scheme .eq. 'std') then
            write(unit,'(A,F6.2)')
     $           '     horizontal velocity [m s-1]            : ',u
         endif
c
         write(unit,'(A,E10.4)')
     $        '     horizontal C gradient [umol kg-1 m-1]  : ',
     $        dC_dh
c
         write(unit,'(A,E10.4)')
     $        '     horizontal dC13 gradient [per mil m-1] : ',
     $        ddC13_dh
      else
         write(unit,'(A)')
     $        '     horizontal advection not included : '
      endif
c
c ---------------------------------------------------------------------
c     write specification for net community production
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)')' specification for ',
     $     'net community production :'
      write(unit,'(A,F6.2)')
     $    '     change in dC13 of organic matter : ',D_dc13_org
c
c ---------------------------------------------------------------------
c     write additional comments
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' Additional comments : ',notes
c
c ---------------------------------------------------------------------
c     write parameters
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' monte carlo parameters : '
c
      write(unit,'(4(a9,1x)/)') 
     $     'parameter',
     $     'initial  ',
     $     '1-sg     ',
     $     'optimum  '
       
      do n = 1,nparams
         write(unit,'(5x,i2,3x,f7.4,3x,f7.4,3x,f7.4)') 
     $        n,p_init(n),sg_p(n),p_opt(n)
      enddo
c
c ---------------------------------------------------------------------
c     write statistics of standard run
c ---------------------------------------------------------------------
c
      if (sol_scheme .eq. 'std') then
         write(unit,9800)
c     
c     
#ifdef reviewer_costfn
         write(unit,'(/A/)') ' run statistics (Reviewer costfn): '
         write(unit,'(A,F8.1)') '    chi**2 (values)        : ',
     $        chisq_std
         write(unit,'(A,F6.4)') '    r**2 (values)          : ',
     $        rsq_std
#else 
         write(unit,'(/A/)') ' run statistics : '
         write(unit,'(A,F8.4)') '    chi**2 (rates)         : ',
     $        chisq_std
         write(unit,'(A,F6.4)') '    r**2 (rates)           : ',
     $        rsq_std
#endif
         write(unit,'(A,F6.2)') '    annual lack of closure : ',
     $        closure_std
         write(unit,'(A,F10.8)') '    cost function (w/o par): ',
     $        costfn_std
      endif
c
c ---------------------------------------------------------------------
c     write statistics of monte carlo
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' monte carlo statistics : '
      write(unit,'(a,i4)')    '    number of mc runs      : ',
     $     nrun
      write(unit,'(A,f10.6)') '    mean cost function     : ',
     $     mean_costfn
      write(unit,'(A,f10.6)') '    minimal cost function  : ',
     $     costfn_min
      write(unit,'(A,f10.6)') '    maximal cost function  : ',
     $     costfn_max

c
c ---------------------------------------------------------------------
c     write integrated fluxes with unweighted uncertainties
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' Integrated fluxes with ',
     $     'uncertainties [mol m-2]'

      colons = 5
      jmax = (nyears-1)/colons + 1
      do ii=1,nyears,colons
         i1 = ii
         i2 = min((ii+colons-1),nyears)
c         
         write(unit,9100) (int(tstart)+ny-1, ny=i1,i2)


         do i=1,nproc_fluxes
            do ny = i1,i2
               write(tmp(ny),'(2x,F7.2,1x,F5.2)') 
     $              int_fluxes_std(ny,i),sg_int_fluxes(ny,i)
            enddo
            write(unit,9200) proc_fluxes_name(i),
     $        (tmp(ny), ny=i1,i2)
         enddo
c
         write(unit,'(//)')
      enddo
      
cc
cc ---------------------------------------------------------------------
cc     write integrated fluxes with weighted uncertainties
cc ---------------------------------------------------------------------
cc
c      write(unit,9800)
cc
c      write(unit,'(/2A/)') ' Integrated fluxes with weighted ',
c     $     'uncertainties [mol m-2]'
c
c      colons = 5
c      jmax = (nyears-1)/colons + 1
c      do ii=1,nyears,colons
c         i1 = ii
c         i2 = min((ii+colons-1),nyears)
cc         
c         write(unit,9100) (int(tstart)+ny-1, ny=i1,i2)
c
c
c         do i=1,nproc_fluxes
c            do ny = i1,i2
c               write(tmp(ny),'(2x,F7.2,1x,F5.2)') 
c     $              int_fluxes_std(ny,i),wsg_int_fluxes(ny,i)
c            enddo
c            write(unit,9200) proc_fluxes_name(i),
c     $        (tmp(ny), ny=i1,i2)
c         enddo
cc
c         write(unit,'(//)')
c      enddo
cc
c ---------------------------------------------------------------------
c     compute and write average fluxes and uncertainty
c ---------------------------------------------------------------------
c
      do i = 1,nproc_fluxes
         avg_int_fluxes(i) = 0.d0
         avg_sg_int_fluxes(i) = 0.d0
      enddo

      do ny = 1,nyears
         do i = 1,nproc_fluxes
            avg_int_fluxes(i) = avg_int_fluxes(i)+int_fluxes_std(ny,i)
            avg_sg_int_fluxes(i) = avg_sg_int_fluxes(i)
     $           +sg_int_fluxes(ny,i)
         enddo
      enddo
c
      do i = 1,nproc_fluxes
         avg_int_fluxes(i) = avg_int_fluxes(i)/nyears
         avg_sg_int_fluxes(i) = avg_sg_int_fluxes(i)/nyears
      enddo
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' Average Integrated fluxes with ',
     $     'uncertainties[mol m-2]'
c
      write(unit,'(A)') '         Average  1-sg Uncertainty'
      do i=1,nproc_fluxes
         write(unit,'(2x,a4,2x,f7.2,1x,f5.2)') 
     $        proc_fluxes_name(i),avg_int_fluxes(i),
     $         avg_sg_int_fluxes(i)
      enddo
c

c
c ---------------------------------------------------------------------
c     compute and write anomaly fluxes
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' Anomalies of Integrated Fluxes with ',
     $     'uncertainties [mol m-2]'

      colons = 5
      jmax = (nyears-1)/colons + 1
      do ii=1,nyears,colons
         i1 = ii
         i2 = min((ii+colons-1),nyears)
c         
         write(unit,9100) (int(tstart)+ny-1, ny=i1,i2)
         do i=1,nproc_fluxes
            do ny = i1,i2
               write(tmp(ny),'(2x,F7.2,1x,F5.2)') 
     $           (int_fluxes_std(ny,i)-avg_int_fluxes(i)),
     $           sg_int_fluxes(ny,i)
            enddo
            write(unit,9300) proc_fluxes_name(i),
     $        (tmp(ny), ny=i1,i2)
         enddo
         write(unit,'(//)')
      enddo
      
      write(unit,9800)
c
c ---------------------------------------------------------------------
c     write integrated rates with unweighted uncertainties
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' Integrated rates with ',
     $     'uncertainties [umol kg-1]'

      colons = 5
      jmax = (nyears-1)/colons + 1
      do ii=1,nyears,colons
         i1 = ii
         i2 = min((ii+colons-1),nyears)
c         
         write(unit,9100) (int(tstart)+ny-1, ny=i1,i2)


         do i=1,nproc_rates
            do ny = i1,i2
               write(tmp(ny),'(2x,F7.1,1x,F5.1)') 
     $              int_rates_std(ny,i),sg_int_rates(ny,i)
            enddo
            write(unit,9300) proc_fluxes_name(i),
     $        (tmp(ny), ny=i1,i2)
         enddo
c
         write(unit,'(//)')
      enddo
cc
cc ---------------------------------------------------------------------
cc     write integrated rates with weighted uncertainties
cc ---------------------------------------------------------------------
cc
c      write(unit,9800)
cc
c      write(unit,'(/2A/)') ' Integrated rates with weighted ',
c     $     'uncertainties [umol kg-1]'
c      colons = 5
c      jmax = (nyears-1)/colons + 1
c      do ii=1,nyears,colons
c         i1 = ii
c         i2 = min((ii+colons-1),nyears)
cc         
c         write(unit,9100) (int(tstart)+ny-1, ny=i1,i2)
c
c
c         do i=1,nproc_rates
c            do ny = i1,i2
c               write(tmp(ny),'(2x,F7.1,1x,F5.1)') 
c     $              int_rates_std(ny,i),wsg_int_rates(ny,i)
c            enddo
c            write(unit,9300) proc_rates_name(i),
c     $        (tmp(ny), ny=i1,i2)
c         enddo
cc
c         write(unit,'(//)')
c      enddo
c
c      write(unit,9800)
cc
c ---------------------------------------------------------------------
c     compute and write average rates and uncertainty
c ---------------------------------------------------------------------
c
      do i = 1,nproc_rates
         avg_int_rates(i) = 0.d0
         avg_sg_int_rates(i) = 0.d0
      enddo

      do ny = 1,nyears
         do i = 1,nproc_rates
            avg_int_rates(i) = avg_int_rates(i)+int_rates_std(ny,i)
            avg_sg_int_rates(i) = avg_sg_int_rates(i)
     $           +sg_int_rates(ny,i)
         enddo
      enddo
c
      do i = 1,nproc_rates
         avg_int_rates(i) = avg_int_rates(i)/nyears
         avg_sg_int_rates(i) = avg_sg_int_rates(i)/nyears
      enddo
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' Average Integrated rates with ',
     $     'uncertainties[mol m-2]'
c
      write(unit,'(A)') '         Average  1-sg Uncertainty'
      do i=1,nproc_rates
         write(unit,'(2x,a4,2x,f7.2,1x,f5.2)') 
     $        proc_rates_name(i),avg_int_rates(i),
     $         avg_sg_int_rates(i)
      enddo
c

c
c ---------------------------------------------------------------------
c     compute and write anomaly rates
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' Anomalies of Integrated Rates with ',
     $     'uncertainties [mol m-2]'

      colons = 5
      jmax = (nyears-1)/colons + 1
      do ii=1,nyears,colons
         i1 = ii
         i2 = min((ii+colons-1),nyears)
c         
         write(unit,9100) (int(tstart)+ny-1, ny=i1,i2)
         do i=1,nproc_rates
            do ny = i1,i2
               write(tmp(ny),'(2x,F7.2,1x,F5.2)') 
     $           (int_rates_std(ny,i)-avg_int_rates(i)),
     $           sg_int_rates(ny,i)
            enddo
            write(unit,9300) proc_rates_name(i),
     $        (tmp(ny), ny=i1,i2)
         enddo

         write(unit,'(//)')
      enddo
      
      write(unit,9800)
cc
c ---------------------------------------------------------------------
c     output to files for plotting
c       changed output to be centered at 0.5
c ---------------------------------------------------------------------
c
      open(30,file="int_rate_mc_anomalies.set")
      do ny = 1,nyears
         do i = 1,nproc_rates
            write(tmp2(i),'(2x,F8.3,1x,F8.3,1x,F8.3)') 
     $           int_rates_std(ny,i)-avg_int_rates(i),
     $           int_rates_std(ny,i)-avg_int_rates(i)
     $           +sg_int_rates(ny,i),
     $           int_rates_std(ny,i)-avg_int_rates(i)
     $           -sg_int_rates(ny,i)           
         enddo
         write(30,9500) 
     $        int(tstart)+ny-0.5,
     $        (tmp2(i),i=1,nproc_rates)
      enddo
      close(30)
c
      open(32,file="int_flux_mc_anomalies.set")
      do ny = 1,nyears
         do i = 1,nproc_fluxes
            write(tmp2(i),'(2x,F8.3,1x,F8.3,1x,F8.3)') 
     $           int_fluxes_std(ny,i)-avg_int_fluxes(i),
     $           int_fluxes_std(ny,i)-avg_int_fluxes(i)
     $           +sg_int_fluxes(ny,i),
     $           int_fluxes_std(ny,i)-avg_int_fluxes(i)
     $           -sg_int_fluxes(ny,i)         
         enddo
 
         write(32,9500) 
     $        int(tstart)+ny-0.5,
     $        (tmp2(i),i=1,nproc_fluxes)
      enddo
      close(32)
 9500 format(1x,f9.2,10(a30))
c
c ---------------------------------------------------------------------
c     format statements
c ---------------------------------------------------------------------
c
 9800 format(/'-------------------------------------------------')
c
 9100 format(2x,5(9x,i4,2x))
 9200 format(2x,a4,5(a15))
 9300 format(2x,a4,5(a15))
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c     
      return
c
c --------------------------------------------------------------------
c     error handling
c --------------------------------------------------------------------
c     
 9900 write(unit,*) '--WRITE_MC_OUTPUT: an error occurred ',
     $     'while writing..'
      error = .true.
c
      return

      end











