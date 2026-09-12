c =====================================================================
c
c     SUBROUTINE USERINPUT
c 
c     PURPOSE: gets user input, all interactive handling is done
c              in this routine
c
c     VARIABLES:
c 
c     IN : none
c     OUT: exp_name        : name of experiment
c          notes           : character permitting the addition of 
c                               text to the model output
c          sim_mode        : character indicating simulation mode
c          sol_scheme      : character indicating solution scheme
c          nrun            : numbers of run for monte carlo simulation and
c                                optimization runs
c          tstart          : start time of integration
c          tend            : end time of integration
c          ts              : time step
c          def_par         : logical indicating default set of params
c                               is used
c 
c          dc13_coeff      : character specifying which dC13 coefficients
c                              are used                         
c          pv_relship      : character indicating the relationship
c                               that is used for calculating the
c                               piston-velocity
c                              <lm> : Liss and Merlivat
c                              <wa> : Wanninkhof
c          D_fco2_corr     : correction term for mixed layer fco2
c          gasex_fact      : air-sea exchange multiplication factor
c                              (for lm usually 1.7447)
c          ws_coeff        : character specifying wind speed coefficient
c          comp_fco2_online: logical specifying on- or offline comp.
c                              of fCO2
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
c          const_lent      : logical indicating use of ent_length
c          lent_const      : const. entrainment length scale [m] 
c                               for <ma> scheme
c          dlent_dh        : ratio of lent to mixed layer depth []
c          ent_fact        : entrainment multiplication factor
c          const_Kz        : logical indicating if constant Kz or
c                              seasonally varying Kz is used
c          const_vertgrad  : logical indicating if constant or variable
c                              vertical gradient is used
c          dCdz_const      : constant vertical DIC gradient 
c          D_dcdz          : correction term for dcdz
c          dd13C_dC_bml    : ratio of dd13C/dz to dC/dz below mixed layer
c          diff_fact       : vertical diffusion mulitplication factor
c          incl_adv        : logical indicating if effect of advection
c                              is included
c          dC_dh           : horizontal gradient of C
c          ddC13_dh        : horizontal gradient of dC13
c          u               : horizontal velocity
c          D_dc13_org      : change in dC13 of organic matter
c                          
c          plot_out        : logical indicating that plot output files
c                               are generated
c            
c
c     REVISION:
c
c     date     author  remarks
c   
c     31.07.96   ng    first implementation
c     12.08.96   ng    included sim_mode
c     21.08.96   ng    included sn sim_mode
c     09.09.96   ng    included const_lent
c     24.09.96   ng    included ws_coeff
c
c =====================================================================
c
      subroutine userinput(exp_name,notes,
     $     sim_mode,sol_scheme,nrun,def_par,
     $     tstart,tend,ts,
     $     dc13_coeff,
     $     pv_relship,D_fco2_corr,gasex_fact,ws_coeff,
     $     comp_fco2_online,
     $     ent_scheme,h_coeff,
     $     h_th,c_th,dc13_th,dt_ent_re,
     $     const_lent,lent_const,dlent_dh,ent_fact,
     $     const_Kz,kz_const,
     $     const_vertgrad,dCdz_const,D_dcdz,
     $     ddC13_dC_bml,diff_fact,
     $     incl_adv,dC_dh,ddC13_dh,u,
     $     D_dc13_org,
     $     plot_out)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "diagboxmod.h"
# include "simulation.h"
# include "monte_carlo.h"
# include "gasex_params.h"
# include "diffent_params.h"
# include "adv_params.h"
# include "ncp_params.h"
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      character inp1*1
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
c     
c --------------------------------------------------------------------
c     get name of experiment and notes
c --------------------------------------------------------------------
c     
      write(*,8000)
 8000 format(/' Please give name of experiment ',/
     $     '    (is used also for',
     $     ' naming the output files) (max length = 10) : ',$)
      read (*,*) exp_name
c
      write(*,8020)
 8020 format(/' Please add notes to this experiment ',
     $     '(max length = 30) : ',$)
      read (*,*) notes
c     
c --------------------------------------------------------------------
c     get mode of experiment
c --------------------------------------------------------------------
c     
 100  write(*,8100)
 8100 format(/' Please specify simulation mode : ',/
     $        '      <st> : standard simulation (DEFAULT) ',/
     $        '      <mc> : Monte Carlo simulations ',/
     $        '      <op> : parameter optimization ',/
     $        '      <sn> : evaluation of parameter sensitivities',/
     $        ' please type character*2 : ',$)
      read(*,*) sim_mode
c
      if ((sim_mode .ne. 'st') .and. (sim_mode .ne. 'mc') .and.
     $     (sim_mode .ne. 'op') .and. (sim_mode .ne. 'sn')) then
         write(*,*) 'wrong input. please type again...'
         goto 100
      endif
c     
c --------------------------------------------------------------------
c     get solution scheme
c --------------------------------------------------------------------
c     
 110  write(*,8110)
 8110 format(/' Please specify solution scheme : ',/
     $        '      <std> : original standard, solve for dC/dt ',/
     $        '      <adv> : solve for u velocity ',/
     $        '      <phy> : solve for Kz +  u ',/
     $        '      <kz>  : solve for Kz ',/
     $        ' please type character*2 : ',$)
      read(*,*) sol_scheme
c
      if ((sol_scheme .ne. 'std') .and. (sol_scheme .ne. 'adv') .and.
     $     (sol_scheme .ne. 'phy') .and. (sol_scheme .ne. 'kz')) then
         write(*,*) 'wrong input. please type again...'
         goto 110
      endif
c     
c --------------------------------------------------------------------
c     get numbers of run if monte carlo or optimization
c --------------------------------------------------------------------
c     
      if ((sim_mode .eq. 'mc') .or. (sim_mode .eq. 'mi')) then
 120     write(*,8120)
 8120    format(' please specify numbers of run for Monte Carlo',
     $        ' analysis : ',$)
         read(*,*) nrun
c
         if ((nrun .lt. 1) .or. (nrun .gt. nrunmax)) then
            write(*,*) 'input out of boundaries. please type again...'
            goto 120
         endif
      elseif (sim_mode .eq. 'op') then
 130     write(*,8130)
 8130    format(' please specify numbers of run for random initial',
     $        ' parameter optimization : ',$)
         read(*,*) nrun
c
         if ((nrun .lt. 1) .or. (nrun .gt. nrunmax)) then
            write(*,*) 'input out of boundaries. please type again...'
            goto 130
         endif
      endif
c     
c --------------------------------------------------------------------
c     get choice of parameters
c --------------------------------------------------------------------
c     
      write(*,8200)
 8200 format(/' Please specify the choice of parameters : ',/
     $        '    <d> : default parameters (DEFAULT) ',/
     $        '    <i> : interactive choice of parameters ',/
     $        ' please type character*1 : ',$)
      read(*,*) inp1
c
      if (inp1 .eq. 'i') then
         def_par = .false.
      elseif (inp1 .ne. 'd') then
         write(*,*) 'Wrong input. it is assumed that you meant',
     $        ' default parameters'
         def_par = .true.
      else
         def_par = .true.
      endif
c     
c --------------------------------------------------------------------
c     the following sequence is only done for interactive 
c       parameter choice
c --------------------------------------------------------------------
cc
      if (.not. def_par) then
c     
c --------------------------------------------------------------------
c     get simulation parameters
c --------------------------------------------------------------------
c
         write(*,8250)
 8250    format(//' Parameters for simulation ')
c
         write(*,8260)
 8260    format(/' please specify the start time ',
     $        '(DEFAULT = 1983.0) : ',$)
         read(*,*) tstart
c
 270     write(*,8270) 
 8270    format(/' please specify the ending time ',
     $        '(DEFAULT = 1999.0) : ',$)
         read(*,*) tend
c
         if (tend .le. tstart) then
            write(*,*) 'Wrong input. end time must be greater than',
     $           ' start time'
            write(*,*) '      please type again...'
            goto 270
         endif
c
 280     write(*,8280)
 8280    format(/' please specify the time step ',
     $         '(DEFAULT = 1/400 yr = 0.0025) : ',$)
         read(*,*) ts
c
         if (dint((tend-tstart)/ts) .ne. (tend-tstart)/ts) then
            write(*,*) 'Wrong input. (end time - start time)/ts ',
     $           'must be an integer '
            write(*,*) '      please type again...'
            goto 280
         elseif (int((tend-tstart)/ts) .gt. nstepmax) then
            write(*,*) 'Wrong input. too small timestep '
            write(*,*) '      please type again...'
            goto 280            
         endif
c     
c --------------------------------------------------------------------
c     get air-sea exchange parameters
c --------------------------------------------------------------------
c
         write(*,8300)
 8300    format(//' Parameters for air-sea exchange ')
c     
 320     write(*,8320)
 8320    format(/' Please specify the piston-velocity',
     $        ' relationship with wind-speed : ',/
     $        '     <lm> : Liss and Merlivat linear rel.ship',/
     $        '     <wa> : Wanninkhof quadratic rel.ship',/
     $        ' please type character*2 : ',$)
c     
         read(*,*) pv_relship
c     
         if ((pv_relship .ne. 'lm') .and. (pv_relship .ne. 'wa')) then
            write(*,*) 'Wrong input. Please type again...'
            goto 320
         endif
c     
         write(*,8340)
 8340    format(/' Please specify correction to calculated',
     $        ' fCO2 [uatm] (DEFAULT = 0) : ',$)
         read(*,*) D_fco2_corr
c     
         write(*,8360)
 8360    format(/' Please specify the air-sea exchange multiplication ',
     $        'factor : ')
         if (pv_relship .eq. 'lm') then
            write(*,8362)
 8362       format('   for Liss-Merlivat (DEFAULT = 1.7447) : ',$)
         elseif (pv_relship .eq. 'wa') then
            write(*,8364)
 8364       format('   for Wanninkhof (DEFAULT = 1.0) : ',$)
         endif
         read(*,*) gasex_fact
c     
c --------------------------------------------------------------------
c     on/offline calc of fCO2
c --------------------------------------------------------------------
c
         write(*,8370)
 8370    format(/' Please specify how fCO2(oc) should be calculated',/
     $           '  <n> : online ',/
     $           '  <f> : offline',/
     $           ' please type character*1 : ',$)
         read(*,*) inp1
c
         if (inp1 .eq. 'n') then
            comp_fco2_online = .true.
         else
            comp_fco2_online = .false.
         endif
c     
c --------------------------------------------------------------------
c     get vertical diffusion and entrainment parameters
c --------------------------------------------------------------------
c
         write(*,8500)
 8500    format(//' Parameters for vertical diffusion and entrainment')
c     
 520     write(*,8520)
 8520    format(/' Please specify the entrainment scheme : ',/
     $        '    <in> : integrated vertical gradient scheme ',
     $        '(DEFAULT) ',/
     $        '    <mc> : mass conservation scheme ',/
     $        '    <ep> : episodic event scheme ',/
     $        '    <ma> : matear scheme',/
     $        ' please type character*2 : ',$)
c     
         read(*,*) ent_scheme
c     
         if ((ent_scheme .ne. 'in') .and. (ent_scheme .ne. 'mc') .and.
     $        (ent_scheme .ne. 'ep') .and. (ent_scheme .ne. 'ma')) then
            write(*,*) 'Wrong input. please type again...'
            goto 520
         endif
c     
         if (ent_scheme .eq. 'ep') then
 540        write(*,8540)
 8540       format(' please specify the recurrence period of ',
     $           'the entrainment process (DEFAULT = 8 days) : ',$)
c     
            read(*,*) dt_ent_re
c     
            if ((dt_ent_re .lt. 0) .or. (dt_ent_re .gt. 365.)) then
               write(*,*) 'input out of boundaries. ',
     $              'please type again...'
               goto 540
            endif
         else
            dt_ent_re = 0.d0
         endif
c
         if (ent_scheme .eq. 'ma') then
            write(*,8542)
 8542       format(' do you wish do have a constant entrainment',
     $           ' length scale (y/n) (DEFAULT = y) : ',$)
            read(*,*) inp1
            if (inp1 .eq. 'y') then
               const_lent = .true.
            else
               const_lent = .false.
            endif
c
            if (const_lent) then
               write(*,8545)
 8545          format(' please specify an entrainment length scale ',
     $              '(DEFAULT = 1m) : ',$)
               read(*,*) lent_const
            else
               write(*,8546)
 8546          format(' please specify ratio of entrainment length ',
     $              'to mixed layer depth ',
     $              '(DEFAULT = 0.1) : ',$)
               read(*,*) dlent_dh
            endif
         endif
c
         if (ent_scheme .eq. 'mc') then
 550        write(*,8550)
 8550       format(' please specify the thermocline values : ')
            write(*,8552)
 8552       format('     give the depth where thermocline values ',
     $           'are taken (DEFAULT = 200 m) : ',$)
            read(*,*) h_th
c
            write(*,8554)
 8554       format('     give C concentration at this depth ',
     $           '(DEFAULT = 2066 umol/kg) : ',$)
            read(*,*) c_th
c
            write(*,8556)
 8556       format('     give dC13 value at this depth ',
     $           '(DEFAULT = 1.28 per mil) : ',$)
            read(*,*) dc13_th
         endif
c
         write(*,8557)
 8557    format(/' please specify entrainment multiplication ',
     $        'factor : ',$)
         read(*,*) ent_fact
c     
         write(*,8560)
 8560    format(/' Please specify if you want to use seasonally ',
     $        'constant Kz (y) (DEFAULT) ,'/
     $        '      or seasonally varying (n) ',
     $        ' (type <y> or <n>) : ',$)
         read(*,*) inp1
c     
         if (inp1 .eq. 'n') then 
            const_Kz = .false.
         else
            const_Kz = .true.
         endif
c
         if (const_kz) then
            write(*,8565)
 8565       format(/' please specify value of kz ',
     $           '(DEFAULT = 2.d-5 m2 s-1) : ',$)
            read(*,*) kz_const
         endif
c
         if (ent_scheme .ne. 'mc') then
            write(*,8570)
 8570       format(/' please specify if you want to use seasonally ',
     $           'constant vertical gradients (y) ',/
     $           '      or seasonally varying (n) (DEFAULT) ',
     $           '(type <n> or <y>) : ',$)
            read(*,*) inp1
c
            if (inp1 .eq. 'y') then
               const_vertgrad = .true.
            else
               const_vertgrad = .false.
            endif
c
            if (const_vertgrad) then
               write(*,8575)
 8575          format(/' please specify value of const dC/dz ',
     $              '(DEFAULT = 0.45 umol kg-1 m-1) : ',$)
               read(*,*) dCdz_const
            else
c
               write(*,8577)
 8577          format(/' please specify correction term for ',
     $              'dC/dz (DEFAULT = 0.0 umol kg-1 m-1) : ',$)
               read(*,*) D_dcdz
            endif
c
 580        write(*,8580)
 8580       format(/' Please specify the ratio of the vertical ',
     $           ' dC13 and vertical C gradients ',/
     $           '     below the mixed layer (per mil kg umol-1) ',
     $           '(DEFAULT = -0.0048) : ',$)
c     
            read(*,*) ddC13_dC_bml
c     
            if ((ddC13_dC_bml .gt. 0.d0) .or. 
     $           (ddC13_dC_bml .lt. -1.0d0)) then
               write(*,*) 'input out of boundaries. ',
     $              'please type again...'
               goto 580
            endif
         endif
c     
         write(*,8600)
 8600    format(/' Please specify the vertical diffusion ',
     $        ' multiplication factor (DEFAULT = 1) : ',$)
c     
         read(*,*) diff_fact
c
c --------------------------------------------------------------------
c     get horizontal advection parameters
c --------------------------------------------------------------------
c    
         write(*,8700)
 8700    format(//' Parameters for horizontal advection')
c     
         if ((sol_scheme .ne. 'adv') .and. (sol_scheme .ne. 'phy')) then
            write(*,8720)
 8720       format(/' Please specify if you want to include ',
     $           'horizontal advection'/
     $           '        (type <y> or <n>) (DEFAULT = n) : ',$)
c     
            read(*,*) inp1
c     
            if (inp1 .eq. 'y') then
               incl_adv = .true.
            else
               incl_adv = .false.
            endif
         else
            incl_adv = .true.
         endif
c     
         if (incl_adv) then
            write(*,8740)
 8740       format(/' please specify the horizontal C gradient ',
     $           '(DEFAULT = 1.1e-5 umol kg-1 m-1) : ',$)
            read(*,*) dC_dh
c     
            write(*,8760)
 8760       format(/' please specify the horizontal dC13 gradient ',
     $           '(DEFAULT = -1.2d-7 per mil m-1) : ',$)
            read(*,*) ddc13_dh
c     
            if ((sol_scheme .ne. 'adv') .and. 
     $            (sol_scheme .ne. 'phy')) then
                write(*,8780)  
 8780           format(/' please specify the horizontal ',
     $           'advection velocity ',
     $           '(DEFAULT = 0.05 m s-1) : ',$)
                read(*,*) u
             endif
         endif
c
c --------------------------------------------------------------------
c     get net community production parameters
c --------------------------------------------------------------------
c   
         write(*,8800)
 8800    format(//' Parameters for net community production')
c     
         write(*,8820)
 8820    format(/' please specify correction term dC13 value of',
     $           ' organic matter ',/
     $           '          (DEFAULT = 0.0 per mil) : ',$)
c     
         read(*,*) D_dc13_org

c --------------------------------------------------------------------
c     end of section of interactive parameter input
c --------------------------------------------------------------------
c    
      endif
c
c --------------------------------------------------------------------
c     get parameters for plotting
c --------------------------------------------------------------------
c    
      write(*,9000)
 9000 format(//' Parameters for plotting')
c     
      write(*,9020)
 9020 format(/' do you like plot output (type <y> or <n>) : ',$)
      read(*,*) inp1
c     
      if (inp1 .eq. 'y') then
         plot_out = .true.
      else
         plot_out = .false.
      endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      write(*,9900) 
 9900 format(//'--INPUT: terminated '/)
      return
      end

