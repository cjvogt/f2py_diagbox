c =====================================================================
c
c     SUBROUTINE WRITE_SN_OUTPUT
c
c     PURPOSE: writes the output of the sensitivity runs to the unit
c                specified by input unit
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
c          p_init          : inital parameter
c          sg_p            : uncertainty of parameters
c
c          int_rates_std   : integrated rates of change of standard run
c          int_fluxes_std  : integrated fluxes of standard run
c          nyears          : number of analysis years
c          chisq_std       : chi**2 between sim and obs rates of std run
c          rsq_std         : r**2 between simulated and obs rates of std run
c          closure_std     : annual lack of closure of std run
c          costfn_std      : value of cost function of std run
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
c 
c =====================================================================
c
      subroutine write_sn_output(unit,
     $     notes,exp_name,
     $     tstart,tend,ts,     
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
     $     p_init,sg_p,
c
     $     int_rates_std,int_fluxes_std,nyears,
     $     d_int_rates_dp,d_int_fluxes_dp,
     $     chisq_std,rsq_std,closure_std,
     $     costfn_std,
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
# include "mc_statistics.h"
# include "sensitivity.h"

      logical error
      integer unit
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      integer i,np,j
      character date*17
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
     $            '---------',
     $            ' Results of sensitivity studies of the ',
     $            '    diagnostic box model',
     $            '    of the seasonal carbon cycle at ',
#ifdef BATS
     $            '    the U.S. JGOFS station BATS ',/
#elif BERM
     $            '    the Hydrostation 'S' ',/
#endif
     $            '    written by Nicolas Gruber, Aug. 1996 ',
     $            '    this is version 1.0 ',/
     $            '    based on the program model.F developed',
     $            '    during 1989 to 1992',
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
      write(unit,'(A,F5.1)')'   start time         : ',tstart
      write(unit,'(A,F5.1)')'   end time           : ',tend
      write(unit,'(A,F5.1)')'   time step (days)   : ',ts
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' specifications for location :'
#ifdef BATS
      write(unit,'(A)')     '   U.S. JGOFS station BATS '
#elif BERM
      write(unit,'(A)')     '   Hydrostation S '
#endif
c
c ---------------------------------------------------------------------
c     write specification for dc13
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)')' specifications for dC13 coefficients :'
c
      if (dc13_coeff .eq. 'o1') then
         write(unit,'(2A)')
     $        '     dc13 coefficients from original data with 1st ',
     $        'order fit'
      elseif (dc13_coeff .eq. 'o2') then
         write(unit,'(2A)')
     $        '     dc13 coefficients from original data with 2nd ',
     $        'order fit'
      elseif (dc13_coeff .eq. 'a2') then
         write(unit,'(2A)')
     $        '     dc13 coefficients from artificial data with 2nd ',
     $        'order fit'
      elseif (dc13_coeff .eq. 'n1') then
         write(unit,'(2A)')
     $        '     dc13 coefficients from adjusted data with 1st ',
     $        'order fit'
      elseif (dc13_coeff .eq. 'n2') then
         write(unit,'(2A)')
     $        '     dc13 coefficients from adjusted data with 2nd ',
     $        'order fit'
      elseif (dc13_coeff .eq. 'm2') then
         write(unit,'(2A)')
     $        '     dc13 coefficients from manipulated adj. data ',
     $        'with 2nd order fit'
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
         write(unit,'(A,F6.2)')
     $        '     horizontal velocity [m s-1]            : ',u
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
      write(unit,'(/A/)') ' investigated parameters : '
c
      write(unit,'(3(a9,1x)/)') 
     $     'parameter',
     $     'initial  ',
     $     '1-sg     '
       
      do np = 1,nparams
         write(unit,'(5x,i2,3x,g10.4,3x,g10.4)') 
     $        np,p_init(np),sg_p(np)
      enddo
c
c ---------------------------------------------------------------------
c     write statistics of standard run
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' run statistics of standard run : '
c
      write(unit,'(A,F8.4)') '    chi**2 (rates)         : ',
     $     chisq_std
      write(unit,'(A,F6.4)') '    r**2 (rates)           : ',
     $     rsq_std
      write(unit,'(A,F6.2)') '    annual lack of closure : ',
     $     closure_std
      write(unit,'(A,F10.6)') '    cost function          : ',
     $     costfn_std
c
c ---------------------------------------------------------------------
c     write integrated fluxes
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' Integrated fluxes [mol m-2]'
      write(unit,9100) (i, i=1,nyears)
c
      do j=1,nproc_fluxes
         write(unit,9200) proc_fluxes_name(j),
     $        (int_fluxes_std(i,j), i=1,nyears)
      enddo
c
c ---------------------------------------------------------------------
c     write integrated rates
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' Integrated rates [umol kg-1]'
      write(unit,9100) (i, i=1,nyears)
c
      do j=1,nproc_rates
         write(unit,9300) proc_rates_name(j),
     $        (int_rates_std(i,j), i=1,nyears)
      end do
c
c ---------------------------------------------------------------------
c     write sensitivities of integrated fluxes
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' parameter sensitivities for ',
     $     'integrated fluxes [mol m-2] d_int_fluxes_dp * sg_p : '
c
c ---------------------------------------------------------------------
c     loop over all processes
c ---------------------------------------------------------------------
c
      do j = 1,nproc_fluxes
c
c ---------------------------------------------------------------------
c     write results in tabular form
c ---------------------------------------------------------------------
c
         write(unit,'(//a,a/)') ' process : ',proc_fluxes_name(j)
c
         write(unit,9150) (i, i=1,nyears)
         write(unit,9420) (int_fluxes_std(i,j),i=1,nyears)
         do np = 1,nparams
            write(unit,9400) 
     $           np,p_init(np),sg_p(np),
     $           (d_int_fluxes_dp(np,i,j)*sg_p(np),i=1,nyears)
         enddo
c
c ---------------------------------------------------------------------
c     end of loop over all processes
c ---------------------------------------------------------------------
c
      enddo
c
c ---------------------------------------------------------------------
c     write sensitivities of integrated rates 
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' parameter sensitivities for ',
     $     'integrated rates [mol m-3] d_int_rates_dp * sg_p : '
c
c ---------------------------------------------------------------------
c     loop over all processes
c ---------------------------------------------------------------------
c
      do j = 1,nproc_rates
c
c ---------------------------------------------------------------------
c     write results in tabular form
c ---------------------------------------------------------------------
c
         write(unit,'(//a,a/)') ' process : ',proc_rates_name(j)
c
         write(unit,9150) (i, i=1,nyears)
         write(unit,9520) (int_rates_std(i,j),i=1,nyears)
         do np = 1,nparams
            write(unit,9500) 
     $           np,p_init(np),sg_p(np),
     $           (d_int_rates_dp(np,i,j)*sg_p(np),i=1,nyears)
         enddo
c
c ---------------------------------------------------------------------
c     end of loop over all processes
c ---------------------------------------------------------------------
c
      enddo
c
c ---------------------------------------------------------------------
c     draw a final line
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
c ---------------------------------------------------------------------
c     format statements
c ---------------------------------------------------------------------
c
 9800 format(/'-------------------------------------------------')
c
 9100 format(6x,5(1x,a10))
 9150 format(8x,'parameter',8x,5(1x,a10))
 9200 format(2x,a4,5(1x,F10.4))
 9300 format(2x,a4,5(3x,F8.2))
 9400 format(2x,i2,2(1x,g10.4),5(1x,F10.4))
 9420 format(/2x,'#    value     1-sg ',4x,5(1x,F10.4)/)
 9500 format(2x,i2,2(1x,g10.4),5(3x,F8.4))
 9520 format(/2x,'#    value     1-sg ',4x,5(3x,F8.4)/)
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
 9900 write(unit,*) '--WRITE_OUTPUT: an error occurred while writing..'
      error = .true.
c
      return

      end











