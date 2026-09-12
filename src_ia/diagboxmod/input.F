c =====================================================================
c
c     SUBROUTINE INPUT.F
c 
c     PURPOSE: this routine initializes the values needed for
c                doing the simulation
c
c     VARIABLES:
c 
c     IN : time         : actual time
c          ts           : time step
c          dc13_coeff      : character specifying which dC13 coefficients
c                              are used       
c                            <o1> : original data with 1st order fit   
c                            <o2> : original data with 2nd order fit   
c                            <a2> : artficial data with 2nd order fit   
c          h_coeff      : character specifying which coefficients for
c                            mixed layer depth should be taken
c                            <b> : bats    <s> : station S
c
c          ws_coeff     : character specifying wind speed coefficient
c                            <ih> : Isemer and Hasse
c                            <ae> : AEROCE
c          comp_fco2_online: logical specifying how fCO2 should be calc.
c
c     OUT: C_new        : mixed layer C concentration at time t+ts
c          C13_new      : mixed layer C13 concentration at time t+ts
c          h_new        : mixed layer depth at time t+ts
c
c          C            : mixed layer C at time t+0.5*ts
c          C13          : mixed layer C13 at time t+0.5*ts
c          h            : mixed layer depth at time t+0.5ts
c          temp         : mixed layer temperature at time t+0.5*ts          
c          sal          : mixed layer salinity at time t+0.5*ts    
c          ws           : wind-speed at time t+0.5*ts             
c          fCO2_o       : CO2 fugacity in the mixed layer at time t+0.5*ts
c          pco2_a       : pCO2 in the atmosphere at time t+0.5*ts          
c          dc13_a       : dC13 in the atmosphere at time t+0.5*ts           
c          dCdz_var     : vertical C gradient below the mixed layer 
c                            at time t+0.5*ts  
c          kz_var       : value of seasonally varying kz at time t+0.5*ts
c          salk         : mixed layer salk at t+0.5*ts
c          error        : logical indicating error
c 
c     REVISION:
c
c     date     author  remarks
c   
c     2.08.96   ng    first implementation
c     4.09.96   ng    included dc13_coeff
c     5.09.96   ng    included possibility to do backwards integration
c    24.09.96   ng    included ws_coeff
c    15.10.99   ng    adapted for interannual variability runs, time
c                       is now in years, input data includes spline
c                       component, added error, fCO2 is now calculated
c                       interactively from alk, dic, t, and s.
c    06.06.00   ng    added option for fCO2
c    07.06.00   ng    avoid mld less than 15m, added salk
c
c =====================================================================
c
      subroutine input(time,ts,dc13_coeff,h_coeff,ws_coeff,
     $     comp_fco2_online,
     $     C_new,C13_new,h_new,
     $     C,c13,h,
     $     temp,sal,ws,
     $     salk,fCO2_o,pco2_a,dc13_a,dCdz_var,kz_var,error)
c
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
# include "diagboxmod.h"
# include "simulation.h"
# include "processes.h"
# include "slabs.h"
# include "harm_coeff.h"
# include "input.h"
# include "diffent_params.h"
# include "gasex_params.h"
      logical error
c
c ---------------------------------------------------------------------
c     external functions
c ---------------------------------------------------------------------
c
      double precision harmonic,delta_to_c13,spline
      double precision k1_mb,k2_mb
      double precision kb_dicks, k1si_sni,Kw_m95drcph
      double precision tb_uppst
      double precision kp1_ym95, kp2_ym95, kp3_ym95
      double precision solco2kg,co2smpl
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      double precision dc13_tmp
      double precision std_sal
#ifdef BATS_89_93
      parameter (std_sal = 36.6d0)
#else
      parameter (std_sal = 35.0d0)
#endif
c
c ---------------------------------------------------------------------
c     fCO2 parameters
c ---------------------------------------------------------------------
c
      double precision tp,tsi
      double precision k1,k2,kb,kw
      double precision k1si,kp1,kp2,kp3,k2divk1
      double precision tb,kh
      double precision tco2,talk
      double precision co2_guess,co2_iter
      double precision UF
      parameter(UF = 1.d6)
      integer maxiter
      parameter(maxiter=30)
      double precision tol
      parameter(tol=1.0d-5)    
c  
      integer iter,int_error
c     
c =====================================================================
c     begin of executable code
c =====================================================================
c
c
c ---------------------------------------------------------------------
c     regular integration
c      *_new variables are at t+delta_t
c      other variables are at t+1/2*delta_t
c ---------------------------------------------------------------------
c

      C_new =      harmonic(time+ts,c_sdic   ,7)  
     $                     + spline(time+ts,'sdic ',error)
c

      C     =      harmonic(time+0.5d0*ts,c_sdic   ,7)
     $                     + spline(time+0.5d0*ts,'sdic ',error)
c
      dc13_tmp =   harmonic(time+ts,c_dc13,7)
     $                     + spline(time+ts,      'dc13 ',error)
      c13_new =    delta_to_c13(dc13_tmp,C_new)
c
      dc13_tmp =   harmonic(time+0.5d0*ts,c_dc13,7) 
     $                     + spline(time+0.5d0*ts,'dc13 ',error)
      c13      =   delta_to_c13(dc13_tmp,C)
c
      h_new =   harmonic(time+ts,c_h,7)
     $                     + spline(time+ts,      'h    ',error)
c
      if (h_new .lt. 15.) then
         h_new = 15.
      endif
c         
      h     =   harmonic(time+0.5d0*ts,c_h,7)
     $                     + spline(time+0.5d0*ts,'h    ',error)
c
       if (h .lt. 15.) then
         h = 15.
      endif
     
c      
      salk =       harmonic(time+0.5d0*ts,c_salk  ,7) 
     $                     + spline(time+0.5d0*ts,'salk ',error)
      temp =       harmonic(time+0.5d0*ts,c_temp  ,7)
     $                     + spline(time+0.5d0*ts,'temp ',error)
      sal  =       harmonic(time+0.5d0*ts,c_sal   ,7)
     $                     + spline(time+0.5d0*ts,'sal  ',error)
      ws   =       harmonic(time+0.5d0*ts,c_ws  ,7)
     $                     + spline(time+0.5d0*ts,'ws   ',error)
      pco2_a =     harmonic(time+0.5d0*ts,c_pco2a,7)      
     $                     + spline(time+0.5d0*ts,'pco2a',error)
      dc13_a =     harmonic(time+0.5d0*ts,c_dc13a,7)
     $                     + spline(time+0.5d0*ts,'dc13a',error)
      kz_var =     harmonic(time+0.5d0*ts,c_kz_var,7)
     $                     + spline(time+0.5d0*ts,'kz   ',error)
      dCdz_var =   harmonic(time+0.5d0*ts,c_dCdz,  5)
     $                     + spline(time+0.5d0*ts,'dcdz ',error)
#ifdef debugging
c
c --------------------------------------------------------------------
c     for debugging
c --------------------------------------------------------------------
 
c      write(*,*)
c      write(*,*) '--INPUT: time = ',time
c      write(*,*) '  temp (harm + spline) : ',temp,
c     $     harmonic(time+0.5d0*ts,c_temp  ,5),
c     $     spline(time+0.5d0*ts,'temp ',error)
c      write(*,*) '  sal (harm + spline) : ',sal,
c     $     harmonic(time+0.5d0*ts,c_sal  ,5),
c     $     spline(time+0.5d0*ts,'sal  ',error)
c      write(*,*) '  ws (harm + spline) : ',ws, 
c     $     harmonic(time+0.5d0*ts,c_ws  ,5),
c     $     spline(time+0.5d0*ts,'ws   ',error)
c      write(*,*) '  h  (harm + spline) : ',h,
c     $     harmonic(time+0.5d0*ts,c_h,7),
c     $     spline(time+0.5d0*ts,'h    ',error)
c      write(*,*) '  h_new : ',h_new

#endif
c
c --------------------------------------------------------------------
c     compute fCO2, add salinity back to dic and alk
c --------------------------------------------------------------------
c   
      if (comp_fco2_online) then
c
c --------------------------------------------------------------------
c     remove salinity normalization
c --------------------------------------------------------------------
c   
         tco2 = C*sal/std_sal
         talk = salk*sal/std_sal    
c     
c ---------------------------------------------------------------------
c     determine the dissociation constants, use Mehrbach refitted
c ---------------------------------------------------------------------
c
         k1 = k1_mb(temp,sal)
         k2 = k2_mb(temp,sal)
         
         kb = kb_dicks(temp,sal)
         kw = Kw_m95drcph(temp,sal)
         k1si = k1si_sni()
         kp1 = kp1_ym95(temp,sal)
         kp2 = kp2_ym95(temp,sal)
         kp3 = kp3_ym95(temp,sal)
         
         k2divk1 = k2/k1
c     
c ---------------------------------------------------------------------
c     initialize total borate concentration
c       multiply with 1.d6 in order to get mumol/kg from mol/kg
c     initialize other nutrient conc to zero
c ---------------------------------------------------------------------
c
         tb = tb_uppst(sal)*1.d6
         
         tsi = 0.d0
         tp = 0.d0
c     
c ---------------------------------------------------------------------
c     calculate first_guess co2 concentration
c ---------------------------------------------------------------------
c
         co2_guess = co2smpl(talk,tco2,k2divk1,0.95d0)
c
c ---------------------------------------------------------------------
c     calculate co2 concentration iteratively using roger finks
c       routine
c ---------------------------------------------------------------------
c
         call getco2(co2_guess,tol,maxiter,
     $        k1,k2,kb,
     $        k1si,kp2,kp3,kw,
     $        tco2,tb,tsi,tp,
     $        talk,UF,
     $        co2_iter,kh,iter,int_error)
c     
c ---------------------------------------------------------------------
c     calculate fco2 in the ocean
c ---------------------------------------------------------------------
c
         fco2_o = co2_iter / solco2kg(temp,sal)
      
#ifdef debugging
c
c --------------------------------------------------------------------
c     for debugging
c --------------------------------------------------------------------
c 
c      write(*,*) 
c      write(*,*) '--INPUT: time = ',time
c      write(*,*) '  temp, sal  : ', temp,sal
c      write(*,*) '  sdic, salk : ', c,salk
c      write(*,*) '  dic, alk   : ', tco2,talk
c      write(*,*) '  co2_iter   : ', co2_iter
c      write(*,*) '  sol co2    : ', solco2kg(temp,sal)
c
#endif
c
c --------------------------------------------------------------------
c     compute fco2_o from coefficients
c --------------------------------------------------------------------
c        
      else
        fco2_o =     harmonic(time+0.5d0*ts,c_fco2_o,7)
     $        + spline(time+0.5d0*ts,'fco2o',error)
      endif
c
c --------------------------------------------------------------------
c     return to main program
c --------------------------------------------------------------------
c        
      return
      end

