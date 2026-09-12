c     c--------------------------------------------------------------------------c
c     c-                                                                        -c
c     c- Author :  Roger Fink,                                                  -c
c     c-           Physics Institute,                                           -c
c     c-           Climate and Environmental Physics,                           -c
c     c-           CH-3012 Bern                                                 -c
c     c- Date   :  1996/07/16                                                   -c
c     c- Version:  2.00                                                         -c
c     c-                                                                        -c
c     c--------------------------------------------------------------------------c
c     c-                                                                        -c
c     c- procedure getco2                                                       -c
c     c- ****************                                                       -c
c     c-                                                                        -c
c     c- Calculation of the CO2 concentration (i.e. H2CO3+CO2)                  -c
c     c-                                                                        -c
c     c-                                                                        -c
c     c- Input:                                                                 -c
c     c- ------                                                                 -c
c     c-                                                                        -c
c     c-  Guess   : user's estimate of root                                     -c
c     c-  Tol     : tolerance in answer                                         -c
c     c-  MaxIter : maximum number of iterations                                -c
c     c-                                                                        -c
c     c-  Dissociation constants:                                               -c
c     c-  K1, K2  : first and second of carbonic acid                           -c
c     c-  Kb      : boric acid                                                  -c
c     c-  Ksi     : silicium acid                                               -c
c     c-  Kp2, Kp3: second and third of phosphoric acid                         -c
c     c-  Kw      : dissociation of water                                       -c
c     c-                                                                        -c
c     c-  Total concentrations:                                                 -c
c     c-  TCO2    : total dissolved inorganic carbon                            -c
c     c-  TB      : total borate                                                -c
c     c-  TSi     : total silicate                                              -c
c     c-  TP      : total phosphate                                             -c
c     c-                                                                        -c
c     c-  TAlk    : total alkalinity                                            -c
c     c-                                                                        -c
c     c-  UF      : unit factor (see comments below)                            -c
c     c-                                                                        -c
c     c-                                                                        -c
c     c- Output:                                                                -c
c     c- -------                                                                -c
c     c-                                                                        -c
c     c- Root     : calculated roots                                            -c
c     c- kH       : total hydrogen ion concentration (SWS)                      -c
c     c-                                                                        -c
c     c- Iter     : number of iterations it took to find root                   -c
c     c- Error    : flags if something went wrong                               -c
c     c-                                                                        -c
c     c- Errors: 0: No error                                                    -c
c     c-         1: Iter >= MaxIter                                             -c
c     c-         2: The slope was zero at some point                            -c
c     c-         3: Tol <= 0                                                    -c
c     c-         4: MaxIter < 0                                                 -c
c     c-                                                                        -c
c     c-                                                                        -c
c     c- All variables are double precision, except "Iter" and "Error", which   -c
c     c- are integer.                                                           -c
c     c-                                                                        -c
c     c-                                                                        -c
c     c- rem.:                                                                  -c
c     c- It is often useful (e.g. if this subroutine is used in a model) to pass-c
c     c- to "Guess" and to "Root" the same argument of the program part, which  -c
c     c- uses this subroutine. This can be done here without any problems. (In  -c
c     c- this case the value of this argument will be changed from the "Guess"  -c
c     c- to the "Root" value.)                                                  -c
c     c-                                                                        -c
c     c--------------------------------------------------------------------------c
c     c--------------------------------------------------------------------------c
c     c- remark to the units:                                                   -c
c     c- --------------------                                                   -c
c     c-                                                                        -c
c     c- The concentration units should be given in mol/kg-solution             -c
c     c- (eq/kg-solution for TA).                                               -c
c     c- If different units are used, they have to be converted to              -c
c     c- mol/kg-solution.                                                       -c
c     c- In both cases set UF=1.d0.                                             -c
c     c-                                                                        -c
c     c- Alternatively, the conversion to mol/kg-solution may be scipped, which -c
c     c- can save CPU time if e.g. CO2 have to be calculated at every time step -c
c     c- of an integration.                                                     -c
c     c- Only the contribution to the water alkalinity has to be corrected for  -c
c     c- the `wrong' units, all other `wrong' concentrations cancel each other. -c
c     c- This is performed by UF. E.g.:                                         -c
c     c-                                                                        -c
c     c- Input units:  mol/m3 ---> UF = density of water in units of kg/m3      -c
c     c-              umol/kg ---> UF = 1.d6                                    -c
c     c-                                                                        -c
c     c-                                                                        -c
c     c- Note that the total hydrogen ion concentration kH, which is related to -c
c     c- the SWS pH-scale, satisfies for numerical reasons the unit relation    -c
c     c- "[kH]=[dissociation constants]".                                       -c
c     c-                                                                        -c
c     c--------------------------------------------------------------------------c
c     c-                                                                        -c
c     c- The dissociation constants have to referre to the SWS pH scale         -c
c     c- (Alkalinity terms!).                                                   -c
c     c-                                                                        -c
c     c--------------------------------------------------------------------------c
c     c-                                                                        -c
c     c-  The code for the Newton-Raphson's algorithm, which have been changed  -c
c     c-  much for this purpose, has been originaly taken from                  -c
c     c-                                                                        -c
c     c-     Turbo Pascal Numerical Methods Toolbox                             -c
c     c-     (C) Copyright 1986 Borland International.                          -c
c     c-                                                                        -c
c     c-      Version Date (Borland): 26 January 1987                           -c
c     c-                                                                        -c
c     c--------------------------------------------------------------------------c
      subroutine getco2(
c
c     Input:
c     ------
c     
c          Startwert des CO2 und Abbruchkriterien:
     $     Guess,
     $     Tol,
     $     MaxIter,
c     
c          Dissoziationskonstanten:
     $     K1, K2, Kb, Ksi, Kp2, Kp3,
c     
c          Ionenprodukt des Wassers:
     $     Kw,
c     
c          Massenbilanz:
     $     TCO2, TB, TSi, TP,
c     
c          Ladungsbilanz:
     $     TAlk,
c     
c          Unit factor:
     $     UF,
c     
c     
c     Output:
c     -------
c     
c          solution of CO2 and total hydrogen ion concentration (SWS):
     $     Root,
     $     kH,
     $     Iter,
     $     Error
     $     )
c     
c     
c     
      implicit none
c     
c     global variables:
c     -----------------
      integer MaxIter, Iter, Error
      double precision
     $     Guess,
     $     Tol,
     $     K1, K2, Kb, Ksi, Kp2, Kp3,
     $     Kw,
     $     TCO2, TB, TSi, TP,
     $     TAlk,
     $     UF,
     $     Root,
     $     kH
c     
c     
c
c     local variables:
c     ----------------
c
c     constant to check if the derivative of the function is not to small
c     or if at the estimated root satisfies the root condition:
      double precision NearZero
      parameter(NearZero  = 1.0d-15)
c
c     Flags that a root has been found:
      logical found
c
c     Iteration variables:
      double precision
     $     co2Old, OldY, OldDeriv,
     $     co2New, NewY, NewDeriv
c
c
c     alkalinity:
      double precision Ac, Ab, Aw, Asi, Ap
c
      double precision X, h2, h3, K_, S_
c
c     Partielle Ableitungen:
      double precision AcX, AcCO2, AbkH, AsikH, ApkH, AwkH, kHX, XCO2,
     $     h2kH, h3kH
c
c
c
c
c
c     initials:
c     =======================================================================
c
      found = .false.
      Iter = 0
      Error = 0
c
      co2Old = Guess
c
      K_=dsqrt(K1/K2)
      S_=dsqrt(K1*K2)
c
c
c     calculating the value of the function (for co2=co2Old) and
c     it's derivative:
c     -------------------------------------------------------------
c
c     Rechnerische Hilfsgroessen:
      X=-K_/2.d0+dsqrt((K_*K_/4.d0)+TCO2/co2Old-1.d0)
      kH=S_/X
      h2=1.d0+kH/Kp2+Kp3/kH
      h3=1.d0+kH/Kp3+(kH*kH)/(Kp2*Kp3)
c
c     Berechnung der Alkalinitaeten:
      Ac=(K_+2.d0*X)*X*co2Old
      Ab=TB*Kb/(kH+Kb)
      Asi=TSi*Ksi/(kH+Ksi)
      Ap=TP*(1.d0/h2+2.d0/h3)
      Aw=(Kw/kH-kH)*UF
c
c      write(*,*) 'Alkalinities  before the iteration : '
c      write(*,'(a,x,f7.2)') 'Carbonate-Alk : ',Ac
c      write(*,'(a,x,f7.2)') 'Borate-Alk    : ',Ab
c      write(*,'(a,x,f7.2)') 'Silicate-Alk  : ',Asi
c      write(*,'(a,x,f7.2)') 'Phosphate-Alk : ',Ap
c      write(*,'(a,x,f7.2)') 'water-Alk     : ',Aw
c
c     Berechnung der partiellen Ableitungen:
      AcX=(K_+4.d0*X)*co2Old
      AcCO2=(K_+2.d0*X)*X
      XCO2=-TCO2/((2.d0*X+K_)*(co2Old*co2Old))
      AbkH=-Ab/(kH+Kb)
      AsikH=-Asi/(kH+Ksi)
      h2kH=1.d0/Kp2-Kp3/(kH*kH)
      h3kH=(1.d0+2.d0*kH/Kp2)/Kp3
      ApkH=-TP*(h2kH/(h2*h2)+2.d0*h3kH/(h3*h3))
      AwkH=(-Kw/(kH*kH)-1.d0)*UF
      kHX=-kH/X
c
      OldY = Ac+Ab+Asi+Ap+Aw-TAlk
      OldDeriv = AcCO2+(AcX+kHX*(AbkH+AsikH+ApkH+AwkH))*XCO2
c
c
c
c
      if (dabs(OldDeriv) .le. NearZero) Error = 2
      if (Tol .le. 0.d0) Error = 3
      if (MaxIter .le. 0) Error = 4
c
c
c
c
c     performe the iteration:
c     =======================
      do Iter=1,MaxIter
c
         co2New = co2Old - OldY / OldDeriv
c
c
c        calculating the value of the function (for co2=co2New) and
c        it's derivative:
c        -------------------------------------------------------------
c        Rechnerische Hilfsgroessen:
         X=-K_/2.d0+dsqrt((K_*K_/4.d0)+TCO2/co2New-1.d0)
         kH=S_/X
         h2=1.d0+kH/Kp2+Kp3/kH
         h3=1.d0+kH/Kp3+(kH*kH)/(Kp2*Kp3)
c     
c        Berechnung der Alkalinitaeten:
         Ac=(K_+2.d0*X)*X*co2New
         Ab=TB*Kb/(kH+Kb)
         Asi=TSi*Ksi/(kH+Ksi)
         Ap=TP*(1.d0/h2+2.d0/h3)
         Aw=(Kw/kH-kH)*UF
c     
c        Berechnung der partiellen Ableitungen:
         AcX=(K_+4.d0*X)*co2New
         AcCO2=(K_+2.d0*X)*X
         XCO2=-TCO2/((2.d0*X+K_)*(co2New*co2New))
         AbkH=-Ab/(kH+Kb)
         AsikH=-Asi/(kH+Ksi)
         h2kH=1.d0/Kp2-Kp3/(kH*kH)
         h3kH=(1.d0+2.d0*kH/Kp2)/Kp3
         ApkH=-TP*(h2kH/(h2*h2)+2.d0*h3kH/(h3*h3))
         AwkH=(-Kw/(kH*kH)-1.d0)*UF
         kHX=-kH/X
c     
         NewY = Ac+Ab+Asi+Ap+Aw-TAlk
         NewDeriv = AcCO2+(AcX+kHX*(AbkH+AsikH+ApkH+AwkH))*XCO2
c     
c
c        testing, if a root has been found:
c        ==================================
         found = (dabs(co2New - co2Old) .le. dabs(co2Old*Tol))
c
cc        infos for debugging:
c         write(6,*) 'Iter = ',Iter
c         write(6,*) 'OldY = ',OldY
c         write(6,*) 'NewY = ',NewY
c         write(6,*) 'co2Old = ',co2Old
c         write(6,*) 'co2New = ',co2New
c         write(6,*)
c         write(6,*)
c
         co2Old = co2New
         OldY = NewY
         OldDeriv = NewDeriv
c
         if (found .or. error .ne. 0 .or. Iter .ge. MaxIter) goto 5555
c
c        
c        testing, if the derivative is nearly zero:
c        ==========================================
         if ( .not. found) then
            if (dabs(OldDeriv) .le. NearZero) Error = 2
         endif
c
      enddo
 5555 continue
c
c
c
c     assignements for the output:
c     ============================
      Root = co2Old
      if (.not. found .and. Error .eq. 0 .and. Iter .gt. MaxIter) then
         Error = 1
      endif
c
c      write(*,*) 'Alkalinities  : '
c      write(*,'(a,x,f7.2)') 'Carbonate-Alk : ',Ac
c      write(*,'(a,x,f7.2)') 'Borate-Alk    : ',Ab
c      write(*,'(a,x,f7.2)') 'Silicate-Alk  : ',Asi
c      write(*,'(a,x,f7.2)') 'Phosphate-Alk : ',Ap
c      write(*,'(a,x,f7.2)') 'water-Alk     : ',Aw
      
      
c
c
      end 
c     end getco2
c
c
c
