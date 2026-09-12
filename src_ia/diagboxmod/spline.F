c =====================================================================
c
c     DOUBLE PRECISION FUNCTION SPLINE
c 
c     PURPOSE: this function interpolates the output of a spline function
c            stored in a file to a given time. The file is chosen
c            according to the variable given as an argument
c

c     VARIABLES:
c 
c     IN : time         : time (years)
c          var          : character identifying variable
c     OUT: spline       : interpolated value
c          error        : logical value indicating problem
c 
c     REVISION:
c
c     date     author  remarks
c   
c     15.10.99   ng    first implementation
c     06.06.00   ng    added new spline files (combined BERM/BATS)
c
c =====================================================================
c
      double precision function spline(time,var,error)
      implicit none
c
c ---------------------------------------------------------------------
c     global variables
c ---------------------------------------------------------------------
c
      double precision time
      character*5 var
      logical error
c
c ---------------------------------------------------------------------
c     functions
c ---------------------------------------------------------------------
c
      double precision iplin
c
c ---------------------------------------------------------------------
c     local variables
c ---------------------------------------------------------------------
c
      logical first,eof,ex
      data first /.true./
      integer ndatamax,ndata
      parameter (ndatamax = 10000)
      double precision dic(ndatamax),dc13(ndatamax)
      double precision alk(ndatamax)
      double precision ws(ndatamax),h(ndatamax)
      double precision temp(ndatamax),sal(ndatamax)
      double precision fco2o(ndatamax)
      double precision dc13a(ndatamax),pco2a(ndatamax)
      double precision kz(ndatamax),dcdz(ndatamax)
c
      double precision t_dic(ndatamax),t_dc13(ndatamax)
      double precision t_alk(ndatamax)
      double precision t_ws(ndatamax),t_h(ndatamax)
      double precision t_temp(ndatamax),t_sal(ndatamax)
      double precision t_fco2o(ndatamax)
      double precision t_dc13a(ndatamax),t_pco2a(ndatamax)
      double precision t_kz(ndatamax),t_dcdz(ndatamax)        
c
      character*45 sdic_spline_file,dc13_spline_file
      character*45 salk_spline_file
      character*45 ws_spline_file,h_spline_file
      character*45 temp_spline_file,sal_spline_file
      character*45 fco2o_spline_file
      character*45 dc13a_spline_file,pco2a_spline_file
      character*45 kz_spline_file,dcdz_spline_file

c     print chosen spline files to terminal
c
#ifdef const_splines
      data sdic_spline_file /'splines/sdic_spl_stas_83-99.dat'/
      data salk_spline_file /'splines/salk_spl_stas_83-99.dat'/
      data dc13_spline_file /'splines/dc13_spl_stas_83-99.dat'/
      data ws_spline_file /'splines/ws_spl_stas_83-99.dat'/
      data h_spline_file /'splines/h_spl_stas_83-99.dat'/
      data temp_spline_file /'splines/temp_spl_stas_83-99.dat'/
      data sal_spline_file /'splines/sal_spl_stas_83-99.dat'/
      data fco2o_spline_file /'splines/fco2o_spl_stas_83-99.dat'/
      data dc13a_spline_file /'splines/dc13a_spl_stas_83-99.dat'/
      data pco2a_spline_file /'splines/pco2a_spl_stas_83-99.dat'/
      data kz_spline_file /'splines/kz_spl_stas_83-99.dat'/
      data dcdz_spline_file /'splines/dcdz_spl_stas_83-99.dat'/

#elif BERM_BATS_83_98
c     used up to August '01
      data sdic_spline_file /'splines/spline_data.sdic_berm_bats'/
      data salk_spline_file /'splines/spline_data.salk_berm_bats'/
      data dc13_spline_file /'splines/spline_data.dc13_berm_bats'/
      data temp_spline_file /'splines/spline_data.temp_berm_bats'/
      data sal_spline_file  /'splines/spline_data.sal_berm_bats'/
      data fco2o_spline_file /'splines/spline_data.pco2_berm_bats'/
      data ws_spline_file   /'splines/ws_85-96_ncep_res_spline.out'/
      data h_spline_file /'splines/mld_ts_stas_81-98_res_spline.out'/
      data dc13a_spline_file /'splines/dc13a_spl_stas_83-99.dat'/
      data pco2a_spline_file /'splines/pco2a_spl_stas_83-99.dat'/
      data kz_spline_file /'splines/kz_spl_stas_83-99.dat'/
      data dcdz_spline_file /'splines/dcdz_spl_stas_83-99.dat'/
#else
c
c     new spline fits, updated May 3, 2002
c
      data sdic_spline_file 
     $     /'splines/spline_data.salk_berm_bats_extended'/
      data salk_spline_file 
     $     /'splines/spline_data.sdic_berm_bats_extended'/
      data dc13_spline_file 
     $     /'splines/spline_data.dc13_berm_bats_extended'/
      data temp_spline_file 
     $     /'splines/spline_data.temp_berm_bats_extended'/
      data sal_spline_file  
     $     /'splines/spline_data.sal_berm_bats_extended'/
      data fco2o_spline_file 
     $     /'splines/spline_data.pco2_berm_bats_020501'/
      data dc13a_spline_file /'splines/spline_atmdc13_ljo_bme.out'/
      data pco2a_spline_file /'splines/spline_atmco2_ljo_bme.out'/
      data dcdz_spline_file 
     $     /'splines/ddicdz_bats_1-156_res_spline.out'/
      data h_spline_file 
     $     /'splines/mld_ts_stas_bats_1-156_res_spline.out'/
c
c     not updated on May 3, 2002
c
      data ws_spline_file   
     $     /'splines/ws_ncep_monthly_spline.out'/

      data kz_spline_file /'splines/kz_spl_stas_83-99.dat'/
#endif      
      
      integer ndic,nalk,ndc13,nws,nh,ntemp,nsal
      integer ndc13a,npco2a,nkz,ndcdz,i,nfco2o
c
      integer ierr
c
      save t_h,h,nh
      save dic, t_dic,ndic, dc13,t_dc13, ndc13
      save alk,t_alk,nalk, ws,t_ws,nws
      save temp,t_temp,ntemp, sal,t_sal,nsal
      save fco2o,t_fco2o,nfco2o
      save dc13a, t_dc13a, ndc13a, pco2a,t_pco2a,npco2a
      save kz,t_kz,nkz, dcdz,t_dcdz,ndcdz
c
c
c =====================================================================
c     begin of executable code
c =====================================================================
c
      if (first) then
         
         eof = .false.
         write(*,*) '--SPLINE: getting splinefit data from files...'
        write(*,*) '--SPLINE: chosen spline files:'
        write(*,*) '  sdic  : ', trim(sdic_spline_file)
        write(*,*) '  salk  : ', trim(salk_spline_file)
        write(*,*) '  dc13  : ', trim(dc13_spline_file)
        write(*,*) '  temp  : ', trim(temp_spline_file)
        write(*,*) '  sal   : ', trim(sal_spline_file)
        write(*,*) '  fco2o : ', trim(fco2o_spline_file)
        write(*,*) '  ws    : ', trim(ws_spline_file)
        write(*,*) '  h     : ', trim(h_spline_file)
        write(*,*) '  dc13a : ', trim(dc13a_spline_file)
        write(*,*) '  pco2a : ', trim(pco2a_spline_file)
        write(*,*) '  kz    : ', trim(kz_spline_file)
        write(*,*) '  dcdz  : ', trim(dcdz_spline_file)
c         
c --------------------------------------------------------------------
c     check for all files
c --------------------------------------------------------------------
c             
         inquire(file=sdic_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) sdic_spline_file
            goto 9900
         endif
         inquire(file=dc13_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) dc13_spline_file
            goto 9900
         endif         
         inquire(file=salk_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) salk_spline_file
            goto 9900
         endif         
         inquire(file=temp_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) temp_spline_file
            goto 9900
         endif         
         inquire(file=sal_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) sal_spline_file
            goto 9900
         endif         
         inquire(file=fco2o_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) fco2o_spline_file
            goto 9900
         endif         
         inquire(file=ws_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) ws_spline_file
            goto 9900
         endif         
         inquire(file=h_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) h_spline_file
            goto 9900
         endif         
         inquire(file=dc13a_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) dc13a_spline_file
            goto 9900
         endif         
         inquire(file=kz_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) kz_spline_file
            goto 9900
         endif   
         inquire(file=dcdz_spline_file,exist = ex)
         if (.not. ex) then
            write(*,8000) dcdz_spline_file
            goto 9900
         endif
c
8000     format('--SPLINE: file ',a30,' does not exist.'/
     $             ' program stops.')
c     
c --------------------------------------------------------------------
c     get spline coefficients from files
c        start with sdic
c --------------------------------------------------------------------
c            
         open(20,file=sdic_spline_file,form="formatted",
     $         access="sequential")
         ndic = 0
         do while (.not. eof)
            ndic = ndic + 1
            read(20,*,end=1000) t_dic(ndic),dic(ndic)
         enddo
 1000    ndic = ndic - 1
         close(20)
c
c --------------------------------------------------------------------
c     get spline coefficients salk
c --------------------------------------------------------------------
c    
         open(20,file=salk_spline_file,form="formatted",
     $         access="sequential")
         nalk = 0
         do while (.not. eof)
            nalk = nalk + 1
            read(20,*,end=1010) t_alk(nalk),alk(nalk)
         enddo
 1010    nalk = nalk - 1         
         close(20)
c
c --------------------------------------------------------------------
c     get spline coefficients dc13
c --------------------------------------------------------------------
c    
         open(20,file= dc13_spline_file,form="formatted",
     $         access="sequential")
         ndc13 = 0
         do while (.not. eof)
            ndc13 = ndc13 + 1
            read(20,*,end=1020) t_dc13(ndc13),dc13(ndc13)
         enddo
 1020    ndc13 = ndc13 - 1         
         close(20)
c     
c --------------------------------------------------------------------
c     get spline coefficients temp
c --------------------------------------------------------------------
c    
         open(20,file=temp_spline_file,form="formatted",
     $         access="sequential")
         ntemp = 0
         do while (.not. eof)
            ntemp = ntemp + 1
            read(20,*,end=1030) t_temp(ntemp),temp(ntemp)
         enddo
 1030    ntemp = ntemp - 1         
         close(20)
c     
c --------------------------------------------------------------------
c     get spline coefficients salt
c --------------------------------------------------------------------
c    
         open(20,file=sal_spline_file,form="formatted",
     $         access="sequential")
         nsal = 0
         do while (.not. eof)
            nsal = nsal + 1
            read(20,*,end=1040) t_sal(nsal),sal(nsal)
         enddo
 1040    nsal = nsal - 1         
         close(20)
c     
c --------------------------------------------------------------------
c     get spline coefficients ws
c --------------------------------------------------------------------
c    
         open(20,file=ws_spline_file,form="formatted",
     $         access="sequential")
         nws = 0
         do while (.not. eof)
            nws = nws + 1
            read(20,*,end=1050) t_ws(nws),ws(nws)
         enddo
 1050    nws = nws - 1                  
         close(20)
c     
c --------------------------------------------------------------------
c     get spline coefficients fco2oc
c --------------------------------------------------------------------
c    
         open(20,file=fco2o_spline_file,form="formatted",
     $         access="sequential")
         nfco2o = 0
         do while (.not. eof)
            nfco2o = nfco2o + 1
            read(20,*,end=1055) t_fco2o(nfco2o),fco2o(nfco2o)
         enddo
 1055    nfco2o = nfco2o - 1                  
         close(20)
c
c --------------------------------------------------------------------
c     get spline coefficients mixed layer depth h
c --------------------------------------------------------------------
c    
         open(20,file=h_spline_file,form="formatted",
     $         access="sequential")
         nh = 0
         do while (.not. eof)
            nh = nh + 1
            read(20,*,end=1060) t_h(nh),h(nh)
         enddo
 1060    nh = nh - 1         
         close(20)
c     
c --------------------------------------------------------------------
c     get spline coefficients pco2a
c --------------------------------------------------------------------
c    
         open(20,file=pco2a_spline_file,form="formatted",
     $         access="sequential")
         npco2a = 0
         do while (.not. eof)
            npco2a = npco2a + 1
            read(20,*,end=1070) t_pco2a(npco2a),pco2a(npco2a)
         enddo
 1070    npco2a = npco2a - 1         
         close(20)
c     
c --------------------------------------------------------------------
c     get spline coefficients dc13a
c --------------------------------------------------------------------
c    
         open(20,file=dc13a_spline_file,form="formatted",
     $         access="sequential")
         ndc13a = 0
         do while (.not. eof)
            ndc13a = ndc13a + 1
            read(20,*,end=1080) t_dc13a(ndc13a),dc13a(ndc13a)
         enddo
 1080    ndc13a = ndc13a - 1         
         close(20)
c     
c --------------------------------------------------------------------
c     get spline coefficients kz
c --------------------------------------------------------------------
c    
         open(20,file=kz_spline_file,form="formatted",
     $         access="sequential")
         nkz = 0
         do while (.not. eof)
            nkz = nkz + 1
            read(20,*,end=1090) t_kz(nkz),kz(nkz)
         enddo
 1090    nkz = nkz - 1        
         close(20)
c     
c --------------------------------------------------------------------
c     get spline coefficients dcdz
c --------------------------------------------------------------------
c    
         open(20,file=dcdz_spline_file,form="formatted",
     $         access="sequential")
         ndcdz = 0
         do while (.not. eof)
            ndcdz = ndcdz + 1
            read(20,*,end=1100) t_dcdz(ndcdz),dcdz(ndcdz)
         enddo
 1100    ndcdz = ndcdz - 1         
         close(20)
c     
c --------------------------------------------------------------------
c    give confirmation
c --------------------------------------------------------------------
c             
         write(*,*) '--SPLINE: finished getting data from file.'
         write(*,*) '--SPLINE: number of observations for'
         write(*,*) '          sdic  : ',ndic
         write(*,*) '          salk  : ',nalk
         write(*,*) '          dc13  : ',ndc13
         write(*,*) '          temp  : ',ntemp
         write(*,*) '          salt  : ',nsal
         write(*,*) '          fco2o : ',nfco2o
         write(*,*) '          ws    : ',nws
         write(*,*) '          mld   : ',nh
         write(*,*) '          pco2a : ',npco2a
         write(*,*) '          dc13a : ',ndc13a
         write(*,*) '          kz    : ',nkz
         write(*,*) '          dc/dz : ',ndcdz
c     
c --------------------------------------------------------------------
c    start and end time
c --------------------------------------------------------------------
c             
         write(*,*) '--SPLINE: start and end time'
         write(*,*) '          sdic  : ',t_dic(1),t_dic(ndic),ndic
         write(*,*) '          salk  : ',t_alk(1),t_alk(nalk),nalk
         write(*,*) '          dc13  : ',t_dc13(1),t_dc13(ndc13),
     $        ndc13
         write(*,*) '          temp  : ',t_temp(1),t_temp(ntemp),
     $        ntemp
         write(*,*) '          salt  : ',t_sal(1),t_sal(nsal),nsal
         write(*,*) '          fco2o : ',t_fco2o(1),t_fco2o(nfco2o),
     $        nfco2o
         write(*,*) '          ws    : ',t_ws(1),t_ws(nws),nws
         write(*,*) '          mld   : ',t_h(1),t_h(nh),nh
         write(*,*) '          pco2a : ',t_pco2a(1),t_pco2a(npco2a),
     $        npco2a
         write(*,*) '          dc13a : ',t_dc13a(1),t_dc13a(ndc13a),
     $        ndc13a
         write(*,*) '          kz    : ',t_kz(1),t_kz(nkz),nkz
         write(*,*) '          dc/dz : ',t_dcdz(1),t_dcdz(ndcdz),
     $        ndcdz
c     
c --------------------------------------------------------------------
c     end of pass for first call
c --------------------------------------------------------------------
c     
         first = .false.
      endif
c     
c --------------------------------------------------------------------
c     give info
c --------------------------------------------------------------------
c   
#ifdef debugging  
c      write(*,*) '--SPLINE: value of var : ',var
c      write(*,*) '          value of time: ',time
c     
#endif
c
c --------------------------------------------------------------------
c     determine value at time t, use linear interpolation in
c         time
c --------------------------------------------------------------------
c     
      if (var .eq. 'sdic ') then
c         write(*,*) '--> sdic'
         spline = iplin(time,t_dic,dic,ndic,ierr)
         if (ierr .ne. 0) goto 9900
      elseif (var .eq. 'salk ') then
c         write(*,*) '--> salk'
         spline = iplin(time,t_alk,alk,nalk,ierr)
         if (ierr .ne. 0) goto 9900
      elseif (var .eq. 'dc13 ') then
c         write(*,*) '--> dc13'
         spline = iplin(time,t_dc13,dc13,ndc13,ierr)
         if (ierr .ne. 0) goto 9900
       elseif (var .eq. 'temp ') then
c         write(*,*) '--> temp'
         spline = iplin(time,t_temp,temp,ntemp,ierr)
         if (ierr .ne. 0) goto 9900
       elseif (var .eq. 'sal  ') then
c         write(*,*) '--> sal'
         spline = iplin(time,t_sal,sal,nsal,ierr)
         if (ierr .ne. 0) goto 9900
       elseif (var .eq. 'fco2o') then
c          write(*,*) '--> fco2o'
        spline = iplin(time,t_fco2o,fco2o,nfco2o,ierr)
         if (ierr .ne. 0) goto 9900
       elseif (var .eq. 'ws   ') then
c         write(*,*) '--> ws'
         spline = iplin(time,t_ws,ws,nws,ierr)
         if (ierr .ne. 0) goto 9900
       elseif (var .eq. 'h    ') then
c         write(*,*) '--> h'
         spline = iplin(time,t_h,h,nh,ierr)
         if (ierr .ne. 0) goto 9900
       elseif (var .eq. 'pco2a') then
c         write(*,*) '--> pco2a'
         spline = iplin(time,t_pco2a,pco2a,npco2a,ierr)
         if (ierr .ne. 0) goto 9900
       elseif (var .eq. 'dc13a') then
c         write(*,*) '--> dc13a'
         spline = iplin(time,t_dc13a,dc13a,ndc13a,ierr)
         if (ierr .ne. 0) goto 9900
       elseif (var .eq. 'kz   ') then
c          write(*,*) '--> kz'
        spline = iplin(time,t_kz,kz,nkz,ierr)
         if (ierr .ne. 0) goto 9900
       elseif (var .eq. 'dcdz ') then
c         write(*,*) '--> dcdz'
         spline = iplin(time,t_dcdz,dcdz,ndcdz,ierr)
         if (ierr .ne. 0) goto 9900
      else 
         write(*,*) '--SPLINE: variable not available'
         goto 9900
      endif 
c
c      if (var .eq. 'h    ') then
c         write(*,*) '--SPLINE: time,value : ',time,spline
c         write(*,*) '  values of t_h : ',(t_h(i),i=1,nh)
c         write(*,*) '  values of h   : ',(h(i),i=1,nh)
c      endif
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
 9900 write(*,*) '--SPLINE: an error occured in spline.'
      error = .true.
      return
      end
