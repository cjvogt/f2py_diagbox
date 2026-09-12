c =====================================================================
c
c     SUBROUTINE WRITE_OP_OUTPUT
c
c     PURPOSE: writes the results of the optimization to the unit
c                specified by unit
c     
c
c     VARIABLES:
c 
c     IN:  unit            : unit number to write to
c
c          p_init          : value of parameters at the beginning
c          p_opt           : value of parameters at the optimum
c          sq_p            : uncertainties of parameters
c          chisq_init      : value of chi**2 at the beginning
c          rsq_init        : value of r**2 at the beginning
c          closure_init    : value of closure at the beginning
c          amplitude_init  : value of amplitude at the beginning
c          c_min_init      : value of C_min at the beginning
c          c_max_init      : value of C_max at the beginning
c          costfn_init     : value of costfn at the beginning
c          chisq           : chi**2 between simulated and obs rates
c          rsq             : r**2 between simulated and obs rates
c          closure         : annual lack of closure
c          amplitude       : annual amplitude
c          c_min           : minimum C
c          c_max           : maximum C
c          amplitude_obs   : observed amplitude
c          c_min_obs       : observed minimum C
c          c_max_obs       : observed maximum C
c          costfn          : value of cost function (without parameter
c                               terms) at the end
c          iter            : number of iterations needed
c          fret            : value of cost function (with parameter
c                               terms) at the end
c          hessian         : hessian matrix
c          covm            : inverse of hessian matrix, covariance matrix
c
c     OUT: error           : logical indicating that an error occured
c
c     REVISIONS:
c 
c     date      author  remarks
c     
c     03.08.96   ng     first implementation
c     08.08.96   ng     included statistics and updated on entrainment
c     05.09.96   ng     adapted to new def of cost function
c 
c =====================================================================
c
      subroutine write_op_output(unit,
     $     exp_name,
     $     p_opt,p_init,sg_p,
     $     chisq_init,rsq_init,closure_init,
     $     amplitude_init,
     $     c_min_init,c_max_init,costfn_init,
     $     chisq,rsq,closure,amplitude,
     $     c_min,c_max,
     $     amplitude_obs,
     $     c_min_obs,c_max_obs,costfn,
     $     iter,fret,hessian,covm,
     $     error)
c
c --------------------------------------------------------------------
c     global variables
c --------------------------------------------------------------------
c     
      implicit none
# include "diagboxmod.h"
# include "conj_grad.h"
# include "statistics.h"
# include "op_statistics.h"
# include "sg_params.h"

      logical error
      integer unit
c
c --------------------------------------------------------------------
c     local variables
c --------------------------------------------------------------------
c     
      character date*17
      integer np,i,j
      double precision costfntmp,par_term_init,par_sum_init
      double precision chisq_term_init
      double precision clos_term_init
      double precision ampl_term_init
      double precision corr(nparams,nparams)
c
c --------------------------------------------------------------------
c     external functions
c --------------------------------------------------------------------
c 
      double precision cost_function
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
     $            ' Results of parameter optimization of the ',/
     $            '    seasonal diagnostic box model',/
     $            '    of the seasonal carbon cycle at ',/
#ifdef BATS
     $            '    the U.S. JGOFS station BATS ',//
#elif BERM
     $            '    the Hydrostation 'S' ',//
#endif
     $            '    written by Nicolas Gruber, Aug. 1996 ',/
     $            '    this is version 1.0 ',//
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
c     write parameters, denormalize them
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)')' optimization of parameters :'
c
      do np = 1,nparams
         write(unit,'(A,i2,1x,a,g10.4,a,g10.4)')
     $        '    parameter : ',np,';  initial : ',p_init(np),
     $        ';  optimized : ',sg_p(np)*p_opt(np) + p_init(np)
      enddo
c
c ---------------------------------------------------------------------
c     write statistics
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' optimization of run statistics : '
c
      write(unit,'(A,F8.4,a,f8.4)') '    chi**2  : initial : ',
     $     chisq_init,';  optimized : ',chisq
      write(unit,'(A,F8.4,a,f8.4)') '    r**2    : initial : ',
     $     rsq_init  ,';  optimized : ',rsq
      write(unit,'(A,F8.2,a,f8.2)') '    closure : initial : ',
     $     closure_init, ';  optimized : ',closure
      write(unit,'(A,F8.2,a,f8.2)') '    ampl.   : initial : ',
     $     amplitude_init, ';  optimized : ',amplitude
      write(unit,'(A,F8.2,a,f8.2)') '    min C   : initial : ',
     $     c_min_init - c_min_obs,
     $                   ';  optimized : ',c_min - c_min_obs
      write(unit,'(A,F8.2,a,f8.2)') '    max C   : initial : ',
     $     c_max_init - c_max_obs,
     $                   ';  optimized : ',c_max - c_max_obs
c
      write(unit,'(/A,F10.8,a,f10.8)') 
     $     '    cost function (w/o par)  : initial : ',
     $     costfn_init,';  optimized : ',costfn
      write(unit,'(A,f10.8,a,f10.8)')
     $     '    cost function (with par) : initial : ',
     $     costfn_init,';  optimized : ',fret
      write(unit,'(A,i2)')
     $     '    number of iterations     : ',iter
c
c ---------------------------------------------------------------------
c     write detailed statistics
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' details of cost function : '
c
c
c --------------------------------------------------------------------
c     calculate cost function : part of parameter deviation
c -------------------------------------------------------------------- 
c
      par_sum = 0.d0
      do np = 1,nparams
c         par_sum = par_sum + ((p_opt(np) - p_init(np))/sg_p(np))**2
         par_sum = par_sum + p_opt(np)**2
      enddo
c
      costfntmp = cost_function(chisq,closure,par_sum,nparams,
     $     amplitude,amplitude_obs,
     $     c_min,c_max,c_min_obs,c_max_obs,
     $     chisq_term,clos_term,par_term,ampl_term)
c
      par_sum_init = 0.d0
      costfntmp = cost_function(chisq_init,closure_init,
     $     par_sum_init,nparams,
     $     amplitude_init,amplitude_obs,
     $     c_min_init,c_max_init,c_min_obs,c_max_obs,
     $     chisq_term_init,clos_term_init,
     $     par_term_init,ampl_term_init)
c
      write(unit,'(A,F10.8,a,f10.8)') '    chi**2 term  : initial : ',
     $     chisq_term_init,';  optimized : ',chisq_term
      write(unit,'(A,F10.8,a,f10.8)') '    closure term : initial : ',
     $     clos_term_init, ';  optimized : ',clos_term
      write(unit,'(A,F10.8,a,f10.8)') '    par term     : initial : ',
     $     par_term_init  ,';  optimized : ',par_term
      write(unit,'(A,F10.8,a,f10.8)') '    ampl term    : initial : ',
     $     ampl_term_init  ,';  optimized : ',ampl_term
c
c ---------------------------------------------------------------------
c     write hessian matrix
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' hessian matrix ',
     $     '(d2J/(dx_i dx_j)) : '
c
      write(unit,'(10x,10(8x,i2,1x))')(j,j=1,nparams)
c

      do i = 1,nparams
         write(unit,'(a6,i2,2x,10(g10.4,1x))')
     $        'row : ',i,(hessian(i,j),j=1,nparams)
      enddo      
c
c ---------------------------------------------------------------------
c     write covariance matrix
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' covariance matrix ',
     $     '(inverse of hessian matrix) : '
c
      write(unit,'(10x,10(8x,i2,1x))')(j,j=1,nparams)
c

      do i = 1,nparams
         write(unit,'(a6,i2,2x,10(g10.4,1x))')
     $        'row : ',i,(covm(i,j),j=1,nparams)
      enddo      
c
c ---------------------------------------------------------------------
c     write correlation matrix
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/A/)') ' correlation matrix :'
c
      write(unit,'(10x,10(3x,i2,1x))')(j,j=1,nparams)
c
      do i = 1,nparams
         do j = 1,nparams
            if ((covm(i,i) .ne. 0.d0) .and. 
     $           (covm(j,j) .ne. 0.d0)) then
               corr(i,j) = covm(i,j)/
     $              (dsqrt(covm(i,i))*dsqrt(covm(j,j)))
            else
               corr(i,j) = -99.99
            endif
         enddo
         write(unit,'(a6,i2,2x,10(f7.3,1x))')
     $        'row : ',i,(corr(i,j),j=1,nparams)
      enddo      
c
c
c ---------------------------------------------------------------------
c     write parameter std deviations and relative uncertainties
c ---------------------------------------------------------------------
c
      write(unit,9800)
c
      write(unit,'(/2A/)') ' estimated one-std deviation uncertainty',
     $     ' and relative uncertainty'
c
      write(unit,'(a/)') 'parameter      1-sg        rel. uncert.'
c
      do np = 1,nparams
         if ((sg_p(np)*p_opt(np) + p_init(np)) .ne. 0.d0) then
c
            write(unit,'(i2,2x,g10.4,2x,g10.4,2x,g10.4)')
     $           np,
     $           sg_p(np)*p_opt(np) + p_init(np),
     $           dsqrt(covm(np,np)),
     $           dsqrt(covm(np,np))/(sg_p(np)*p_opt(np) + p_init(np))
         else
            write(unit,'(i2,2x,g10.4,2x,g10.4,2x,g10.4)')
     $           np,
     $           sg_p(np)*p_opt(np) + p_init(np),
     $           dsqrt(covm(np,np)),
     $           -99.99
         endif 
      enddo      
c
      write(unit,9800)
c
c ---------------------------------------------------------------------
c     format statements
c ---------------------------------------------------------------------
c
 9800 format(/'-------------------------------------------------')
c
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











