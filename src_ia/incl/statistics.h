c==================== include file statistics.h ==============================
c
c variables for the run statistics
c
c REVISIONS:
c
c date     author    remarks
c
c 08.8.96  ng        first implementation
c 05.9.96  ng        included c_min,c_max,c_min_obs,c_max_obs
c
c variables:
c
c name          type   description
c-----------------------------------------------------------------------------
c 
c chisq         dp s   chi**2 between simulated and observed rate of change
c rsq           dp s   r**2 of fit between simulated and observed rates
c closure       dp s   annual lack of closure
c amplitude     dp s   simulated annual amplitude
c amplitude_obs dp s   observed annual amplitude
c c_min         dp s   simulated minimum C
c c_max         dp s   simulated maximum C
c c_min_obs     dp s   observed minimum C
c c_max_obs     dp s   observed maximum C
c par_sum       dp s   value of parameter misfit
c costfn        dp s   value of cost function
c chisq_term    dp s   value of chi**2 term in cost function
c clos_term     dp s   value of closure term in cost function
c ampl_term     dp s   value of amplitude term in cost function
c par_term      dp s   value of parameter term in cost function
c
c-------------------------------------------------------------------------
c
      double precision chisq,rsq,closure,costfn,amplitude
      double precision amplitude_obs,par_sum
      double precision c_min,c_max,c_min_obs,c_max_obs
      double precision chisq_term,clos_term,ampl_term,par_term



