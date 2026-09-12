c==================== include file op_statistics.h ==========================
c
c variables for the optimization run statistics
c
c REVISIONS:
c
c date     author    remarks
c
c 08.8.96  ng        first implementation
c
c variables:
c
c name           type   description
c-----------------------------------------------------------------------------
c 
c chisq_init     dp s   chi**2 between simulated and observed rate of change at
c                            the beginning of the optimization
c rsq_init       dp s   r**2 of fit between simulated and observed rates at
c                            the beginnning of the optimization
c closure_init   dp s   annual lack of closure at the beginning of the optim.
c amplitude_init dp s   annual amplitude at the beginning of the optimization
c costfn_init    dp s   value of cost function at the beginning of the optim.
c c_min_init     dp s   value of minimum C at the beginning of the optim.
c c_max_init     dp s   value of maximum C at the beginning of the optim.
c
c-------------------------------------------------------------------------
c
      double precision chisq_init,rsq_init,closure_init,costfn_init
      double precision amplitude_init,c_min_init,c_max_init



