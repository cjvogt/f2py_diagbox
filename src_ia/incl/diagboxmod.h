c===================== include file diagboxmod.h =========================
c
c general variables of diagnostic box model
c
c REVISIONS:
c
c date     author    remarks
c
c 31.7.96  ng        first implementation
c  2.9.96  ng        included inverse mode
c
c variables:
c
c name          type   description
c---------------------------------------------------------------------------
c               
c exp_name      c  s   name of experiment 
c notes         c  s   notes to add to experiment on output
c sim_mode      c  s   character specifying simulation mode
c                       <st>  : standard simulation
c                       <mc>  : monte carlo simulation
c                       <op>  : parameter optimization
c                       <sn>  : sensitivity mode
c                       <in>  : inverse calculation of u and dC13/dx
c sol_scheme    c  s    <std> : original standard, solving for dC/dt
c               c  s    <adv> : solve for u velocity
c               c  s    <phy> : solve for Kz + u
c               c  s    <kz>  : solve for Kz
c
c dc13_coeff    c  s   character specifying which dc13 coefficients are
c                        used
c                       <o1>  : original data with harmonic fit of 1st order
c                       <o2>  : original data with harmonic fit of 2nd order
c                       <a2>  : artificial curve with 2nd order harm fit
c def_par       l  s   logical specifying use of default parameters for run
c incl_adv      l  s   logical specifying whether advection is included
c                        or not
c plot_out      l  s   logical indicating that plots should be generated
c
c===========================================================================
c
      character exp_name*10,notes*30,dc13_coeff*2
      logical def_par,incl_adv,plot_out
      character sim_mode*2,sol_scheme*3
      
