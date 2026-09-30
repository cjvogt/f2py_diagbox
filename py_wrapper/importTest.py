def main():
    import numpy as np
    import diagbox
    import diagbox_plots as diag_plot
    import os, sys
    import numpy as np
    import matplotlib.pyplot as plt

    incl_adv = False
    plot_out = False
    # c     
    # c --------------------------------------------------------------------
    # c     simulation control parameters
    # c --------------------------------------------------------------------
    # c     
    tstart = 1983.0           #1983.0d0
    tend = 2001.0             #2001.d0
    ts = 0.0025               #0.0025d0 
    # c     
    # c --------------------------------------------------------------------
    # c     air-sea gasexchange parameters
    # c --------------------------------------------------------------------
    # c     
    pv_relship = 'wa'

    D_fco2_corr = 0.0          #0.0d0 
    sg_D_fco2_corr = 10.0      #10.d0 
    sg_D_fco2_corr = 1.0       #1.d0 

    gasex_fact = 1.0           #1.0d0 
    sg_gasex_fact = 0.1        #0.1d0

    ws_coeff = 'ih'

    # c     
    # c --------------------------------------------------------------------
    # c     diffusion - entrainment parameters
    # c --------------------------------------------------------------------
    # c     h_coeff = 'a'

    ent_scheme = 'ma'
    const_kz = False

    kz_const = 2.0e-5          #2.d-5
    sg_kz_const = 1.0e-5       #1.d-5

    const_vertgrad = False
    dCdz_const = 0.45          #0.45d0

    D_dcdz = 0.0               #0.0d0
    sg_D_dcdz = 0.1            #0.1d0

    h_th = 200.0               #200.d0
    c_th =  2066.31            #2066.31d0
    dc13_th = 1.28             #1.28d0

    const_lent = True
    lent_const = 12.0          #12.d0
    sg_lent_const = 6.0        #6.d0
    sg_lent_const = 3.0        #3.d0

    dlent_dh = 0.11
    sg_dlent_dh = 0.4

    ddC13_dC_bml = -0.0048     #-0.0048d0 
    sg_ddC13_dC_bml = 0.0008   #0.0008d0

    ent_fact = 1.0             #1.0d0
    diff_fact = 1.0            #1.0d0
    sg_diff_fact = 0.5         #0.5d0
    sg_ent_fact = 0.5          #0.5d0
    # c     
    # c --------------------------------------------------------------------
    # c     advection parameters
    # c --------------------------------------------------------------------
    # c     
    dC_dh = 1.1e-5             #1.1d-5
    sg_dC_dh = 0.3e-5          #0.3d-5

    ddC13_dh = -1.0e-7         #-1.0d-7
    sg_ddC13_dh = 1.0e-7       #1.0d-7

    u = 0.05                   #0.05d0
    sg_u = 0.05                #0.05d0
    # c     
    # c --------------------------------------------------------------------
    # c     net community production parameters
    # c --------------------------------------------------------------------
    # c     
    D_dc13_org = 0.0           #0.0d0
    sg_D_dc13_org = 3.0        #3.0d0

    #Provide name for this experiment
    exp_name = ""

    #Provide notes for this experiment
    notes = ""

    #Specify simulation mode. Options: st, mc, op, sn
    sim_mode = "st"

    #Specify solution scheme. Options: std, adv, phy, kz
    sol_scheme = "std"

    #SKIPPED MONTE CARLO FOR NOW

    #Specify piston-velocity. Options: lm, ma
    pv_relship = "lm"

    if pv_relship == "lm":
        gasex_fact = 1.7447
    elif pv_relship == "wa":
        gasex_fact = 1.0


    D_fco2_corr = 0 #DEFAULT

    #Specific how fCO2(oc) should be calculated. Options: True or False
    compl_foc2_online = True

    #Specify entrainment scheme: Options: in, mc, ep, ma
    ent_scheme = "in"

    dt_ent_re = 0.0         #0.d0

    #Constant entrainment
    const_lent = True

    #Entrainment length scale
    lent_const = 1 #1 meter

    #Entrainment length to mixed-layer depth
    dlent_dh = 0.1

    #Entrainment multiplication
    ent_fact = 1

    #Specify if you want to use seasonally constant Kz or varying. Options True or False
    const_Kz = True

    kz_const = 2e-5         #2.d-5 m2 s-1

    #Specify if you want to use seasonal or constant veritcal gradients: Options True (Constant) or False (Seasonal)
    const_vertgrad = False

    plot_out = False

    nstep = int((tend-tstart)/ts)

    #dc13_coeff, h_coeff = 1
    #error = False
    default_sim = diagbox.run_simulation(sol_scheme, tstart, tend, ts, nstep, 1, pv_relship, 
                                         D_fco2_corr, gasex_fact, ws_coeff, compl_foc2_online,
                                         ent_scheme, 1, h_th, c_th, dc13_th, dt_ent_re, const_lent, lent_const,
                                         dlent_dh, ent_fact, const_Kz, kz_const, const_vertgrad, dCdz_const,
                                          D_dcdz, ddC13_dC_bml, diff_fact, incl_adv, dC_dh, ddC13_dh, u, D_dc13_org, 0)

    out = default_sim
    assert len(out) == 35
# ...then the (u_var, d_c, ...) = out unpack and the plotting code
    
    # print(diagbox.run_simulation.__doc__)
    # print(default_sim)    
    # for each in default_sim:
    #     print(each)


    

    (u_var, d_c, d_c13, fluxes, c_obs, dc13_obs, h_obs, temp_obs, sal_obs, salk_obs,
    fco2_o_obs, pco2_a_obs, fco2_a_obs, pco2_a_dry, dc13_a_obs, ws_obs, dcdz_obs,
    kz_obs, c_sim, int_rates, int_fluxes, nyears, chisq, rsq, closure, amplitude,
    amplitude_obs, c_min, c_max, c_min_obs, c_max_obs, costfn, avg_fract,
    dc13_prc, avg_kex) = out

    diag_plot.plot_run_simulation(out, tstart, tend, ts, nstep)

    # # valid length: last nonzero entry of c_sim
    # n = np.flatnonzero(c_sim)[-1] + 1
    # print("valid points:", n, " nstep:", nstep, " stride ~", nstep / n)

    # t = tstart + (tend - tstart) * np.arange(n) / n       # verify against the Fortran output stride
    # clean = lambda a: np.where(np.isclose(a, -99.99, atol=0.01), np.nan, a)

    # # sanity: where does the simulated DIC spike?
    # i = np.argmax(c_sim[:n])
    # print(f"c_sim max {c_sim[i]:.1f} at t={t[i]:.3f};  c_obs max {c_obs[:n].max():.1f}")

    # fig, ax = plt.subplots(4, 1, figsize=(10, 12), sharex=True)

    # ax[0].plot(t, c_obs[:n], label="obs")
    # ax[0].plot(t, c_sim[:n], label="sim", lw=0.8)
    # ax[0].set_ylabel("DIC"); ax[0].legend()

    # ax[1].plot(t, c_sim[:n] - c_obs[:n]); ax[1].axhline(0, c="k", lw=0.5)
    # ax[1].set_ylabel("sim - obs")

    # ax[2].plot(t, temp_obs[:n], label="T"); ax[2].set_ylabel("temp (°C)")
    # ax2 = ax[2].twinx(); ax2.plot(t, h_obs[:n], c="C1", label="MLD"); ax2.set_ylabel("MLD (m)")

    # for k in range(7):
    #     ax[3].plot(t, d_c[k, :n], lw=0.6, label=f"term {k}")
    # ax[3].set_ylabel("d_c terms"); ax[3].set_xlabel("year"); ax[3].legend(ncol=7, fontsize=7)

    # plt.tight_layout(); plt.show()

    # # annual budget
    # yrs = tstart + np.arange(nyears)
    # plt.figure(figsize=(9, 4))
    # for k in range(7):
    #     plt.plot(yrs, int_rates[:nyears, k], marker="o", label=f"term {k}")
    # plt.legend(ncol=4, fontsize=8); plt.xlabel("year"); plt.ylabel("annual integrated rate"); plt.show()

if __name__ == "__main__":
    main()
