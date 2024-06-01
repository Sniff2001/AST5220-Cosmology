import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const

if __name__ == "__main__" or True:
    a = np.loadtxt("cosmology.txt").T
    x = a[0]
    t = a[1]
    eta = a[2]
    Hp = a[3]
    dHpdx = a[4]
    ddHpddx = a[5]
    OmegaB = a[6]
    OmegaCDM = a[7]
    OmegaLambda = a[8]
    OmegaR = a[9]
    OmegaNu = a[10]
    OmegaK = a[11]
    luminosity_distance = a[12]
    z = 1 / np.exp(x) - 1
    mtoMpc = 1 / (const.mega * const.parsec)
    HpMpc = Hp / (const.kilo * const.hecto * mtoMpc)

    supernovadata = np.loadtxt("data/supernovadata.txt").T
    z_data = supernovadata[0]
    luminosity_distance_data_Gpc = supernovadata[1]
    error_data_Gpc = supernovadata[2]

    """
    chidata = np.loadtxt("results_supernovafitting.txt").T
    chi2 = chidata[0]
    hs = chidata[1]
    OmegaMs = chidata[2]
    OmegaKs = chidata[3]
    """

    #### Planck data (2018) arXiv:1807.06209 [astro-ph.CO]
    OmegaK_Planck = 0.001
    OmegaK_Planckerr = 0.002
    OmegaM_Planck = 0.315
    OmegaM_Planckerr = 0.007
    OmegaLambda_Planck = 1 - OmegaK_Planck - OmegaM_Planck
    OmegaLambda_Planckerr = OmegaM_Planckerr + OmegaK_Planckerr
    H0_Planck = 67.37
    H0_Planckerr = 0.54

    #### S. Dhawan et al. https://doi.org/10.1051/0004-6361/201731501
    H0_Dhawan = 72.8
    H0_Dhawanerr = 1.6

    mr_eq = lambda m,r: r/m
    mc_eq = lambda m,c: np.cbrt(m/c)
    acc_exp = lambda m,c: np.cbrt(m/(2*c))

    def add_to_plot(ax, x, y, xscale="linear", yscale="linear", label=None, xlabel="x", ylabel="y", title="", **kwargs):
        ax.plot(x, y, label=label, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        if title is not None:
            ax.set_title(title)
        ax.grid(ls=":")


    # Quick test (or demonstration if you will)
    # fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2, tight_layout=True)
    # add_to_plot(ax1, x, Hp/(const.kilo*mtoMpc), yscale="log", ylabel=r"$\mathcal{H}$")
    # add_to_plot(ax2, x, eta*mtoMpc, yscale="log", ylabel=r"$\eta$")
    # add_to_plot(ax3, x, eta*Hp/const.c, ylabel=r"$\eta \mathcal{H} / c$")
    # add_to_plot(ax4, x, OmegaR + OmegaNu, label=r"$\Omega_{relativistic}$")
    # add_to_plot(ax4, x, OmegaB + OmegaCDM, label=r"$\Omega_M$")
    # add_to_plot(ax4, x, OmegaLambda, label=r"$\Omega_\Lambda$", ylabel=r"$\Omega$")
    # ax4.legend()
    # add_to_plot(ax5, OmegaB + OmegaCDM, OmegaLambda, xlabel=r"$\Omega_M$", ylabel=r"$\Omega_\Lambda$")
    # fig.show()

    fig, ((ax1), (ax2), (ax3)) = plt.subplots(3, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    add_to_plot(ax1, x, 1 / Hp * dHpdx, xlabel=None, ylabel=r"$\mathcal{H}^{-1} \mathcal{H}'(x)$")
    ax1.axhline(-1, ls="--", alpha=0.4, c="C1", label="Radiation Domination")
    ax1.axhline(-0.5, ls="--", alpha=0.4, c="C2", label="Matter Domination")
    ax1.axhline(1, ls="--", alpha=0.4, c="C3", label="Cosmological Constant Domination")
    ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.2), frameon=False, ncol=3, fontsize=8)
    add_to_plot(ax2, x, 1 / Hp * ddHpddx, xlabel=None, ylabel=r"$\mathcal{H}^{-1} \mathcal{H}''(x)$")
    ax2.axhline(1, ls="--", alpha=0.4, c="C1")
    ax2.axhline(0.25, ls="--", alpha=0.4, c="C2")
    ax2.axhline(1, ls="-.", alpha=0.4, c="C3")
    add_to_plot(ax3, x, eta * Hp / const.c, ylabel=r"$\eta \mathcal{H} c^{-1}$")
    ax3.axhline(1, ls="--", alpha=0.4, c="C1")
    fig.show()

    fig, ax = plt.subplots()
    add_to_plot(ax, x, HpMpc, yscale="log", ylabel=r"$\mathcal{H} \quad \mathrm{[10^2 \: km \: s^{-1} \: Mpc^{-1}]}$",
                title="Conformal Hubble Parameter")
    #ACC_EXP = acc_exp(OmegaB[-1] + OmegaCDM[-1], OmegaLambda[-1])
    #ax.axvline(np.log(ACC_EXP), color="C5", ls="--",
    #           label="Acceleration starts")
    fig.show()

    fig, ((ax1), (ax2)) = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    add_to_plot(ax1, x, t / (const.giga * const.year), yscale="log", xlabel=None, ylabel=r"$t \quad \mathrm{[Gyr]}$", title=rf"$t_0 = {t[-1] / (const.giga * const.year):.2f} $ Gyr $\qquad \eta_0 c^{{-1}} = {eta[-1] / (const.c * const.giga * const.year):.2f} $ Gyr")
    add_to_plot(ax2, x, eta / const.c / (const.giga * const.year), yscale="log", ylabel=r"$\eta c^{-1} \quad \mathrm{[Gyr]}$")
    fig.show()
    print(f"Age of the universe is {t[-1] / (const.giga * const.year):.2f} Gyr.")

    fig, ax = plt.subplots()
    MR = mr_eq(OmegaB[-1] + OmegaCDM[-1], OmegaR[-1] + OmegaNu[-1])
    MC = mc_eq(OmegaB[-1] + OmegaCDM[-1], OmegaLambda[-1])
    ACC_EXP = acc_exp(OmegaB[-1] + OmegaCDM[-1], OmegaLambda[-1])
    print(f"Radiation-Matter Eq.: x = {np.log(MR):.2f}, z = {1/MR - 1:.2f}, t = {np.interp(np.log(MR), x, t / (const.kilo * const.year)):.2f} kyr.")
    print(f"Matter-Dark Energy Eq.: x = {np.log(MC):.2f}, z = {1 / MC - 1:.2f}, t = {np.interp(np.log(MC), x, t / (const.giga * const.year)):.2f} Gyr.")
    print(f"Acceleration starts: x = {np.log(ACC_EXP):.2f}, z = {1 / ACC_EXP - 1:.2f}, t = {np.interp(np.log(ACC_EXP), x, t / (const.giga * const.year)):.2f} Gyr.")
    add_to_plot(ax, x, OmegaR + OmegaNu, label=r"$\Omega_{relativistic}$", color="C0")
    ax.axvline(np.log(MR), color="C3", ls="--",
               label="Radiation-Matter Eq.")
    add_to_plot(ax, x, OmegaB + OmegaCDM, label=r"$\Omega_M$", color="C1")
    ax.axvline(np.log(MC), color="C4", ls="--",
               label="Matter-Dark Energy Eq.")
    add_to_plot(ax, x, OmegaLambda, label=r"$\Omega_\Lambda$", ylabel=r"$\Omega$", color="C2")
    ax.axvline(np.log(ACC_EXP), color="C5", ls="--",
               label="Acceleration starts")

    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), frameon=False, ncol=3)
    fig.show()

    fig, ax = plt.subplots()
    ax.invert_xaxis()
    ax.set_ylim(3, 8)

    add_to_plot(ax, z, luminosity_distance / (const.giga * const.parsec) / np.insert(z[:-1], -1, z[-2]), label="Simulation", xscale="log",
                xlabel=r"$z$", ylabel=r"Luminosity Distance $d_L/z$ $\mathrm{"r"[Gpc]}$")
    ax.errorbar(z_data, luminosity_distance_data_Gpc / z_data, yerr=error_data_Gpc / z_data, label="Supernova Data", capsize=3, fmt='.',
                ecolor='k', mfc='k', mec='k')
    ax.set_xlim(8e-3, 1.5)
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), frameon=False, ncol=2)
    fig.show()

    fig, ax = plt.subplots()
    ax.grid(ls=":")
    ax.plot(OmegaB + OmegaCDM, 1 - OmegaB - OmegaCDM, ls="--", color="k", label="Flat Universe")
    OmegaMs_1sigma = OmegaMs[chi2 - np.min(chi2) < 3.53]
    OmegaK_1sigma = OmegaKs[chi2 - np.min(chi2) < 3.53]
    OmegaLambda_1sigma = 1 - OmegaMs_1sigma - OmegaK_1sigma
    ax.scatter(OmegaMs_1sigma, OmegaLambda_1sigma, color="C0", label=r"$1\sigma$")
    ax.legend()
    ax.set_xlabel(r"$\Omega_{m,0}$")
    ax.set_ylabel(r"$\Omega_{\Lambda,0}$")
    fig.show()

    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, tight_layout=True, figsize=[8, 4.8])
    H0s = 100 * hs
    OmegaLambdas = 1 - OmegaMs - OmegaKs
    probability_distribution = lambda x: 1 / (np.std(x) * np.sqrt(2 * np.pi)) * np.exp(-0.5 * ((x - np.mean(x)) / np.std(x)) ** 2)

    ax1.hist(H0s[200:], density=True, bins=30)
    ax1.plot(np.sort(H0s[200:], axis=None), probability_distribution(np.sort(H0s[200:], axis=None)))
    ax1.axvline(H0_Planck, color="k", ls="--", label="Planck (2018)")
    ax1.axvspan(H0_Planck - H0_Planckerr, H0_Planck + H0_Planckerr, color="C2",
                alpha=0.2)
    ax1.axvline(H0_Dhawan, color="k", ls="dotted", label="S. Dhawan et al. (2018)")
    ax1.axvspan(H0_Dhawan - H0_Dhawanerr, H0_Dhawan + H0_Dhawanerr, color="C3",
                alpha=0.2)
    ax1.set_xlabel(r"$H_0 \quad \mathrm{[km\: s^{-1} \: Mpc^{-1}]}$")
    handles, labels = ax1.get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.55, 1.025), frameon=False, ncols=2)
    #ax1.legend(loc='upper right', frameon=False)
    ax2.hist(OmegaLambdas[200:], density=True, bins=30)
    ax2.plot(np.sort(OmegaLambdas[200:], axis=None), probability_distribution(np.sort(OmegaLambdas[200:], axis=None)))
    ax2.axvline(OmegaLambda_Planck, color="k", ls="--", label="Planck (2018)")
    ax2.axvspan(OmegaLambda_Planck - OmegaLambda_Planckerr, OmegaLambda_Planck + OmegaLambda_Planckerr, color="C2", alpha=0.2)
    #ax2.legend(loc='upper right', frameon=False)
    ax2.set_xlabel(r"$\Omega_{\Lambda,0}$")
    ax3.hist(OmegaKs[200:], density=True, bins=30)
    ax3.plot(np.sort(OmegaKs[200:], axis=None), probability_distribution(np.sort(OmegaKs[200:], axis=None)))
    ax3.axvspan(OmegaK_Planck - OmegaK_Planckerr, OmegaK_Planck + OmegaK_Planckerr, color="C2", alpha=0.2)
    ax3.axvline(OmegaK_Planck, color="k", ls="--", label="Planck (2018)")
    ax3.set_xlabel(r"$\Omega_{k,0}$")
    #ax3.legend(loc='upper left', frameon=False)
    ax4.hist(OmegaMs[200:], density=True, bins=30)
    ax4.plot(np.sort(OmegaMs[200:], axis=None), probability_distribution(np.sort(OmegaMs[200:], axis=None)))
    ax4.set_xlabel(r"$\Omega_{m,0}$")
    ax4.axvspan(OmegaM_Planck - OmegaM_Planckerr, OmegaM_Planck + OmegaM_Planckerr, color="C2", alpha=0.2)
    ax4.axvline(OmegaM_Planck, color="k", ls="--", label="Planck (2018)")
    #ax4.legend(loc='upper right', frameon=False)
    fig.show()
