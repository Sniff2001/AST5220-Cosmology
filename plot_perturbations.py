import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import pandas as pd

if __name__ == "__main__":
    a = np.loadtxt("perturbations_k0.01.txt").T
    x = a[0]
    Theta0_0 = a[1]
    Theta1_0 = a[2]
    Theta2_0 = a[3]
    Phi_0 = a[4]
    Psi_0 = a[5]
    Pi_0 = a[6]
    delta_cdm_0 = a[7]
    delta_b_0 = a[8]
    v_cdm_0 = a[9]
    v_b_0 = a[10]
    #source_0 = a[11]
    #source5_0 = a[12]
    #source50_0 = a[13]
    #source500_0 = a[14]
    a = np.loadtxt("perturbations_k0.1.txt").T
    Theta0_1 = a[1]
    Theta1_1 = a[2]
    Theta2_1 = a[3]
    Phi_1 = a[4]
    Psi_1 = a[5]
    Pi_1 = a[6]
    delta_cdm_1 = a[7]
    delta_b_1 = a[8]
    v_cdm_1 = a[9]
    v_b_1 = a[10]
    # source_1 = a[11]
    # source5_1 = a[12]
    # source50_1 = a[13]
    # source500_1 = a[14]
    a = np.loadtxt("perturbations_k0.001.txt").T
    Theta0_2 = a[1]
    Theta1_2 = a[2]
    Theta2_2 = a[3]
    Phi_2 = a[4]
    Psi_2 = a[5]
    Pi_2 = a[6]
    delta_cdm_2 = a[7]
    delta_b_2 = a[8]
    v_cdm_2 = a[9]
    v_b_2 = a[10]
    # source_2 = a[11]
    # source5_2 = a[12]
    # source50_2 = a[13]
    # source500_2 = a[14]
    z = 1 / np.exp(x) - 1

    def add_to_plot(ax, x, y, xscale="linear", yscale="linear", label=None, xlabel="x", ylabel="y", title="", **kwargs):
        ax.plot(x, y, label=label, **kwargs)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        if title is not None:
            ax.set_title(title)
        ax.grid(ls=":")

    fig, ax = plt.subplots(tight_layout=True)
    add_to_plot(ax, x, delta_cdm_2, yscale="log", label=r"$k = 0.001 \:\mathrm{Mpc^{-1}}$", ylabel="Overdensity",
                title=r"$\delta_{CDM}, \delta_b$", zorder=3)
    add_to_plot(ax, x, delta_cdm_0, yscale="log", label=r"$k = 0.01 \:\mathrm{Mpc^{-1}}$", ylabel="Overdensity",
                title=r"$\delta_{CDM}, \delta_b$", zorder=4)
    add_to_plot(ax, x, delta_cdm_1, yscale="log", label=r"$k = 0.1 \:\mathrm{Mpc^{-1}}$", ylabel="Overdensity",
                title=r"$\delta_{CDM}, \delta_b$", zorder=5)
    add_to_plot(ax, x, delta_b_2, yscale="log", ylabel="Overdensity", title=r"$\delta_{CDM}, \delta_b$", ls="--", c="C0", zorder=6)
    add_to_plot(ax, x, delta_b_0, yscale="log", ylabel="Overdensity", title=r"$\delta_{CDM}, \delta_b$", ls="--", c="C1", zorder=7)
    add_to_plot(ax, x, delta_b_1, yscale="log", ylabel="Overdensity", title=r"$\delta_{CDM}, \delta_b$", ls="--", c="C2", zorder=8)
    add_to_plot(ax, x, 4*Theta0_2, yscale="log", ylabel="Overdensity", title=r"$\delta_{CDM}, \delta_b$", ls=":",
                c="C0", zorder=2)
    add_to_plot(ax, x, 4*Theta0_0, yscale="log", ylabel="Overdensity", title=r"$\delta_{CDM}, \delta_b$", ls=":",
                c="C1", zorder=1)
    add_to_plot(ax, x, 4*Theta0_1, yscale="log", ylabel="Overdensity", title=r"$\delta_{CDM}, \delta_b, \delta_\gamma$", ls=":",
                c="C2", zorder=0)
    ax.legend()
    fig.show()
    fig, ax = plt.subplots(tight_layout=True)
    add_to_plot(ax, x, v_cdm_2, yscale="log", label=r"$k = 0.001 \:\mathrm{Mpc^{-1}}$", ylabel="Velocity",
                title=r"$v_{CDM}, v_b$", zorder=3)
    add_to_plot(ax, x, v_cdm_0, yscale="log", label=r"$k = 0.01 \:\mathrm{Mpc^{-1}}$", ylabel="Velocity",
                title=r"$v_{CDM}, v_b$", zorder=4)
    add_to_plot(ax, x, v_cdm_1, yscale="log", label=r"$k = 0.1 \:\mathrm{Mpc^{-1}}$", ylabel="Velocity",
                title=r"$v_{CDM}, v_b$", zorder=5)
    add_to_plot(ax, x, v_b_2, yscale="log", ylabel="Velocity", title=r"$v_{CDM}, v_b$", ls="--", c="C0", zorder=6)
    add_to_plot(ax, x, v_b_0, yscale="log", ylabel="Velocity", title=r"$v_{CDM}, v_b$", ls="--", c="C1", zorder=7)
    add_to_plot(ax, x, v_b_1, yscale="log", ylabel="Velocity", title=r"$v_{CDM}, v_b$", ls="--", c="C2", zorder=8)
    add_to_plot(ax, x, -3*Theta1_2, yscale="log", ylabel="Velocity", title=r"$v_{CDM}, v_b$", ls=":", c="C0", zorder=2)
    add_to_plot(ax, x, -3*Theta1_0, yscale="log", ylabel="Velocity", title=r"$v_{CDM}, v_b$", ls=":", c="C1", zorder=1)
    add_to_plot(ax, x, -3*Theta1_1, yscale="log", ylabel="Velocity", title=r"$v_{CDM}, v_b, v_\gamma$", ls=":", c="C2", zorder=0)
    ax.legend()
    fig.show()
    fig, ax = plt.subplots(tight_layout=True)
    add_to_plot(ax, x, Theta0_2, label=r"$k = 0.001 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Theta_0$")
    add_to_plot(ax, x, Theta0_0, label=r"$k = 0.01 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Theta_0$")
    add_to_plot(ax, x, Theta0_1, label=r"$k = 0.1 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Theta_0$")
    ax.legend()
    fig.show()
    fig, ax = plt.subplots(tight_layout=True)
    add_to_plot(ax, x, Theta1_2, label=r"$k = 0.001 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Theta_1$")
    add_to_plot(ax, x, Theta1_0, label=r"$k = 0.01 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Theta_1$")
    add_to_plot(ax, x, Theta1_1, label=r"$k = 0.1 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Theta_1$")
    ax.legend()
    fig.show()
    fig, ax = plt.subplots(tight_layout=True)
    add_to_plot(ax, x, Theta2_2, label=r"$k = 0.001 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Theta_2$")
    add_to_plot(ax, x, Theta2_0, label=r"$k = 0.01 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Theta_2$")
    add_to_plot(ax, x, Theta2_1, label=r"$k = 0.1 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Theta_2$")
    ax.legend()
    fig.show()
    fig, ax = plt.subplots(tight_layout=True)
    add_to_plot(ax, x, Phi_2, label=r"$k = 0.001 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Phi$")
    add_to_plot(ax, x, Phi_0, label=r"$k = 0.01 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Phi$")
    add_to_plot(ax, x, Phi_1, label=r"$k = 0.1 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Phi$")
    ax.legend()
    fig.show()
    fig, ax = plt.subplots(tight_layout=True)
    add_to_plot(ax, x, Phi_2 + Psi_2, label=r"$k = 0.001 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Phi + \Psi$")
    add_to_plot(ax, x, Phi_0 + Psi_0, label=r"$k = 0.01 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Phi + \Psi$")
    add_to_plot(ax, x, Phi_1 + Psi_1, label=r"$k = 0.1 \:\mathrm{Mpc^{-1}}$", ylabel="Amplitude", title=r"$\Phi + \Psi$")
    ax.legend()
    fig.show()

