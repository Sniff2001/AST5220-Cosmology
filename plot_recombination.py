import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
import pandas as pd

if __name__ == "__main__":
    a = np.loadtxt("recombination.txt").T
    x = a[0]
    t = a[1]
    Xe = a[2]
    ne = a[3]
    tau = a[4]
    dtaudx = a[5]
    ddtauddx = a[6]
    gtilde = a[7]
    dgtildedx = a[8]
    ddgtildeddx = a[9]
    sound_horizon = a[10]
    saha = a[11]
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

    fig, ax = plt.subplots()
    last_scattering = x[np.argmin(np.abs(tau - 1))]
    recombination = x[np.argmin(np.abs(Xe - 0.1))]
    last_scattering_vis = x[np.argmax(gtilde)]
    add_to_plot(ax, x, tau, label=r"$\tau (x)$", yscale="log")
    add_to_plot(ax, x, -dtaudx, label=r"$-\tau ' (x)$", yscale="log")
    add_to_plot(ax, x[500:-500], ddtauddx[500:-500], label=r"$\tau '' (x)$", ylabel="Amplitude", yscale="log")
    ax.axvline(last_scattering, color="k", ls="--", label="Last Scattering")
    ax.legend()
    fig.show()

    fig, ax = plt.subplots()
    add_to_plot(ax, x, Xe, label=r"$X_e$", yscale="log", ylabel=r"$X_e$")
    add_to_plot(ax, x, saha, label=r"$X_e$ (Only Saha)", yscale="log", ylabel=r"$X_e$", ls="--", zorder=0)
    ax.axvline(recombination, color="k", ls="--", label="Recombination")
    ax.set_ylim(bottom=1e-4, top=2e0)
    ax.legend()
    fig.show()

    fig, ax = plt.subplots()
    add_to_plot(ax, x, gtilde/np.max(gtilde), label=r"$\tilde{g}$")
    add_to_plot(ax, x, dgtildedx/np.max(dgtildedx), label=r"$\frac{\mathrm{d}\tilde{g}}{\mathrm{d}x}$", ls="--")
    add_to_plot(ax, x, ddgtildeddx/np.max(ddgtildeddx), label=r"$\frac{\mathrm{d}^2\tilde{g}}{\mathrm{d}x^2}$", ylabel=r"Amplitude", ls="dotted")
    ax.set_xlim(left=-8, right=-5)
    ax.legend()
    fig.show()

    values = np.array([last_scattering, last_scattering_vis, recombination, x[np.nanargmin(np.abs(saha - 0.1))]])
    values_z = 1 / np.exp(values) - 1
    values_t = np.array([t[np.argmin(np.abs(tau - 1))], t[np.argmax(gtilde)], t[np.argmin(np.abs(Xe - 0.1))], t[np.nanargmin(np.abs(saha - 0.1))]]) / (const.kilo * const.year)
    data = np.vstack((np.vstack((values, values_z)), values_t))
    table = pd.DataFrame(data=data, index=["x", "z", "t [kyr]"], columns = ["Last Scattering (OD)", "Last Scattering (Vis)", "Recombination", "Recombination (Saha)"])
    print(table)
    print("Free electrons today: ", Xe[-1])
    print(f"Sound horizon at decoupling (OD): {sound_horizon[np.argmin(np.abs(tau - 1))] / (const.mega * const.parsec):.2f} Mpc")
    print(f"Sound horizon at decoupling (Vis): {sound_horizon[np.argmax(gtilde)] / (const.mega * const.parsec):.2f} Mpc")

