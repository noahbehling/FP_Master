import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp
import os


def lin(x, A, B):
    return A * x + B


def line(x, B):
    return B


# Erster Auswertungspart: Linearverstärker
R1 = [1000, 1000, 15000]
R2 = [100000, 150000, 100000]
borders = [12, 13, 14]


for i in range(3):
    f, U_a, delta_t = np.genfromtxt(
        f"data/linearverstaerker_{i + 1}.txt", delimiter=", ", unpack=True,
    )
    log_f = np.log10(f)
    V = U_a / 0.05
    log_V = np.log10(V)
    ln = np.linspace(log_f[0] - 1, log_f[borders[i]])
    ln2 = np.linspace(log_f[borders[i]], max(log_f) + 1)
    plt.plot(log_f[: borders[i]], log_V[: borders[i]], "k+", label="Messwerte")
    plt.plot(log_f[borders[i] :], log_V[borders[i] :], "C5+", label="Messwerte")

    params, cov = curve_fit(line, log_f[: borders[i]], log_V[: borders[i]])

    unc = np.sqrt(np.diag(cov))
    err = unp.uarray(params, unc)
    print(f"Fit-Parameter d. Plateau-Fits: B = {err}")
    print(f"Leerlaufverstärkung = {10**err}")
    V_0 = 10 ** err
    plt.axhline(
        params, xmax=0.65 if i == 2 else 0.5, label="Fit: Plateau",
    )

    params, cov = curve_fit(lin, log_f[borders[i] :], log_V[borders[i] :])
    unc = np.sqrt(np.diag(cov))
    err = unp.uarray(params, unc)
    print(f"Fit-Parameter d. abfallenden Fits: A = {err[0]}, B = {err[1]}")
    print(f"Grenzfrequenz f = {10**((np.log10(1/np.sqrt(2)) - err[1])/err[0])}")
    print(
        f"Bandbreitenprodukt = {10**((np.log10(1/np.sqrt(2)) - err[1])/err[0]) * V_0}"
    )

    plt.plot(ln2, lin(ln2, *params), label="Fit: abfallender Bereich", color="C2")
    # plt.semilogx()
    # plt.semilogy()
    plt.legend(loc=0, framealpha=0.4)
    plt.ylim(0, 1.2 if i == 2 else 2.5)
    plt.xlabel(r"$\log_{10}$(f/kHz)")
    plt.ylabel(r"$\log_{10}$(V)")
    plt.grid(True)
    plt.tight_layout()
    plt.title(f"Messreihe {i+1}")
    plt.tight_layout()
    plt.savefig(f"plots/linearverstaerker_{i+1}.pdf")
    plt.clf()

    fig, ax = plt.subplots(1)
    ax.plot(f, delta_t * 10 ** (-6) * 2 * f * 10 ** 3, "k+", label="Messwerte")
    plt.xlabel("f/kHz")
    ax.yaxis.set_major_formatter(plt.FormatStrFormatter("%g $\pi$"))
    ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(base=1.0))
    plt.ylabel(r"$\Delta \phi$/rad")
    plt.ylim(0.8, 2.2)
    plt.semilogx()
    plt.legend(loc=0)
    plt.grid(True)
    plt.title(f"Messreihe {i+1}")
    plt.tight_layout()
    plt.savefig(f"plots/linearverstaerker_phase_{i + 1}.pdf")
    plt.clf()

print("\n Integrator: \n")

f, U_a = np.genfromtxt("data/integrator.txt", delimiter=",", unpack=True)
log_f = np.log(f)
log_U_a = np.log(U_a)

plt.plot(log_f[:-2], log_U_a[:-2], "k+", label="Messwerte")
plt.plot(
    log_f[-2:], log_U_a[-2:], "+", c="gray", label="Nicht verwendete Messwerte",
)
params, cov = curve_fit(lin, log_f[:-2], log_U_a[:-2])
unc = np.sqrt(np.diag(cov))
err = unp.uarray(params, unc)
print(f"Fit-Parameter für den Integrator: A = {err[0]} und B = {err[1]}")

ln = np.linspace(-6, 0)
plt.plot(ln, lin(ln, *params), label="Fit")
plt.xlabel(r"$\log_{10}$(f/kHz)")
plt.ylabel(r"$\log_{10}$($U_{\mathrm{a}}$/V)")
plt.legend(loc=0)
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/integrator.pdf")
plt.clf()

print("\n Differenzierer: \n")
f, U_a = np.genfromtxt("data/differenzierer.txt", delimiter=",", unpack=True)
log_f = np.log(f)
log_U_a = np.log(U_a)

plt.plot(log_f, log_U_a, "k+", label="Messwerte")
params, cov = curve_fit(lin, log_f, log_U_a)
unc = np.sqrt(np.diag(cov))
err = unp.uarray(params, unc)
print(f"Fit-Parameter für den Integrator: A = {err[0]} und B = {err[1]}")

ln = np.linspace(-4, 1)
plt.plot(ln, lin(ln, *params), label="Fit")
plt.xlabel(r"$\log_{10}$(f/kHz)", fontsize=15)
plt.ylabel(r"$\log_{10}$($U_{\mathrm{a}}$/V)", fontsize=15)
plt.legend(loc=0)
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/differenzierer.pdf")
plt.clf()
