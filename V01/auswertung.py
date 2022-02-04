import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp


def lin(x, a, b):
    return a * x + b


def N(x, N0, a, b):
    return N0 * np.exp(-a * x) + b


N_start = unp.uarray(3648888, np.sqrt(3648888))
N_stop = 17157
N_channel = 8128 - 65
zeit = 180899
T_sz = 11e-6

print("Verzögerung:")
del_t, counts = np.genfromtxt("data/verzoegerung_comma.txt", delimiter=",", unpack=True)

del_t = np.array(del_t)
counts = np.array(counts)
counts = unp.uarray(counts, np.sqrt(counts))

ma_left = np.logical_and(-10 <= del_t, del_t <= -1.5)
ma_right = np.logical_and(1.5 <= del_t, del_t <= 10)

max_count = np.sum(counts[(-1.5 <= del_t) & (1.5 >= del_t)]) / np.count_nonzero(
    (-1.5 <= del_t) & (1.5 >= del_t)
)
plt.axhline(
    unp.nominal_values(max_count) / 2,
    -15,
    15,
    color="k",
    linestyle="--",
    linewidth=0.8,
    alpha=0.7,
)
print("Gemittelte maximale Anzahl an Counts", max_count)


params_left, cov_left = curve_fit(
    lin, del_t[ma_left], unp.nominal_values(counts)[ma_left]
)
unc = np.sqrt(np.diag(cov_left))
err = unp.uarray(params_left, unc)
t_half_left = (max_count / 2 - err[1]) / err[0]
print("Fit-Parameter linke Flanke: ", err, "\n Werte half maximum links: ", t_half_left)

ln = np.linspace(-10, -1, 10000)
plt.plot(ln, lin(ln, *params_left), "C0", label="Fits der Flanken")
plt.axvline(
    unp.nominal_values(t_half_left),
    0,
    360,
    color="k",
    linestyle="--",
    linewidth=0.8,
    alpha=0.7,
)

params_right, cov_right = curve_fit(
    lin, del_t[ma_right], unp.nominal_values(counts)[ma_right]
)
unc = np.sqrt(np.diag(cov_right))
err = unp.uarray(params_right, unc)
t_half_right = (max_count / 2 - err[1]) / err[0]
plt.axvline(
    unp.nominal_values(t_half_right),
    0,
    360,
    color="k",
    linestyle="--",
    label="Halbwertsbreite",
    linewidth=0.8,
    alpha=0.7,
)
print(
    "Fit-Parameter rechte Flanke: ", err, "\n Werte half maximum rechts: ", t_half_right
)
print("Halbwertsbreite: ", t_half_right - t_half_left, " ns")


ln = np.linspace(1, 10, 10000)
plt.plot(ln, lin(ln, *params_right), "C0")

plt.errorbar(
    del_t,
    unp.nominal_values(counts),
    xerr=None,
    yerr=unp.std_devs(counts),
    fmt="o",
    c="C1",
    ms=5,
    lw=1,
    label="Daten",
    zorder=0,
)
plt.xlabel(r"$\Delta t$/ns")
plt.ylabel(r"$N$")
plt.legend(loc=0)
plt.tight_layout()
plt.savefig("plots/HWB.pdf")

print("\n Untergrundratenbestimmung:")

f = N_start / zeit
print("Myonenfrequenz: ", f, " pro Sekunde")
N_fehl = N_start * T_sz * f * np.e ** (f * T_sz)
print("Untergrundrate: ", N_fehl)
N_fehl_norm = N_fehl / N_channel
print("Untergrundrate pro Channel: ", N_fehl_norm)


def N_fit(x, N0, a):
    return N(x, N0, a, unp.nominal_values(N_fehl_norm))


print("\n MCA Kalibrierung:")

del_t, channel1, counts1, channel2, counts2, channel3, counts3 = np.genfromtxt(
    "data/puls_comma.txt", delimiter=",", unpack=True
)
channel_norm = (channel1 * counts1 + channel2 * counts2 + channel3 * counts3) / (
    counts1 + counts2 + counts3
)

ln = np.linspace(0, 8200, 10000)
params, cov = curve_fit(lin, channel_norm, del_t)
unc = np.sqrt(np.diag(cov))
err = unp.uarray(params, unc)

plt.clf()
plt.plot(ln, lin(ln, *params), label="Linearer Fit")
print(f"a = {err[0]}, b = {err[1]}")
plt.plot(channel_norm, del_t, "k+", label="Messwerte")
plt.xlabel("Kanalnummer")
plt.ylabel(r"$\Delta t_{Puls}$" + r"/$\mu$s")
plt.legend(loc=0)
plt.tight_layout()
plt.savefig("plots/Kalibrierung.pdf")

print("\n Lebensdauer Bestimmung:")
channels, counts = np.genfromtxt("data/channel_daten.txt", delimiter=",", unpack=True)
counts = counts.astype("int32")
channels = channels.astype("int32")
N_VKA = np.sum(counts)
print("Signale im VKA: ", N_VKA)


t = err[0] * channels + err[1]


plt.clf()
params, cov = curve_fit(N_fit, unp.nominal_values(t)[65:], counts[65:])
unc = np.sqrt(np.diag(cov))
err = unp.uarray(params, unc)
ln = np.linspace(0, 12, 100000)
print(*err)


plt.errorbar(
    unp.nominal_values(t),
    counts,
    xerr=None,
    yerr=np.sqrt(counts),
    fmt="+",
    elinewidth=0.2,
    c="k",
    label="Messwerte",
    zorder=0,
)

# plt.plot(
#     unp.nominal_values(t)[:65],
#     counts[:65],
#     color="lightgray",
#     marker="+",
#     label="Nicht einbezogene Messwerte",
# )
plt.plot(ln, N_fit(ln, *params), c="r", label="Exponentieller Fit")
plt.legend(loc=0)
plt.xlabel(r"t/$\mu$s")
plt.ylabel(r"$N$")
plt.tight_layout()
plt.savefig("plots/lebensdauer.pdf")
print(f"Lebensdauer = {1 / err[1]}")
