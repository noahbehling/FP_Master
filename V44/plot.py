import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, root
from uncertainties import ufloat
import uncertainties.unumpy as unp
import matplotlib as mpl
from scipy.signal import find_peaks
from scipy.stats import sem


json_file_path = "data/Ergebnisse.json"
try:
    with open(json_file_path, "r") as json_file:
        Ergebnisse = json.load(json_file)
except FileNotFoundError as err:
    Ergebnisse = {}

angle, intensity = np.genfromtxt("data/detectorscn.UXD", unpack=True)
# Gaußfunktion als Ausgleichskurve
def gauss(x, a, b, sigma, mu):
    return (
        a
        / np.sqrt(2 * np.pi * sigma ** 2)
        * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))
        + b
    )


# mask = (angle > -0.1) & (angle != 0.3)
mask = angle == angle  # dummy mask
p0 = [10 ** 6, 0, 10 ** (-2), 10 ** (-2)]
params, pcov = curve_fit(gauss, angle[mask], intensity[mask], p0=p0)

# Ausgleichskurven Parameter abspeichern
if not "Detektorscan" in Ergebnisse:
    Ergebnisse["Detektorscan"] = dict()
if not "Ausgleichsrechnung" in Ergebnisse["Detektorscan"]:
    Ergebnisse["Detektorscan"]["Ausgleichsrechnung"] = dict()
for i, name in enumerate(["a", "b", "sigma[degree]", "mu[degree]"]):
    Ergebnisse["Detektorscan"]["Ausgleichsrechnung"][
        name
    ] = f"{ufloat(params[i],np.absolute(pcov[i][i])**0.5) : .2u}"

a = ufloat(params[0], np.absolute(pcov[0][0]) ** 0.5)
b = ufloat(params[1], np.absolute(pcov[1][1]) ** 0.5)
s = ufloat(params[2], np.absolute(pcov[2][2]) ** 0.5)
m = ufloat(params[3], np.absolute(pcov[3][3]) ** 0.5)
print(a / unp.sqrt(2 * np.pi * s ** 2) * unp.exp(-((m - m) ** 2) / (2 * s ** 2)) + b)
angle_linspace = np.linspace(np.min(angle), np.max(angle), 1000)
intensity_gauss = gauss(angle_linspace, *params)

# Intensitäts Maximum
I_max = np.max(intensity_gauss)
Ergebnisse["Detektorscan"]["I_max_gemessen"] = f"{np.max(intensity) : .4e}"
Ergebnisse["Detektorscan"]["I_max_gauss"] = f"{I_max : .4e}"

# Halbwertsbreite anhand der Ausgleichsrechnung
# Full Width Half Maximum
left_FWHM = root(lambda x: gauss(x, *params) - (I_max / 2), x0=-0.01).x[0]
right_FWHM = root(lambda x: gauss(x, *params) - (I_max / 2), x0=0.1).x[0]
uleft_FWHM = (
    a / unp.sqrt(2 * np.pi * s ** 2) * unp.exp(-((left_FWHM - m) ** 2) / (2 * s ** 2))
    + b
)
uright_FWHM = (
    a / unp.sqrt(2 * np.pi * s ** 2) * unp.exp(-((right_FWHM - m) ** 2) / (2 * s ** 2))
    + b
)
print(uright_FWHM)
FWHM = right_FWHM - left_FWHM
print(FWHM)
Ergebnisse["Detektorscan"]["Halbwertsbreite[degree]"] = f"{FWHM : .4e}"
# Wikipedia: Die Halbwertsbreite einer Normalverteilung ist das ungefähr 2,4-Fache (genau 2*sqrt(2*ln(2))) der Standardabweichung.
Ergebnisse["Detektorscan"][
    "Halbwertsbreite_gauss[degree]"
] = f"{params[2]*(2*np.sqrt(2*np.log(2))) : .4e}"
print(pcov[2][2] ** (0.5) * (2 * np.sqrt(2 * np.log(2))))
# Detektor Scan Plotten
plt.plot(angle_linspace, intensity_gauss, "r--", label="Ausgleichskurve")
# plt.plot([left_FWHM, right_FWHM], [I_max/2, I_max/2], 'g--', label='Halbwertsbreite')
plt.vlines(left_FWHM, 0, 1e6, "g", label="FWHM")
plt.vlines(right_FWHM, 0, 1e6, "g")
plt.plot(angle, intensity, "bx", label="Messdaten")
plt.xlabel(r"$\alpha \:/\: \si{\degree}$")
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.ylabel(r"$I \:/\:$ Hits / s")
plt.legend()
plt.grid(alpha=0.6)
plt.tight_layout(pad=0.15, h_pad=1.08, w_pad=1.08)
plt.savefig("build/plot_detektorscan.pdf")
plt.clf()

print("Plot: Z-Scan...")
z, intensity = np.genfromtxt("data/z-scan1.UXD", unpack=True)

# Strahlbreite Ablesen
i_d = [23, 29]
d0 = np.abs(z[i_d[0]] - z[i_d[1]])  # mm

if not "Z-Scan" in Ergebnisse:
    Ergebnisse["Z-Scan"] = dict()
Ergebnisse["Z-Scan"]["d_0[mm]"] = d0

# Z Scan Plotten
plt.axvline(z[i_d[0]], color="red", linestyle="dashed", label="Strahlgrenzen")
plt.axvline(z[i_d[1]], color="red", linestyle="dashed")
plt.plot(z, intensity, "bx", label="Messdaten")
plt.xlabel(r"$z \:/\: \si{\milli\meter}$")
plt.ylabel(r"$I \:/\:$ Hits / s")
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.legend()
plt.grid(alpha=0.4)
plt.tight_layout(pad=0.15, h_pad=1.08, w_pad=1.08)
plt.savefig("build/plot_zscan.pdf")
# plt.show()
plt.clf()


##################
## Rocking-Scan 1
##################
print("Plot: Rocking-Scan...")
angle, intensity = np.genfromtxt("data/rockingcurve.UXD", unpack=True)

# Geometriewinkel ablesen
i_g = [8, -9]
a_g = np.mean(np.abs(angle[i_g]))

D = 20  # mm

# Geometriewinkel brechnen aus Strahlbreite und Probenlänge
a_g_berechnet = np.rad2deg(np.arcsin(d0 / D))

if not "Rockingscan" in Ergebnisse:
    Ergebnisse["Rockingscan"] = dict()
Ergebnisse["Rockingscan"]["alpha_g_l[degree]"] = angle[i_g[0]]
Ergebnisse["Rockingscan"]["alpha_g_r[degree]"] = angle[i_g[1]]
Ergebnisse["Rockingscan"]["alpha_g[degree]"] = a_g
Ergebnisse["Rockingscan"]["alpha_g_berechnet[degree]"] = a_g_berechnet


# Rocking Scan Plotten
plt.axvline(angle[i_g[0]], color="red", linestyle="dashed", label="Geometriewinkel")
plt.axvline(angle[i_g[1]], color="red", linestyle="dashed")
plt.plot(angle, intensity, "bx", label="Messdaten")
plt.xlabel(r"$\alpha \:/\: \si{\degree}$")
plt.ylabel(r"$I \:/\:$ Hits / s")
plt.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
plt.legend()
plt.grid(alpha=0.4)
plt.tight_layout(pad=0.15, h_pad=1.08, w_pad=1.08)
plt.savefig("build/plot_rockingscan.pdf")
# plt.show()
plt.clf()

# Reflektivitätsscan:
a_refl, I_refl = np.genfromtxt("data/scan1.UXD", unpack=True)
# Diffuser Scan
a_diff, I_diff = np.genfromtxt("data/scn2.UXD", unpack=True)

if (a_refl != a_diff).any():
    print("Der Diffuse Scan und der Reflektivitätsscan passen nicht zueinander!")

# Winkel sollten gleich sein
a = a_refl

# Anfang und Ende abschneiden
a_min = 0.0
a_max = 2.5
mask = (a >= a_min) & (a <= a_max)
a = a[mask]
I_refl = I_refl[mask]
I_diff = I_diff[mask]

#################
## Berechnungen
#################
# Eingehende Intensität als das Maximum vom Detektorscan
# aber mit 5 multipliziert weil nun statt 1s 5s pro Winkel gemessen wurden
I_0 = float(Ergebnisse["Detektorscan"]["I_max_gauss"]) * 5

# Reflektivität: R=I_r/I_0
R_refl = I_refl / I_0
R_diff = I_diff / I_0

# diffusen Scan abziehen
# R = R_refl - R_diff

# Geometriewinkel
a_g = float(Ergebnisse["Rockingscan"]["alpha_g[degree]"])

# Strahlbreite
d_0 = float(Ergebnisse["Z-Scan"]["d_0[mm]"])
D = 20  # mm
print(d_0, D)

# Geometriefaktor
# G = np.ones_like(R_refl)
# um Geometriefaktor korrigieren
# R = R_refl - R_diff
# R_refl = R_refl * G
# R_diff = R_diff * G
R = R_refl - R_diff
R_G = np.zeros(np.size(R))
for i in np.arange(np.size(a)):
    if a[i] <= a_g and a[i] > 0:
        R_G[i] = R[i] * np.sin(np.deg2rad(a_g)) / np.sin(np.deg2rad(a[i]))
    else:
        R_G[i] = R[i]


# Ideale Fresnelreflektivität
a_c_Si = 0.223
R_ideal = (a_c_Si / (2 * a)) ** 4

## Peaks finden
# Curve Fit für find_peaks
peaks_mask = (a >= 0.2) & (a < 1.0)


def f(x, b, c):
    return b * x + c


params, pcov = curve_fit(f, a[peaks_mask], np.log(R_G[peaks_mask]))
R_fit = np.exp(f(a[peaks_mask], *params))

# Minima der Kissig-Oszillation finden
i_peaks, peak_props = find_peaks(-(R_G[peaks_mask] - R_fit), distance=7)
i_peaks += np.where(peaks_mask)[0][0]

# Schichtdicke bestimmen
lambda_ = 1.54 * 10 ** (-10)  # m

delta_a = np.diff(np.deg2rad(a[i_peaks]))
delta_a_mean = ufloat(np.mean(delta_a), sem(delta_a))

d = lambda_ / (2 * delta_a_mean)
if not "Messung" in Ergebnisse:
    Ergebnisse["Messung"] = dict()
Ergebnisse["Messung"]["delta_a_mean[degree]"] = f"{delta_a_mean:.2u}"
Ergebnisse["Messung"]["d[m]"] = f"{d:.2u}"

## Parrat Algorithmus

# Speichere R_G und a_i für den interaktiven Plot
np.savetxt("build/R_G.csv", list(zip(a, R_G)), header="a_i,R_G", fmt="%.4f,%.10e")

# a_i Einfallswinkel
# n sind Brechungsindizes
# n1 Luft, n2 Schicht, n3 Substrat
# sigma sind Raugigkeiten
# sigma1 Schicht, sigma2 Substrat
# z1=0, z2 Schichtdicke
# k=2pi/lambda Betrag des Wellenwektors
# Konstanten:
n1 = 1.0
z1 = 0.0
k = 2 * np.pi / lambda_

# Werte durch Anpassung, sodass R_G und R_parr gut passen
delta2 = 0.6 * 10 ** (-6)
delta3 = 6.0 * 10 ** (-6)
sigma1 = 5.5 * 10 ** (-10)  # m
sigma2 = 6.45 * 10 ** (-10)  # m
z2 = 8.6 * 10 ** (-8)  # m
b1 = (delta2 / 40) * 1j
b2 = (delta3 / 200) * 1j


def parrat_rau(a_i, delta2, delta3, sigma1, sigma2, z2):
    n2 = 1.0 - delta2 + b1
    n3 = 1.0 - delta3 + b2

    a_i = np.deg2rad(a_i)

    kz1 = k * np.sqrt(n1 ** 2 - np.cos(a_i) ** 2)
    kz2 = k * np.sqrt(n2 ** 2 - np.cos(a_i) ** 2)
    kz3 = k * np.sqrt(n3 ** 2 - np.cos(a_i) ** 2)

    r12 = (kz1 - kz2) / (kz1 + kz2) * np.exp(-2 * kz1 * kz2 * sigma1 ** 2)
    r23 = (kz2 - kz3) / (kz2 + kz3) * np.exp(-2 * kz2 * kz3 * sigma2 ** 2)

    x2 = np.exp(-2j * kz2 * z2) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    R_parr = np.abs(x1) ** 2 
    return R_parr


params = [delta2, delta3, sigma1, sigma2, z2]

R_parr = parrat_rau(a, *params)

print(params)


# Kritischer Winkel
a_c2 = np.rad2deg(np.sqrt(2 * delta2))
a_c3 = np.rad2deg(np.sqrt(2 * delta3))

Ergebnisse["Messung"]["a_c2[degree]"] = a_c2
Ergebnisse["Messung"]["a_c3[degree]"] = a_c3

R_ideal[a <= a_c3] = np.nan
#if R_parr[a <= a_c3] > 1:
R_parr[a <= a_c3]= 1

############
## Plotten
############
# Reflektivitäts Scan Plotten
print("Plot: Part 2")
mpl.rcParams["lines.linewidth"] = 0.9
mpl.rcParams["axes.grid.which"] = "major"
plt.plot(a, R, "b-", label="gemessene Reflektivität")
# plt.plot(a, R_ideal, '-', color='black', label='theor. Fresnelreflektivität')
plt.plot(a, R_G, "r-", label=r"gem. Reflektivität mit Korrektur durch Geometriefaktor")
plt.plot(a[i_peaks], R_G[i_peaks], ".", label="Minima", alpha=0.8)
# plt.plot(a[peaks_mask],R_fit, '--', label='Peaks Curve Fit')
plt.xlabel(r"$\alpha_\text{i} \:/\: \si{\degree}$")
plt.ylabel(r"$R$")
plt.yscale("log")
plt.grid(alpha=0.4)
plt.legend(loc="upper right", prop={"size": 8})
plt.tight_layout(pad=0.15, h_pad=1.08, w_pad=1.08)
plt.savefig("build/plot_reflektivitaet.pdf")
# plt.show()
plt.clf()


print("Plot: Mess-Scan...")
mpl.rcParams["lines.linewidth"] = 0.9
mpl.rcParams["axes.grid.which"] = "major"
# plt.axvline(a_c2, linewidth=0.6, linestyle='dashed', color='blue', label=r'$\alpha_\text{c,PS},\alpha_\text{c,Si}$')
# plt.axvline(a_c3, linewidth=0.6, linestyle='dashed', color='blue')
# plt.plot(a, R_refl/10, '-', label='Reflektivitätsscan / 10')
# plt.plot(a, R_diff/10, '-', label='Diffuser Scan / 10')
# plt.plot(a, R/10, '-', label='Reflektivitätsscan - Diffuser Scan / 10')
plt.axvline(a_c2, linewidth=0.5, color="green", label=r"$\alpha_\text{c,PS}")
plt.axvline(a_c3, linewidth=0.5, color="pink", label=r"$\alpha_\text{c,Si}$")
plt.plot(a, R_ideal, "-", color="black", label="theor. Fresnelreflektivität")
plt.plot(a, R_parr, "-", label="Parratt-Kurve")
plt.plot(a, R_G, "-", label=r"gem. Reflektivität mit Korrektur durch Geometriefaktor")
#plt.plot(a, np.ones_like(a))
# plt.plot(a[i_peaks], R_G[i_peaks], 'kx', label='Oszillationsminima',alpha=0.8)
# plt.plot(a[peaks_mask],R_fit, '--', label='Peaks Curve Fit')
plt.xlabel(r"$\alpha_\text{i} \:/\: \si{\degree}$")
plt.ylabel(r"$R$")
plt.yscale("log")
plt.grid(alpha=0.4)
plt.legend(loc="upper right", prop={"size": 8})
plt.tight_layout(pad=0.15, h_pad=1.08, w_pad=1.08)
plt.savefig("build/plot_messung.pdf")
# plt.show()
plt.clf()

with open(json_file_path, "w") as json_file:
    json.dump(Ergebnisse, json_file, indent=4)
