import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.optimize import curve_fit, root
from scipy.signal import find_peaks
from scipy.stats import sem
from uncertainties import ufloat

##################
## Daten Einlesen
##################
# Ergebnisse JSON Datei einlesen (am Anfang)
json_file_path = 'data/Ergebnisse.json'
try:
    with open(json_file_path,'r') as json_file:
        Ergebnisse = json.load(json_file)
except FileNotFoundError as err:
    Ergebnisse = {}

if not 'Messung' in Ergebnisse:
    Ergebnisse['Messung'] = dict()

# Reflektivitätsscan:
a_refl, I_refl = np.genfromtxt('data/scan1.UXD', unpack=True)
# Diffuser Scan
a_diff, I_diff = np.genfromtxt('data/scn2.UXD', unpack=True)

if (a_refl!=a_diff).any():
    print('Der Diffuse Scan und der Reflektivitätsscan passen nicht zueinander!')

# Winkel sollten gleich sein
a = a_refl
    
# Anfang und Ende abschneiden
a_min = 0.01 
a_max = 1.6
mask = (a >= a_min) & (a <= a_max)
a = a[mask]
I_refl = I_refl[mask]
I_diff = I_diff[mask]

#################
## Berechnungen
#################
# Eingehende Intensität als das Maximum vom Detektorscan
# aber mit 5 multipliziert weil nun statt 1s 5s pro Winkel gemessen wurden
I_0 = float(Ergebnisse['Detektorscan']['I_max_gauss']) * 5

# Reflektivität: R=I_r/I_0
R_refl = I_refl / I_0
R_diff = I_diff / I_0

# diffusen Scan abziehen
R = R_refl - R_diff

# Geometriewinkel
a_g = float(Ergebnisse['Rockingscan']['alpha_g[degree]'])

# Strahlbreite
d_0 = float(Ergebnisse['Z-Scan']['d_0[mm]'])
D = 20 # mm

# Geometriefaktor
G = np.ones_like(R)
G[a < a_g] = D/d_0 * np.sin( np.deg2rad(a[a < a_g]) )

# um Geometriefaktor korrigieren
R_G = R * G

# Ideale Fresnelreflektivität
a_c_Si = 0.223
R_ideal = (a_c_Si / (2 * a))**4

## Peaks finden
# Curve Fit für find_peaks
peaks_mask = (a>=0.3) & (a<=1.19)
def f(x,b,c):
    return b*x+c

params, pcov = curve_fit(f,a[peaks_mask],np.log(R_G[peaks_mask]))
R_fit = np.exp(f(a[peaks_mask],*params))

# Minima der Kissig-Oszillation finden
i_peaks, peak_props = find_peaks(-(R_G[peaks_mask]-R_fit), distance=7)
i_peaks += np.where(peaks_mask)[0][0]

# Schichtdicke bestimmen
lambda_ = 1.54*10**(-10) # m

delta_a = np.diff(np.deg2rad(a[i_peaks]))
delta_a_mean = ufloat(np.mean(delta_a),sem(delta_a))

d = lambda_ / (2*delta_a_mean)

Ergebnisse['Messung']['delta_a_mean[degree]'] = f'{delta_a_mean:.2u}'
Ergebnisse['Messung']['d[m]'] = f'{d:.2u}'

## Parrat Algorithmus

# Speichere R_G und a_i für den interaktiven Plot
np.savetxt('build/R_G.csv', list(zip(a, R_G)), header='a_i,R_G', fmt='%.4f,%.10e')

# a_i Einfallswinkel
# n sind Brechungsindizes
# n1 Luft, n2 Schicht, n3 Substrat
# sigma sind Raugigkeiten
# sigma1 Schicht, sigma2 Substrat
# z1=0, z2 Schichtdicke
# k=2pi/lambda Betrag des Wellenwektors
# Konstanten:
n1 = 1.
z1 = 0.
k = 2*np.pi/lambda_ 

# Werte durch Anpassung, sodass R_G und R_parr gut passen
delta2 = 0.5*10**(-6)
delta3 = 6.75*10**(-6)
sigma1 = 8.0*10**(-10) # m
sigma2 = 6.3*10**(-10) # m
z2 = 8.55*10**(-8) # m

def parrat_rau(a_i,delta2,delta3,sigma1,sigma2,z2):
    n2 = 1. - delta2
    n3 = 1. - delta3

    a_i = np.deg2rad(a_i)

    kz1 = k * np.sqrt(np.abs(n1**2 - np.cos(a_i)**2))
    kz2 = k * np.sqrt(np.abs(n2**2 - np.cos(a_i)**2))
    kz3 = k * np.sqrt(np.abs(n3**2 - np.cos(a_i)**2))

    r12 = (kz1 - kz2) / (kz1 + kz2) * np.exp(-2 * kz1 * kz2 * sigma1**2)
    r23 = (kz2 - kz3) / (kz2 + kz3) * np.exp(-2 * kz2 * kz3 * sigma2**2)

    x2 = np.exp(-2j * kz2 * z2) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)
    R_parr = np.abs(x1)**2

    return R_parr

params = [delta2,delta3,sigma1,sigma2,z2]

R_parr = parrat_rau(a, *params)


# Kritischer Winkel
a_c2 = np.rad2deg(np.sqrt(2*delta2))
a_c3 = np.rad2deg(np.sqrt(2*delta3))

Ergebnisse['Messung']['a_c2[degree]'] = a_c2
Ergebnisse['Messung']['a_c3[degree]'] = a_c3

R_ideal[a <= a_c3] = np.nan
R_parr[a <= a_c3] = np.nan

############
## Plotten
############
# Reflektivitäts Scan Plotten
print('Plot: Mess-Scan...')
mpl.rcParams['lines.linewidth'] = 0.9
mpl.rcParams['axes.grid.which'] = 'major'
plt.axvline(a_c2, linewidth=0.6, linestyle='dashed', color='blue', label=r'$\alpha_\text{c,PS},\alpha_\text{c,Si}$')
plt.axvline(a_c3, linewidth=0.6, linestyle='dashed', color='blue')
plt.plot(a, R_refl/10, '-', label='Reflektivitätsscan / 10')
plt.plot(a, R_diff/10, '-', label='Diffuser Scan / 10')
plt.plot(a, R/10, '-', label='Reflektivitätsscan - Diffuser Scan / 10')
plt.plot(a, R_ideal, '-',color='pink', label='Fresnelreflektivität von Si')
plt.plot(a, R_parr, '-', label='Theoriekurve (manueller Fit)')
plt.plot(a, R_G, '-', label=r'(Reflektivitätsscan - Diffuser Scan)$\cdot G$')
plt.plot(a[i_peaks], R_G[i_peaks], 'kx', label='Oszillationsminima',alpha=0.8)
# plt.plot(a[peaks_mask],R_fit, '--', label='Peaks Curve Fit')
plt.xlabel(r'$\alpha_\text{i} \:/\: \si{\degree}$')
plt.ylabel(r'$R$')
plt.yscale('log')
plt.legend(loc='upper right',prop={'size': 8})
plt.tight_layout(pad=0.15, h_pad=1.08, w_pad=1.08)
plt.savefig('build/plot_messung.pdf')
#plt.show()
plt.clf()

##########
# Ergebnisse als JSON Datei speichern (am Ende)
with open(json_file_path,'w') as json_file:
    json.dump(Ergebnisse, json_file, indent=4)