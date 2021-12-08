import numpy as np
import matplotlib.pyplot as plt



# Reflektivitätsscan:
angle_refl, intensity_refl = np.genfromtxt('data/scan1.UXD', unpack=True)
# Diffuser Scan
angle_diff, intensity_diff = np.genfromtxt('data/scn2.UXD', unpack=True)
    
#################
## Berechnungen
#################

# Reflektivität: R=I/5I_0; 5 durch 5s Messzeit
I_0 = 9.7272e+05 #Werte jeweils aus vorherigen Berechnungen
R_refl = intensity_refl / (5 *I_0)
R_diff = intensity_diff / (5 *I_0)

R = R_refl - R_diff

# Geometriefaktor anwenden; für Vergleich mir Parratt
a_g = 0.6875658564108561
d_0 = 0.24
D = 20 
s = np.size(R)
G = np.ones(s)
G[angle_refl < a_g] = D/d_0 * np.sin( np.deg2rad(angle_refl[angle_refl < a_g]) )
R_G = R * G

# Parratt-Algorithmus

# a_i Einfallswinkel, Brechuingsindex: n1 Luft, n2 Polystyrol, n3 Silicium
# Rauigkeiten: sigma1 Polystyrol, sigma2 Silizium
# Schichtdicken: d1=0 (Luft), d2: Dicke der Schicht Polystyrol

n1 = 1.
d1 = 0.
wellenlaenge = 1.54e-10 
# Ermittelte Werte durch Anpassung des Parratt
delta2 = 0.6*10**(-6)
delta3 = 6.0*10**(-6)
sigma1 = 5.5*10**(-10) # m
sigma2 = 6.45*10**(-10) # m
d2 = 8.6*10**(-8) # m

def parratt(a_i,delta2,delta3,sigma1,sigma2,d2):
    n2 = 1. - delta2
    n3 = 1. - delta3
    a_i = np.deg2rad(a_i)
    k = 2 * np.pi/wellenlaenge
    kd1 = k * np.sqrt(np.abs(n1**2 - np.cos(a_i)**2))
    kd2 = k * np.sqrt(np.abs(n2**2 - np.cos(a_i)**2))
    kd3 = k * np.sqrt(np.abs(n3**2 - np.cos(a_i)**2))

    r12 = (kd1 - kd2) / (kd1 + kd2) * np.exp(-2 * kd1 * kd2 * sigma1**2)
    r23 = (kd2 - kd3) / (kd2 + kd3) * np.exp(-2 * kd2 * kd3 * sigma2**2)

    x2 = np.exp(-2j * kd2 * d2) * r23
    x1 = (r12 + x2) / (1 + r12 * x2)

    return  np.abs(x1)**2

params = [delta2,delta3,sigma1,sigma2,d2]
print(params)
parratt_curve = parratt(angle_refl, *params)

# Kritischer Winkel
a_Poly = np.rad2deg(np.sqrt(2*delta2))
a_Si = np.rad2deg(np.sqrt(2*delta3))

print("Kritischer Winkel Polystyrol: ", a_Poly)
print("Kritischer Winkel Silicium: ", a_Si)

parratt_curve[angle_refl <= a_Si] = np.nan

#Plot

plt.axvline(a_Poly, linewidth=0.5, color='green', label='a Polystyrol')
plt.axvline(a_Si, linewidth=0.5, color='pink', label='a Silizium')
plt.plot(angle_refl, parratt_curve, '-', label='Parratt-Kurve')
plt.plot(angle_refl, R_G, '-', label='gem. Reflektivität mit Korrektur durch Geometriefaktor')
plt.xlabel('a / °')
plt.ylabel('R')
plt.yscale('log')
plt.grid(alpha=0.4)
plt.legend()
plt.savefig('parrott.pdf')