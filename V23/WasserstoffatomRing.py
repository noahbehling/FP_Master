import matplotlib.pyplot as plt
import numpy as np
from uncertainties import ufloat
import uncertainties.unumpy as unp 
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

def gerade(x,m,n):
 return m*x+n

f = np.genfromtxt(r"data/2_4/3mm.dat", usecols=(0))
amplitude  = np.empty((3,len(f)))
#amplitude[0][:] = np.genfromtxt(r"data/2_4/3mm.dat", usecols=(1))
amplitude[0][:] = np.genfromtxt(r"data/2_4/3mm.dat", usecols=(1))
amplitude[1][:] = np.genfromtxt(r"data/2_4/6mm.dat", usecols=(1))
amplitude[2][:] = np.genfromtxt(r"data/2_4/9mm.dat", usecols=(1))

#peak_kr = find_peaks(amplitude[0],height = 10)
peak_3mm = find_peaks(amplitude[0],height = 5)
peak_6mm = find_peaks(amplitude[1],height = 5)
peak_9mm = find_peaks(amplitude[2],height = 3.5)

print(peak_3mm[0][0])
f_diff = np.empty(3,dtype = np.double)

#f_diff[0] = abs(f[peak_kr[0][0]] - f[peak_kr[0][1]])
f_diff[0] = abs(f[peak_3mm[0][0]] - f[peak_3mm[0][1]])
f_diff[1] = abs(f[peak_6mm[0][0]] - f[peak_6mm[0][1]])
f_diff[2] = abs(f[peak_9mm[0][0]] - f[peak_9mm[0][1]])

Ringdicke = np.array([3,6,9])

popt, pcov = curve_fit(gerade,Ringdicke,f_diff)
a = ufloat(popt[0], np.sqrt(pcov[0][0]))
b = ufloat(popt[1], np.sqrt(pcov[1][1]))

print(a)
print(b)
#print(f[peak_kr[0][0]])
#print(peak_kr[0])


plt.plot(Ringdicke,f_diff,"kx",label="Messwerte der Resonanzen")
plt.plot(Ringdicke,gerade(Ringdicke, *popt),label="Ausgleichsgerade")
plt.xlabel(r"Ringdurchmesser $d$/mm")
plt.ylabel(r"Frequenzdifferenz der Peaks $\Delta f$/$\frac{1}{\text{s}}$")
plt.legend()
plt.grid(alpha=0.4)
plt.savefig(r"build/ausgleichsgerade.pdf")
