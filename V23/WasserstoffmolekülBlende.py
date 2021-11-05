import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

def gerade(x,m,n):
 return m*x+n

f = np.genfromtxt(r"data/3_1/ohne.dat", usecols=(0))
amplitude  = np.empty((5,len(f)))
amplitude[0][:] = np.genfromtxt(r"data/3_1/ohne.dat", usecols=(1))
amplitude[1][:] = np.genfromtxt(r"data/3_1/5mm.dat", usecols=(1))
amplitude[2][:] = np.genfromtxt(r"data/3_1/10mm.dat", usecols=(1))
amplitude[3][:] = np.genfromtxt(r"data/3_1/15mm.dat", usecols=(1))
amplitude[4][:] = np.genfromtxt(r"data/3_1/20mm.dat", usecols=(1))


 



peak_kb = find_peaks(amplitude[0],height = 2)
peak_5mm = find_peaks(amplitude[1],height = 2)
peak_10mm = find_peaks(amplitude[2],height = 1.8)
peak_15mm = find_peaks(amplitude[3],height = 2)
peak_20mm = find_peaks(amplitude[4],height = 1.8)
print(peak_kb[0],peak_5mm[0],peak_10mm[0],peak_15mm[0], peak_20mm[0])
Blendendurchmesser = np.array([0,5,10,15, 20])




plt.plot([0,0,0],f[peak_kb[0]],"kx",label="Resonanzen ohne Blende")
plt.plot([5,5],f[peak_5mm[0]],"bx",label="Resonanzen 5mm Blende")
plt.plot([10,10,10],f[peak_10mm[0]],"gx",label="Resonanzen 10mm Blende")
plt.plot([15,15,15],f[peak_15mm[0]],"rx",label="Resonanzen 15mm Blende")
plt.plot([20,20,20],f[peak_20mm[0]],"x",label="Resonanzen 20mm Blende")
plt.xlabel(r"Blendendurchmesser $d$/mm")
plt.ylabel(r"Frequenz der Peaks $f$/ Hz")
plt.legend()
plt.savefig(r"build/frequenz_blende.pdf")

