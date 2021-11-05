import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

f = np.genfromtxt(r"data/2_5/0.dat", usecols=(0))
amplitude  = np.empty((19,len(f)))
amplitude[0][:] = np.genfromtxt(r"data/2_5/0.dat", usecols=(1))
amplitude[1][:] = np.genfromtxt(r"data/2_5/10.dat", usecols=(1))
amplitude[2][:] = np.genfromtxt(r"data/2_5/20.dat", usecols=(1))
amplitude[3][:] = np.genfromtxt(r"data/2_5/30.dat", usecols=(1))
amplitude[4][:] = np.genfromtxt(r"data/2_5/40.dat", usecols=(1))
amplitude[5][:] = np.genfromtxt(r"data/2_5/50.dat", usecols=(1))
amplitude[6][:] = np.genfromtxt(r"data/2_5/60.dat", usecols=(1))
amplitude[7][:] = np.genfromtxt(r"data/2_5/70.dat", usecols=(1))
amplitude[8][:] = np.genfromtxt(r"data/2_5/80.dat", usecols=(1))
amplitude[9][:] = np.genfromtxt(r"data/2_5/90.dat", usecols=(1))
amplitude[10][:] = np.genfromtxt(r"data/2_5/100.dat", usecols=(1))
amplitude[11][:] = np.genfromtxt(r"data/2_5/110.dat", usecols=(1))
amplitude[12][:] = np.genfromtxt(r"data/2_5/120.dat", usecols=(1))
amplitude[13][:] = np.genfromtxt(r"data/2_5/130.dat", usecols=(1))
amplitude[14][:] = np.genfromtxt(r"data/2_5/140.dat", usecols=(1))
amplitude[15][:] = np.genfromtxt(r"data/2_5/150.dat", usecols=(1))
amplitude[16][:] = np.genfromtxt(r"data/2_5/160.dat", usecols=(1))
amplitude[17][:] = np.genfromtxt(r"data/2_5/170.dat", usecols=(1))
amplitude[18][:] = np.genfromtxt(r"data/2_5/180.dat", usecols=(1))


peak = []
for k in range(0,18):
 peak.append(find_peaks(amplitude[k],height = 3.0))


amp1 = np.empty(18,dtype = np.double)
for k in range(0,18):
 print(k)
 amp1[k] = amplitude[k][280]/amplitude[0][280] * np.sqrt(3/(4*np.pi)) *np.cos(np.pi/4)

amp2 = np.empty(18,dtype = np.double)
for k in range(0,18):
 amp2[k] = amplitude[k][450]/amplitude[0][450] * np.sqrt(3/(4*np.pi))*np.sin(np.pi/4)

def Y_1_0(w):
 t = np.empty(len(w))
 t[:] = abs(np.sqrt(3/(4*np.pi)) * np.cos(np.pi/4))
 return t

def Y_1_1(w):
 return np.abs(np.real(np.sqrt(3/(4*np.pi))*np.sin(np.pi/4)*np.exp(w*1j)))


alpha_grad = np.array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170])
alpha = alpha_grad*2*np.pi/360

#theta = np.arccos(np.cos(alpha)/2 - 1/2)

alpha_th = np.linspace(0,360,num = 200)
alpha_th = alpha_th * 2*np.pi/360

ampY_1_0 = Y_1_0(alpha_th)
ampY_1_1 = Y_1_1(alpha_th)



#ax1.set_rlim(0,0.4)
plt.polar(alpha,amp1,"kx",label="gemessene Druckamplitude")
plt.plot(alpha_th,ampY_1_0,label=r"Kugelflächenfunktion $Y_1^0$")
plt.legend()
plt.savefig(r"build/resonanz1.pdf")
plt.clf()

plt.polar(alpha,amp2,"kx",label="gemessene Druckamplitude")
plt.plot(alpha_th,ampY_1_1,label=r"Kugelflächenfunktion $Y_1^1$")
plt.legend()
plt.savefig(r"build/resonanz2.pdf")
#
