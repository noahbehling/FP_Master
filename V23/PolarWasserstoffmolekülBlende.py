import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

f = np.genfromtxt(r"data/3_2/0.dat", usecols=(0))
amplitude  = np.empty((19,len(f)))
amplitude[0][:] =  np.genfromtxt(r"data/3_2/0.dat", usecols=(1))
amplitude[1][:] =  np.genfromtxt(r"data/3_2/10.dat", usecols=(1))
amplitude[2][:] =  np.genfromtxt(r"data/3_2/20.dat", usecols=(1))
amplitude[3][:] =  np.genfromtxt(r"data/3_2/30.dat", usecols=(1))
amplitude[4][:] =  np.genfromtxt(r"data/3_2/40.dat", usecols=(1))
amplitude[5][:] =  np.genfromtxt(r"data/3_2/50.dat", usecols=(1))
amplitude[6][:] =  np.genfromtxt(r"data/3_2/60.dat", usecols=(1))
amplitude[7][:] =  np.genfromtxt(r"data/3_2/70.dat", usecols=(1))
amplitude[8][:] =  np.genfromtxt(r"data/3_2/80.dat", usecols=(1))
amplitude[9][:] =  np.genfromtxt(r"data/3_2/90.dat", usecols=(1))
amplitude[10][:] = np.genfromtxt(r"data/3_2/100.dat", usecols=(1))
amplitude[11][:] = np.genfromtxt(r"data/3_2/110.dat", usecols=(1))
amplitude[12][:] = np.genfromtxt(r"data/3_2/120.dat", usecols=(1))
amplitude[13][:] = np.genfromtxt(r"data/3_2/130.dat", usecols=(1))
amplitude[14][:] = np.genfromtxt(r"data/3_2/140.dat", usecols=(1))
amplitude[15][:] = np.genfromtxt(r"data/3_2/150.dat", usecols=(1))
amplitude[16][:] = np.genfromtxt(r"data/3_2/160.dat", usecols=(1))
amplitude[17][:] = np.genfromtxt(r"data/3_2/170.dat", usecols=(1))
amplitude[18][:] = np.genfromtxt(r"data/3_2/180.dat", usecols=(1))



#peak = []
#for k in range(0,19):
#    peak = (find_peaks(amplitude[k],height = 1.5))
#    print(peak)

amp1 = np.empty(19,dtype = np.double)
for k in range(0,19):
 amp1[k] = amplitude[k][82]/amplitude[0][82]

amp2 = np.empty(19,dtype = np.double)
for k in range(0,19):
 amp2[k] = amplitude[k][86]/amplitude[0][86]

amp3 = np.empty(19,dtype = np.double)
for k in range(0,19):
 amp3[k] = amplitude[k][200]/amplitude[0][200]
#
##def Y_1_0(w):
## t = np.empty(len(w))
## t[:] = abs(np.sqrt(3/(4*np.pi)) * np.cos(np.pi/4))
## return t
#
##def Y_1_1(w):
## return np.abs(np.real(np.sqrt(3/(4*np.pi))*np.sin(np.pi/4)*np.exp(w*1j)))
#
#
alpha_grad = np.array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180])
alpha = alpha_grad*2*np.pi/360
#
theta = alpha     #np.arccos(np.cos(alpha)/2 - 1/2)
#
#
#
alpha_th = np.linspace(0,180,num = 100)
alpha_th = alpha_th * 2*np.pi/360
#
print(f[82],f[86],f[200],)
#
#
##ampY_1_0 = Y_1_0(alpha_th)
##ampY_1_1 = Y_1_1(alpha_th)
#
plt.figure()
plt.polar(theta,amp1,"kx",label="gemessene Druckamplitude")
plt.legend()
plt.savefig('build/blende_polar1.pdf')
plt.clf()
plt.figure()
plt.polar(theta,amp2,"kx",label="gemessene Druckamplitude")
plt.legend()
plt.savefig('build/blende_polar2.pdf')
plt.figure()
plt.polar(theta,amp3,"kx",label="gemessene Druckamplitude")
plt.legend()
plt.savefig('build/blende_polar3.pdf')
#ax1.plot(alpha_th,ampY_1_0,label=r"$Y_1^0$ Kugelflächenfunktion")
#ax1.set_rlim(0,0.4)

#
#fig2 = plt.figure()
#ax2 = fig2.add_subplot(111,polar="true")
#ax2.plot(theta,amp2,"kx",label="gemessene und normierte Druckamplitude")
##ax2.plot(alpha_th,ampY_1_1,label=r"$Y_1^1$ Kugelflächenfunktion")
#
#fig3 = plt.figure()
#ax3 = fig3.add_subplot(111,polar="true")
#ax3.plot(theta,amp3,"kx",label="gemessene und normierte Druckamplitude")
##ax3.plot(alpha_th,ampY_1_1,label=r"$Y_1^1$ Kugelflächenfunktion")
#
#fig1.legend()
#fig2.legend()
#fig3.legend()
#
##fig1.savefig(r"C:\Users\D-dam\Documents\TU Dortmund\Physik\Fortgeschrittenen Praktikum\V23\Auswertung\WasserstoffmolekülPol_16mmBlendeResonanz1.png")
##fig2.savefig(r"C:\Users\D-dam\Documents\TU Dortmund\Physik\Fortgeschrittenen Praktikum\V23\Auswertung\WasserstoffmolekülPol_16mmBlendeResonanz2.png")
##fig3.savefig(r"C:\Users\D-dam\Documents\TU Dortmund\Physik\Fortgeschrittenen Praktikum\V23\Auswertung\WasserstoffmolekülPol_16mmBlendeResonanz3.png")
#plt.show()
#