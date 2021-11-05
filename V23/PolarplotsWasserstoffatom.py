import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

f = np.genfromtxt(r"data/2_3/0.dat", usecols=(0))
amplitude  = np.empty((19,len(f)))
amplitude[0][:] = np.genfromtxt(r"data/2_3/0.dat", usecols=(1))
amplitude[1][:] = np.genfromtxt(r"data/2_3/10.dat", usecols=(1))
amplitude[2][:] = np.genfromtxt(r"data/2_3/20.dat", usecols=(1))
amplitude[3][:] = np.genfromtxt(r"data/2_3/30.dat", usecols=(1))
amplitude[4][:] = np.genfromtxt(r"data/2_3/40.dat", usecols=(1))
amplitude[5][:] = np.genfromtxt(r"data/2_3/50.dat", usecols=(1))
amplitude[6][:] = np.genfromtxt(r"data/2_3/60.dat", usecols=(1))
amplitude[7][:] = np.genfromtxt(r"data/2_3/70.dat", usecols=(1))
amplitude[8][:] = np.genfromtxt(r"data/2_3/80.dat", usecols=(1))
amplitude[9][:] = np.genfromtxt(r"data/2_3/90.dat", usecols=(1))
amplitude[10][:] = np.genfromtxt(r"data/2_3/100.dat", usecols=(1))
amplitude[11][:] = np.genfromtxt(r"data/2_3/110.dat", usecols=(1))
amplitude[12][:] = np.genfromtxt(r"data/2_3/120.dat", usecols=(1))
amplitude[13][:] = np.genfromtxt(r"data/2_3/130.dat", usecols=(1))
amplitude[14][:] = np.genfromtxt(r"data/2_3/140.dat", usecols=(1))
amplitude[15][:] = np.genfromtxt(r"data/2_3/150.dat", usecols=(1))
amplitude[16][:] = np.genfromtxt(r"data/2_3/160.dat", usecols=(1))
amplitude[17][:] = np.genfromtxt(r"data/2_3/170.dat", usecols=(1))
amplitude[18][:] = np.genfromtxt(r"data/2_3/180.dat", usecols=(1))

# $00 Hz
amp400Hz = np.empty(19,dtype = np.double)
for k in range(0,19):
 amp400Hz[k] = amplitude[k][61]*1/2*np.sqrt(1/np.pi) /amplitude[18][61]  #auf den ersten Wert normiert

def Y_0(w):
 t = np.empty(len(w))
 t[:] = abs(1/2*np.sqrt(1/np.pi))
 return t

#print(amp400Hz)
peak = find_peaks(amplitude[0],height = 10)
#print(peak)

alpha_grad = np.array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180])
alpha = alpha_grad*2*np.pi/360

theta = np.arccos(np.cos(alpha)/2 - 1/2)

theta_th = np.linspace(0,360,num = 200)
theta_th = theta_th * 2*np.pi/360

ampY_0 = Y_0(theta_th)

plt.polar(theta,amp400Hz,"kx",label="gemessene Druckamplitude")
plt.plot(theta_th,ampY_0,label=r"Kugelflächenfunktion $Y_0^0$ ")
plt.legend()
plt.savefig("build/Polar_04.pdf")
#plt.show()

#####2,3 kHz
plt.clf()
amp2300Hz = np.empty(19,dtype = np.double)
for k in range(0,19):
 amp2300Hz[k] = amplitude[k][438]*np.sqrt(3/(4*np.pi))/(amplitude[18][438])  #auf den ersten Wert normiert

def Y_1(w):
 return abs(np.sqrt(3/(4*np.pi))*np.cos(w))

#print(amp2300Hz)
peak = find_peaks(amplitude[0],height = 3)
#print(peak)

alpha_grad = np.array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180])
alpha = alpha_grad*2*np.pi/360

theta = np.arccos(np.cos(alpha)/2 - 1/2)

theta_th = np.linspace(0,360,num = 200)
theta_th = theta_th * 2*np.pi/360



ampY_1 = Y_1(theta_th)

plt.polar(theta,amp2300Hz,"kx",label="gemessene und normierte Druckamplitude")
plt.plot(theta_th,ampY_1,label=r"Kugelflächenfunktion $Y_1^0$")
plt.legend()
plt.savefig("build/Polar_23.pdf")

plt.clf()
amp3700Hz = np.empty(19,dtype = np.double)
for k in range(0,19):
 amp3700Hz[k] = amplitude[k][715]*2*np.sqrt(5/(16*np.pi))/(amplitude[18][715])  #auf den ersten Wert normiert


def Y_2(w):
 return abs(np.sqrt(5/(16*np.pi))*(3*np.cos(w)**2-1))

#print(amp3700Hz)
peak = find_peaks(amplitude[0],height = 10)
#print(peak)

alpha_grad = np.array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180])
alpha = alpha_grad*2*np.pi/360

theta = np.arccos(1/2*np.cos(alpha)-1/2)

theta_th = np.linspace(0,360,num = 200)
theta_th = theta_th * 2*np.pi/360

ampY_2 = Y_2(theta_th)


plt.polar(theta,amp3700Hz,"kx",label="gemessene und normierte Druckamplitude")
plt.plot(theta_th,ampY_2,label=r"Kugelflächenfunktion $Y_2^0$ ")
plt.legend()
plt.savefig(r"build/Polar_37.pdf")


plt.clf()
##7,4 kHz
amp7400Hz = np.empty(19,dtype = np.double)
for k in range(0,19):
 amp7400Hz[k] = amplitude[k][1458]*8/16 * np.sqrt(11/np.pi)/amplitude[18][1458]  #auf den ersten Wert normiert

def Y_5(w):
 return abs(1/16 * np.sqrt(11/np.pi)*(63*np.cos(w)**5 - 70*np.cos(w)**3 + 15 * np.cos(w)))

#print(amp7400Hz)
peak = find_peaks(amplitude[0],height = 10)
#print(peak)

alpha_grad = np.array([0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180])
alpha = alpha_grad*2*np.pi/360

theta = np.arccos(np.cos(alpha)/2 - 1/2)

theta_th = np.linspace(0,360,num = 200)
theta_th = theta_th * 2*np.pi/360

ampY_5 = Y_5(theta_th)

plt.polar(theta,amp7400Hz,"kx",label="gemessene und normierte Druckamplitude")
plt.plot(theta_th,ampY_5,label=r"Kugelflächenfunktion $Y_5^0$")
plt.legend()
plt.savefig(r"build/Polar_74.pdf")