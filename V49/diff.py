import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal

#Laden der Daten aus der Datei "echo_gradient.csv" 

#Die erste Spalte enthält die Zeiten in Sekunden, die zweite Spalte  

#den Realteil und die dritte Spalte den Imaginärteil 

data = np.loadtxt("data/diffus.csv", delimiter=",", skiprows=3, unpack= True) 

times = data[0] 

real = data[1] 

imag = data[2] 
plt.plot(times, real, label=r'Realteil')
plt.plot(times, imag, label=r'Imaginärteil')
plt.xlabel(r'$t / \,$s')
plt.ylabel(r'$U / \,$V')
plt.legend()
plt.tight_layout()
plt.grid(alpha=0.5)
plt.savefig('build/Echo.pdf')
plt.clf()

#Suchen des Echo-Maximums und alle Daten davor abschneiden 

start = np.argmax(real) 

times = times[start:]

real = real[start:] 
imag = imag[start:] 

#Phasenkorrektur - der Imaginärteil bei t=0 muss = 0 sein 

phase = np.arctan2(imag[0], real[0]) 

#Daten in komplexes Array mit Phasenkorrektur speichern 
compsignal = (real*np.cos(phase)+imag*np.sin(phase))+ (-real*np.sin(phase)+imag*np.cos(phase))*1j 
#Offsetkorrektur, ziehe den Mittelwert der letzten 512 Punkte von allen Punkten ab 

compsignal = compsignal - compsignal[-512:-1].mean() 

#Der erste Punkt einer FFT muss halbiert werden 

compsignal[0] = compsignal[0]/2.0 

#Anwenden einer Fensterfunktion (siehe z. Bsp. 

#https://de.wikipedia.org/wiki/Fensterfunktion ) 

#Hier wird eine Gaußfunktion mit sigma = 100 Hz verwendet 

apodisation = 100.0*2*np.pi 

compsignal = compsignal*np.exp(-1.0/2.0*((times-times[0])*apodisation)**2) 

#Durchführen der Fourier-Transformation 

fftdata = np.fft.fftshift(np.fft.fft(compsignal)) 
#Generieren der Frequenzachse 
freqs = np.fft.fftshift(np.fft.fftfreq(len(compsignal), times[1]-times[0])) 

#Speichern des Ergebnisses als txt 

np.savetxt("echo_gradient_fft.txt", np.array([freqs, np.real(fftdata), np.imag(fftdata)]).transpose()) 
#Erstellen eines Plots 
#[466:497]
a =np.real(fftdata).astype(float)
b =freqs.astype(float)
a = scipy.signal.peak_widths( a, [481], rel_height=0.92)
print("a", a)
h = a[0]*4.970673029128921598e+02
print("h", h)
plt.plot(freqs[450:510], np.real(fftdata[450:510]), label=r'Fouriertransformation')
plt.plot([freqs[466], freqs[494]], [0.71042692, 0.71042692], color='red', label='Durchmesser')
#plt.axvline(x=freqs[468], color='red', linestyle='--' )
#plt.axvline(x=freqs[493], color='red', linestyle='--' )
plt.legend(loc='upper right')
plt.xlabel(r'$f$ / Hz')
plt.ylabel(r'Amplitude')
plt.grid(alpha=0.5)
plt.savefig("build/echo_gradient.pdf")
gamma = 2.675 * 10**8
h = 14000
print(2 * np.pi * h /(gamma * 0.0042))