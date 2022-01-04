import numpy as np 
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from scipy.signal import argrelmax

t, A = np.genfromtxt('data/T1.csv', delimiter=';', skip_header=1, unpack=True) 
A = A * 10**-3
print(A)
def f(t, A0, A1, T1):
    return A0 * np.exp(- t/ T1) + A1

params, pcov = curve_fit(f, t, A, p0=(-1, 1, 0.5))

A0 = params[0]
A1 = params[1]
T1 = params[2]
A0_err = np.absolute(pcov[0][0])**0.5
A1_err = np.absolute(pcov[1][1])**0.5
T1_err = np.absolute(pcov[2][2])**0.5

print('A0: ', A0, A0_err)
print('A1: ', A1, A1_err)
print('T1: ', T1, T1_err)

x = np.linspace(np.min(t), np.max(t), 100)
plt.plot(x, f(x, *params), 'r-', label=r'Ausgleichskurve')
plt.plot(t, A, 'bo', label=r'Messwerte')
plt.semilogx()
plt.xlabel(r'$t / \,$ s')
plt.ylabel(r'$U /$\, V')
plt.legend()
plt.tight_layout()
plt.grid(alpha=0.5)
plt.savefig('build/plot_T1.pdf')

plt.clf()

t, M, g, h = np.genfromtxt('data/mg.csv', delimiter=',', skip_header=1, unpack=True)
mask = (t >= 0)
t = t[mask]
M = M[mask]

M_max = argrelmax(M, order =20)[0]
#M_max = np.append(M_max, 0)

def f(t, M0, M1, T2):
    return M0 * np.exp(- t/T2) + M1 

params, pcov = curve_fit(f, t[M_max], M[M_max])
M0 = params[0]
M1 = params[1]
T2 = params[2]
M0_err = np.absolute(pcov[0][0])**0.5
M1_err = np.absolute(pcov[1][1])**0.5
T2_err = np.absolute(pcov[2][2])**0.5
print('M0: ', M0, M0_err)
print('M1: ', M1, M1_err)
print('T2: ', T2, T2_err)

x = np.linspace(np.min(t[M_max]), np.max(t[M_max]), 1000)
plt.plot(x, f(x, *params), 'r-', label=r'Ausgleichskurve')
plt.plot(t, M, 'b-', linewidth=0.5, label=r'gemessene Spannung')
plt.plot(t[M_max], M[M_max], '.',color='black', label=r'Maxima der Spannung')
plt.xlabel(r'$t / \,$s')
plt.ylabel(r'$U /\, $V')
plt.legend()
plt.tight_layout()
plt.grid(alpha =0.5)
plt.savefig('build/plot_T2.pdf')
plt.clf()

t, M = np.genfromtxt('data/diff.csv', delimiter=';', skip_header=1, unpack=True)
t = t * 10**3
M = M * 10**-3
y = (np.log(M) - 2 * t/T2)
x = t**3
nx = x[-1]
ny = y[-1]
def f(x, a, b):
    return a *x + b 


params, pcov = np.polyfit(x, y, deg=1, cov=True)
a = params[0]
b = params[1]
a_err = np.absolute(pcov[0][0])**0.5
b_err = np.absolute(pcov[1][1])**0.5   
print(a, a_err)
print(b, b_err)
plt.plot(x,y, 'bx', label=r'Messwerte')
plt.plot(x, f(x, *params), 'r-', label=r'Ausgleichsgerade')
plt.grid(alpha=0.5)
plt.tight_layout()
plt.xlabel(r'$\tau^3 / \, \symup{s}^3$')
plt.ylabel(r'$\symup{ln}(M(\tau)) - 2\tau/T_2$')
plt.savefig('build/quanti.pdf', bbox_inches='tight')
plt.clf()
t, M = np.genfromtxt('data/diff.csv', delimiter=';', skip_header=1, unpack=True)
#M = np.delete(M, [6, 7])
#t = np.delete(t, [6, 7])

M = M * 10**-3
def f(t, M0, M1, TD):
    return (M0 * np.exp(-2*t/T2) * np.exp(-t**3 * TD) + M1) 

params, pcov = curve_fit(f, t, M, p0=[1, 0.05, 1000])
M0 = params[0]
M1 = params[1]
TD = params[2]
M0_err = np.absolute(pcov[0][0])**0.5
M1_err = np.absolute(pcov[1][1])**0.5
TD_err = np.absolute(pcov[2][2])**0.5
print('M0: ', M0, M0_err)
print('M1: ', M1, M1_err)
print('TD: ', TD, TD_err)

x = np.linspace(np.min(t), np.max(t), 1000)
plt.plot(t, f(t, *params), 'r-', label=r'Ausgleichskurve')
plt.plot(t, M, 'bx', linewidth=0.5, label=r'Messwerte')
#plt.plot(t[M_max], M[M_max], '.',color='black', label='Maxima der Spannung')
plt.xlabel(r'$t / \, $s')
plt.ylabel(r'$U / \, $V')
plt.legend()
plt.tight_layout()
plt.grid(alpha =0.5)
plt.savefig('build/plot_TD.pdf')
gamma = 2.675 * 10**8
from uncertainties import ufloat
uff = ufloat(TD, TD_err)
D =  3/2*uff/( gamma**2 * 0.079**2)
print("D", D)
r = 294.55 * 1.380649*10**-23 /(6 * np.pi * (1002*10**-6) * D )
print(r)
r_anders = (0.0180152 * 0.74/(4/3 * np.pi * 998.21 * 6.02214076 *10**23))**(1/3)
#r_anders = np.cbrt(3 * 1.80152/998.21/(4 * np.pi))
print(r_anders)