import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from math import *
import random
import numpy as np
from scipy.interpolate import make_interp_spline
import time
import smooth_function as sm
import barre

llambda = 0.22 # angle of Cabibbo
Delta_m_atm_square = 0.0025 # eV
m0 = 0.035 # eV
x = np.sqrt(Delta_m_atm_square/np.power(m0,2)-1)

# NuFIT 5.1 (2021) without SKM
# bfp, 1sigma range, 3sigma range
theta_13_exp = np.array([8.62, 8.5, 8.74, 8.25, 8.98])
theta_12_exp = np.array([33.45, 32.7, 34.22, 31.27, 35.87])
theta_13_exp = theta_13_exp*np.pi/180
theta_12_exp = theta_12_exp*np.pi/180
#r_exp = 0.03

pas = 1000
n = 10000000

Phi_22 = []
Phi_23 = []
Phi_32 = []
Phi_22_bis = np.zeros(pas)
Phi_23_bis = np.zeros(pas)
Phi_32_bis = np.zeros(pas)

barre1 = barre.BarreDeProgression(titre='Scattering')

for i in range(n):
    random.seed()

    phi = np.zeros((4,4))
    phi[2,2] = random.uniform(0,2*np.pi)
    phi[2,3] = random.uniform(0,2*np.pi)
    phi[3,2] = random.uniform(0,2*np.pi)

    x_i = np.array([random.uniform(-llambda,llambda) for x in range(5)])
    a = np.mat([[random.uniform(llambda,3) for x in range(4)] for x in range(4)])
    phiA = phi[3,2] + phi[2,3]
    K = (1+np.power(x,2))/(-a[2,3]*a[3,2]*np.exp(1j*(2*phi[2,2]+phiA))+a[2,2]*np.exp(1j*(phi[2,2]+2*phiA)))
    tan_theta_12 = np.abs(1-(-4*a[1,3]*K*(a[3,2]*np.exp(1j*phi[2,3])-a[2,2]*x*np.exp(1j*phi[2,2]))-(1+np.power(x,2))*x_i[1]+x_i[2]+2*x*x_i[3]+np.power(x,2)*x_i[4])/(2*np.power(1+np.power(x,2),3/2))*llambda)
    sin_theta_13 = np.abs((a[1,3]*K*(a[2,2]*np.exp(1j*phi[2,2])+a[3,2]*np.exp(1j*phi[3,2]))+x_i[3]-x*(x_i[2]+x*x_i[3]-x_i[4]))/np.power(1+np.power(x,2),3/2)*llambda)
    #r = 2*llambda/np.power(1+np.power(x,2),3/2)*(x_i[1]*(1+np.power(x,2))+x_i[2]+2*x*x_i[3]+x_4*np.power(x,2))

    if np.sin(theta_13_exp[3]) < sin_theta_13 < np.sin(theta_13_exp[4]) and np.tan(theta_12_exp[3]) < tan_theta_12 < np.tan(theta_12_exp[4]) :
        Phi_22.append(phi[2,2])
        Phi_32.append(phi[3,2])
        Phi_23.append(phi[2,3])
        Phi_22_bis[int(phi[2,2]/(2*np.pi/pas))] += 1
        Phi_23_bis[int(phi[2,3]/(2*np.pi/pas))] += 1
        Phi_32_bis[int(phi[3,2]/(2*np.pi/pas))] += 1

    barre1.maj(i*100/n)

#Phi_22_bis = np.divide(Phi_22_bis, np.sum(Phi_22_bis)*2*np.pi/pas)
#Phi_23_bis = np.divide(Phi_23_bis, np.sum(Phi_23_bis)*2*np.pi/pas)
#Phi_32_bis = np.divide(Phi_32_bis, np.sum(Phi_32_bis)*2*np.pi/pas)

phi_linspace = np.linspace(0, 2*np.pi, pas)

im1 = plt.subplot(231)
im1.hist(Phi_22, bins=50, range=(0,2*np.pi), density=True, label='Density histogram', color='green', alpha=0.3)
#im1.plot(phi_linspace, sm.smooth(Phi_22_bis), label='Smoothed distribution', c='green')
im1.set_xlabel(r"$\varphi_{22}$")
im1.set_ylabel(r"Density of occurences")
im1.grid()
im1.legend()

im2 = plt.subplot(232)
im2.hist(Phi_23, bins=50, range=(0,2*np.pi), density=True, label='Density histogram', color='green', alpha=0.3)
#im2.plot(phi_linspace, sm.final_smooth(Phi_23_bis), label='Smoothed distribution', c='green')
im2.set_xlabel(r"$\varphi_{23}$")
im2.set_ylabel(r"Density of occurences")
im2.grid()
im2.legend()

im3 = plt.subplot(233)
im3.hist(Phi_32, bins=50, range=(0,2*np.pi), density=True, label='Density histogram', color='green', alpha=0.3)
#im3.plot(phi_linspace, sm.final_smooth(Phi_32_bis), label='Smoothed distribution', c='green')
im3.set_xlabel(r"$\varphi_{32}$")
im3.set_ylabel(r"Density of occurences")
im3.grid()
im3.legend()

im4 = plt.subplot(234)
im4.scatter(np.divide(Phi_22,np.pi), np.divide(Phi_23,np.pi), s=1, c='green', alpha=0.3)
im4.set_xlabel(r"$\phi_{22}$")
im4.set_ylabel(r"$\phi_{23}$")
im4.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im4.xaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im4.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im4.yaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im4.grid()

im5 = plt.subplot(235)
im5.scatter(np.divide(Phi_22,np.pi), np.divide(Phi_32,np.pi), s=1, c='green', alpha=0.3)
im5.set_xlabel(r"$\phi_{22}$")
im5.set_ylabel(r"$\phi_{32}$")
im5.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im5.xaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im5.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im5.yaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im5.grid()

im6 = plt.subplot(236)
im6.scatter(np.divide(Phi_23,np.pi), np.divide(Phi_32,np.pi), s=1, c='green', alpha=0.3)
im6.set_xlabel(r"$\phi_{23}$")
im6.set_ylabel(r"$\phi_{32}$")
im6.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im6.xaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im6.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im6.yaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im6.grid()


plt.show()
