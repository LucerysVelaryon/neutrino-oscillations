import matplotlib.pyplot as plt
from math import *
import random
import numpy as np
from scipy.interpolate import make_interp_spline
import time

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

Tan = []
Sin = []
tan_max = 2
sin_max = 1
pas = 1000
Tan_bis = np.zeros(pas)
Sin_bis = np.zeros(pas)
n = 1000000
focus = True

for i in range(n):
    random.seed()

    phi = np.zeros((4,4))
    if focus :
        phi[2,2] = random.uniform(0,np.pi/2)
        phi[2,3] = random.uniform(3*np.pi/2-0.2,3*np.pi/2+0.2)
        phi[3,2] = random.uniform(3*np.pi/2,2*np.pi)
    else :
        phi[2,2] = random.uniform(0,2*np.pi)
        phi[2,3] = random.uniform(0,2*np.pi)
        phi[3,2] = random.uniform(0,2*np.pi)

    x_i = np.array([random.uniform(-llambda,llambda) for x in range(5)])
    a = np.mat([[random.uniform(llambda,3) for x in range(4)] for x in range(4)])
    phiA = phi[3,2] + phi[2,3]
    K = (1+np.power(x,2))/(-a[2,3]*a[3,2]*np.exp(1j*(2*phi[2,2]+phiA))+a[2,2]*np.exp(1j*(phi[2,2]+2*phiA)))
    tan_theta_12 = np.abs(1-(-4*a[1,3]*K*(a[3,2]*np.exp(1j*phi[2,3])-a[2,2]*x*np.exp(1j*phi[2,2]))-(1+np.power(x,2))*x_i[1]+x_i[2]+2*x*x_i[3]+np.power(x,2)*x_i[4])/(2*np.power(1+np.power(x,2),3/2))*llambda)
    sin_theta_13 = np.abs((a[1,3]*K*(a[2,2]*np.exp(1j*phi[2,2])+a[3,2]*np.exp(1j*phi[3,2]))+x_i[3]-x*(x_i[2]+x*x_i[3]-x_i[4]))/np.power(1+np.power(x,2),3/2)*llambda)

    if sin_theta_13**2 <= sin_max and tan_theta_12 <= tan_max :
        Tan.append(tan_theta_12)
        Tan_bis[int(tan_theta_12/(tan_max/pas))] += 1
        Sin.append(sin_theta_13)
        Sin_bis[int(sin_theta_13/(sin_max/pas))] += 1

Sin_bis = np.divide(Sin_bis, np.sum(Sin_bis)*sin_max/pas)
Tan_bis = np.divide(Tan_bis, np.sum(Tan_bis)*tan_max/pas)

tan_linspace = np.linspace(0, tan_max, pas)
sin_linspace = np.linspace(0, sin_max, pas)

def smooth(step_curve, step, occurences):
    curve_copy = np.copy(step_curve)
    j_max = int(pas/100)
    width_max = (np.max(step_curve)-np.min(step_curve))/2
    for i in range(pas) :
        count = 0
        sum = 0
        var_j = j_max
        while not((0 <= i+var_j < pas) and (0 <= i-var_j < pas)) :
            var_j -= 1
        while not((np.abs(step_curve[i+var_j]-step_curve[i]) < width_max) and (np.abs(step_curve[i-var_j]-step_curve[i]) < width_max)) :
            var_j -= 1
        for j in range(-var_j,var_j) :
            count += 1
            sum += step_curve[i+j]
        if var_j != 0 :
            curve_copy[i] = sum/count
    return curve_copy

def recur_smooth(step_curve, step, occurences, nb_recur) :
    if nb_recur == 0 :
        return smooth(step_curve, step, occurences)
    else :
        return recur_smooth(smooth(step_curve, step, occurences), step, occurences, nb_recur-1)

Tan_ter = recur_smooth(Tan_bis, pas, n, 4)
Sin_ter = recur_smooth(Sin_bis, pas, n, 4)

im3 = plt.subplot(221)
im3.scatter(np.power(Sin,2), Tan, s=1, c='green', alpha=0.3)
im3.set_xscale('log')
im3.set_xlabel(r"$\sin^2(\theta_{13})$")
im3.set_ylabel(r"$\tan(\theta_{12})$")
im3.grid()
im3.axhline(y=np.tan(theta_12_exp[0]), c="darkred", label=r"bfp NuFIT 5.1 w/o SKM")
im3.axhspan(ymin=np.tan(theta_12_exp[3]), ymax=np.tan(theta_12_exp[4]), color='darkred', alpha=0.3, label=r"$3\sigma$ NuFIT 5.1 w/o SKM")
im3.axvline(x=np.power(np.sin(theta_13_exp[0]),2), c="darkred")
im3.axvspan(xmin=np.power(np.sin(theta_13_exp[3]),2), xmax=np.power(np.sin(theta_13_exp[4]),2), color='darkred', alpha=0.3)
im3.legend()

im3.text(100, 1.8, r"Low energy model plots, based on pp.65-75 of", fontsize='x-large')
im3.text(100, 1.65, r"S. Marciano's Master Thesis", fontsize='x-large')
im3.text(100, 1.2, r"Ranges and values of the parameters", fontweight='bold')
im3.text(100, 1.05, r"$n=$"+str(int(n/1000))+" thousands random occurences")
im3.text(100, 0.9, r"$\lambda =$"+str(llambda))
im3.text(100, 0.75, r"$x=$"+str(x))
if focus :
    im3.text(100, 0.6, r"$\varphi_{22} \in [0,\pi/2] \ , \ \varphi_{23} \in [3\pi/2-0.2,3\pi/2+0.2] \ , \varphi_{32} \in [3\pi/2,2\pi]$")
else :
    im3.text(100, 0.6, r"$\varphi_{22} \in [0,2\pi] \ , \ \varphi_{23} \in [0,2\pi] \ , \varphi_{32} \in [0,2\pi]$")
im3.text(100, 0.45, r"$x_i\in [-\lambda,\lambda], \ \forall\ i \in \{1,2,3,4\}$")
im3.text(100, 0.3, r"$a_{ij}\in [\lambda,3], \ \forall\ i,j \in \{1,2,3\}$")

im1 = plt.subplot(223)
im1.hist(Tan, bins=50, range=(0,2), density=True, label='Density histogram', color='green', alpha=0.3)
#im1.plot(tan_linspace, Tan_bis)
im1.plot(tan_linspace, Tan_ter, label='Smoothed distribution', c='green')
im1.set_xlabel(r"$\tan(\theta_{12})$")
im1.set_ylabel(r"Density of occurences")
im1.grid()
im1.axvline(x=np.tan(theta_12_exp[0]), c="darkred", label=r"bfp NuFIT 5.1 w/o SKM")
im1.axvspan(xmin=np.tan(theta_12_exp[3]), xmax=np.tan(theta_12_exp[4]), color='darkred', alpha=0.3, label=r"$3\sigma$ NuFIT 5.1 w/o SKM")
im1.legend()

im2 = plt.subplot(224)
im2.hist(Sin, bins=50, range=(0,1), density=True, label='Density histogram', color='green', alpha=0.3)
#im2.plot(sin_linspace, Sin_bis)
im2.plot(sin_linspace, Sin_ter, label='Smoothed distribution', c='green')
im2.set_xlabel(r"$\sin(\theta_{13})$")
im2.set_ylabel(r"Density of occurences")
im2.grid()
im2.axvline(x=np.sin(theta_13_exp[0]), c="darkred", label=r"bfp NuFIT 5.1 w/o SKM")
im2.axvspan(xmin=np.sin(theta_13_exp[3]), xmax=np.sin(theta_13_exp[4]), color='darkred', alpha=0.3, label=r"$3\sigma$ NuFIT 5.1 w/o SKM")
im2.legend()

#plt.suptitle(r"Low energy model plots, based on pp.65-75 of S. Marciano's Master Thesis")
plt.show()
#plt.savefig("graphes/"+str(time.time())+".pdf", format='pdf')
