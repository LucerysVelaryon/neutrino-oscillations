from math import *
import random
import numpy as np
import time
import barre
import csv

llambda = 0.22 # angle of Cabibbo
Delta_m_atm_square = 0.0025 # eV
m0 = 0.035 # eV
x = np.sqrt(Delta_m_atm_square/np.power(m0,2)-1)

# NuFIT 5.1 (2021) IO without SKM
# bfp, 1sigma range, 3sigma range
theta_13_exp = np.array([8.60, 8.48, 8.72, 8.24, 8.98])*np.pi/180
theta_12_exp = np.array([33.45, 32.71, 34.22, 31.27, 35.87])*np.pi/180

n = 20000000
l = []

barre = barre.BarreDeProgression(titre='Computing...')

file = open("data/alltogether.csv", 'a')
csv_writer = csv.writer(file, delimiter=' ')

for i in range(n):
    temps = time.time()
    random.seed()

    precision = 0
    phi_22 = random.uniform(0,2*np.pi)
    phi_23 = random.uniform(0,2*np.pi)
    phi_32 = random.uniform(0,2*np.pi)
    x_1 = random.uniform(-llambda,llambda)
    x_2 = random.uniform(-llambda,llambda)
    x_3 = random.uniform(-llambda,llambda)
    x_4 = random.uniform(-llambda,llambda)
    a_13 = random.uniform(llambda,10*llambda)
    a_22 = random.uniform(llambda,10*llambda)
    a_23 = random.uniform(llambda,10*llambda)
    a_32 = random.uniform(llambda,10*llambda)
    phiA = phi_32 + phi_23
    K = (1+np.power(x,2))/(-a_23*a_32*np.exp(1j*(2*phi_22+phiA))+a_22*np.exp(1j*(phi_22+2*phiA)))
    tan_theta_12 = np.abs(1-(-4*a_13*K*(a_32*np.exp(1j*phi_23)-a_22*x*np.exp(1j*phi_22))-(1+np.power(x,2))*x_1+x_2+2*x*x_3+np.power(x,2)*x_4)/(2*np.power(1+np.power(x,2),3/2))*llambda)
    sin_theta_13 = np.abs((a_13*K*(a_22*np.exp(1j*phi_22)+a_32*np.exp(1j*phi_32))+x_3-x*(x_2+x*x_3-x_4))/np.power(1+np.power(x,2),3/2)*llambda)

    if np.sin(theta_13_exp[1]) < sin_theta_13 < np.sin(theta_13_exp[2]) and np.tan(theta_12_exp[1]) < tan_theta_12 < np.tan(theta_12_exp[2]) :
        precision = 2
    elif np.sin(theta_13_exp[3]) < sin_theta_13 < np.sin(theta_13_exp[4]) and np.tan(theta_12_exp[3]) < tan_theta_12 < np.tan(theta_12_exp[4]) :
        precision = 1

    l = [phi_22, phi_23, phi_32, x_1, x_2, x_3, x_4, a_13, a_22, a_23, a_32, phiA, tan_theta_12, sin_theta_13, precision]
    csv_writer.writerow(l)

    temps = time.time()-temps
    barre.maj(i*100/n, temps)

file.close()
barre.stop()
