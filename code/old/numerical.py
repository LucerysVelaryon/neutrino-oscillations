from math import *
import random
import numpy as np
import time
import barre
import csv
import sympy as sp

llambda = 0.22 # angle of Cabibbo
Delta_m_atm_square = 0.0025 # eV
m0 = 0.035 # eV
mtau = 1776.86e6 # ev
X = np.sqrt(Delta_m_atm_square/np.power(m0,2)-1)

"""
NuFIT 5.1 (2021) IO without SKM
in the following :
angle = array(bfp, 1sigma-, 1sigma+, 3sigma-, 3sigma+)
"""
theta_13_exp = np.array([8.60, 8.48, 8.72, 8.24, 8.98])*np.pi/180 # rad
theta_12_exp = np.array([33.45, 32.71, 34.22, 31.27, 35.87])*np.pi/180 # rad
theta_23_exp = np.array([49.5, 48.3, 50.5, 39.8, 52.1])*np.pi/180 # rad


n = 100000
l = []

barre = barre.BarreDeProgression()

file = open("data/numerical.csv", 'a')
csv_writer = csv.writer(file, delimiter=' ')

for i in range(n):
    temps = time.time()
    random.seed()

    precision = 0
    x = [0]+[random.uniform(-2*llambda,2*llambda) for _ in range(4)]
    a = np.matrix([[0, 0, 0, 0]]+[[0]+[random.uniform(llambda, 10*llambda) for _ in range(3)] for _ in range(3)])
    a[3,3] = 1
    phi_22 = random.uniform(0, 2*np.pi)
    phi_23 = random.uniform(0, 2*np.pi)
    phi_32 = random.uniform(0, 2*np.pi)

    m_nu = sp.Matrix(3, 3, [x[1]*llambda, 1, X,
                            1, x[2]*llambda, x[3]*llambda,
                            X, x[3]*llambda, x[4]*llambda])
    m_l = a[1:,1:]*sp.Matrix(3, 3, [np.power(llambda, 5), np.power(llambda, 3), llambda,
                                    np.power(llambda, 6), np.power(llambda, 2)*np.exp(1j*phi_22), np.exp(1j*phi_23),
                                    np.power(llambda, 6), np.power(llambda, 2)*np.exp(1j*phi_32), 1])

    U_nu = sp.Matrix([list(tup[2][0]) for tup in (m_nu*m_nu.H).eigenvects()]).transpose()
    U_l = sp.Matrix([list(tup[2][0]) for tup in (m_l*m_l.H).eigenvects()]).transpose()
    U_PMNS = U_l.H*U_nu

    sin_theta_13 = np.abs(U_PMNS[0,2])
    tan_theta_12 = np.abs(U_PMNS[0,1])/np.abs(U_PMNS[0,0])
    tan_theta_23 = np.abs(U_PMNS[1,2])/np.abs(U_PMNS[2,2])

    if (np.sin(theta_13_exp[1]) < sin_theta_13 < np.sin(theta_13_exp[2])
        and np.tan(theta_12_exp[1]) < tan_theta_12 < np.tan(theta_12_exp[2])
        and np.tan(theta_23_exp[1]) < tan_theta_23 < np.tan(theta_23_exp[2])) :
        precision = 2
    elif (np.sin(theta_13_exp[3]) < sin_theta_13 < np.sin(theta_13_exp[4])
          and np.tan(theta_12_exp[3]) < tan_theta_12 < np.tan(theta_12_exp[4])
          and np.tan(theta_23_exp[3]) < tan_theta_23 < np.tan(theta_23_exp[4])) :
        precision = 1

    l = [phi_22, phi_23, phi_32]+x[1:]+[a[i,j] for i in range(1,4) for j in range(1,4)]+[sin_theta_13, tan_theta_12, tan_theta_23, precision]
    csv_writer.writerow(l)

    temps = time.time()-temps
    barre.maj(i*100/n, temps)

file.close()
barre.stop()
