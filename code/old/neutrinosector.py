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

file = open("data/neutrinosector.csv", 'a')
csv_writer = csv.writer(file, delimiter=' ')

for i in range(n):
    temps = time.time()
    random.seed()

    precision = 0
    x_1 = random.uniform(-llambda,llambda)
    x_2 = random.uniform(-llambda,llambda)
    x_3 = random.uniform(-llambda,llambda)
    x_4 = random.uniform(-llambda,llambda)
    sin_theta_13 = llambda/np.power(1+np.power(x,2),3/2)*(x*x_2+np.power(x,2)*x_3-x_3-x*x_4)
    tan_theta_12 = 1-(-(1+np.power(x,2))*x_1+x_2+x*(2*x_3+x*x_4))/(2*np.power(1+np.power(x,2),3/2))

    if np.sin(theta_13_exp[1]) < sin_theta_13 < np.sin(theta_13_exp[2]) and np.tan(theta_12_exp[1]) < tan_theta_12 < np.tan(theta_12_exp[2]) :
        precision = 2
    elif np.sin(theta_13_exp[3]) < sin_theta_13 < np.sin(theta_13_exp[4]) and np.tan(theta_12_exp[3]) < tan_theta_12 < np.tan(theta_12_exp[4]) :
        precision = 1

    l = [x_1, x_2, x_3, x_4, tan_theta_12, sin_theta_13, precision]
    csv_writer.writerow(l)

    temps = time.time()-temps
    barre.maj(i*100/n, temps)

file.close()
barre.stop()
