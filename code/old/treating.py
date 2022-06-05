import matplotlib.pyplot as plt
import matplotlib.ticker as tck
from math import *
import random
import numpy as np
import time
import barre
import csv
import pandas as pd
import seaborn as sns

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

df = pd.read_csv("test.csv",
                 names = ['phi_22', 'phi_23', 'phi_32', 'x_1', 'x_2', 'x_3',
                 'x_4', 'a_13', 'a_22', 'a_23', 'a_32', 'phiA', 'tan_theta_12',
                 'sin_theta_13', 'precision'])

mask_1sigma = (df['precision'] == 2)
mask_3sigma = (df['precision'] == 1)
mask_tot = (df['precision'] >= 1)

im1 = plt.subplot(231)
im1.hist([df[mask_1sigma]['phi_22'], df[mask_3sigma]['phi_22']],
         bins=1000, range=(0,2*np.pi), density=True, stacked=True,
         label=[r'NuFIT 5.1 $1\sigma$','NuFIT 5.1 $3\sigma$'],
         color=['red','green'], alpha=0.3)
sns.kdeplot(df[mask_tot]['phi_22'], color = 'green', cut=0)
im1.set_xlabel(r"$\varphi_{22}$")
im1.set_ylabel(r"Density of occurences")
im1.grid()

im2 = plt.subplot(232)
im2.hist([df[mask_1sigma]['phi_23'], df[mask_3sigma]['phi_23']], bins=50, range=(0,2*np.pi), density=True, stacked=True, label=[r'NuFIT 5.1 $1\sigma$','NuFIT 5.1 $3\sigma$'], color=['red','green'], alpha=0.3)
df[mask_tot]['phi_23'].plot.kde(color = 'green', xlim=(-0.2,2*np.pi+0.2))
im2.set_xlabel(r"$\varphi_{23}$")
im2.set_ylabel(r"Density of occurences")
im2.grid()
im2.legend()

im3 = plt.subplot(233)
im3.hist([df[mask_1sigma]['phi_32'], df[mask_3sigma]['phi_32']], bins=50, range=(0,2*np.pi), density=True, stacked=True, label=[r'NuFIT 5.1 $1\sigma$','NuFIT 5.1 $3\sigma$'], color=['red','green'], alpha=0.3)
df[mask_tot]['phi_32'].plot.kde(color = 'green', xlim=(-0.2,2*np.pi+0.2))
im3.set_xlabel(r"$\varphi_{32}$")
im3.set_ylabel(r"Density of occurences")
im3.grid()
im3.legend()

im4 = plt.subplot(234)
sns.kdeplot(x=np.divide(df[mask_1sigma]['phi_22'],np.pi), y=np.divide(df[mask_1sigma]['phi_23'],np.pi),
            color='g', fill=True, cut=0,
            cmap="Greens")
im4.scatter(np.divide(df[mask_1sigma]['phi_22'],np.pi), np.divide(df[mask_1sigma]['phi_23'],np.pi), s=1, c='red', alpha=0.3)
im4.scatter(np.divide(df[mask_3sigma]['phi_22'],np.pi), np.divide(df[mask_3sigma]['phi_23'],np.pi), s=1, c='green', alpha=0.1)
im4.set_xlabel(r"$\phi_{22}$")
im4.set_ylabel(r"$\phi_{23}$")
im4.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im4.xaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im4.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im4.yaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im4.grid()

im5 = plt.subplot(235)
im5.scatter(np.divide(df[mask_1sigma]['phi_22'],np.pi), np.divide(df[mask_1sigma]['phi_32'],np.pi), s=1, c='green', alpha=0.3)
im5.scatter(np.divide(df[mask_3sigma]['phi_22'],np.pi), np.divide(df[mask_3sigma]['phi_32'],np.pi), s=1, c='red', alpha=0.1)
im5.set_xlabel(r"$\phi_{22}$")
im5.set_ylabel(r"$\phi_{32}$")
im5.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im5.xaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im5.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im5.yaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im5.grid()

im6 = plt.subplot(236)
im6.scatter(np.divide(df[mask_1sigma]['phi_23'],np.pi), np.divide(df[mask_1sigma]['phi_32'],np.pi), s=1, c='green', alpha=0.3)
im6.scatter(np.divide(df[mask_3sigma]['phi_23'],np.pi), np.divide(df[mask_3sigma]['phi_32'],np.pi), s=1, c='red', alpha=0.1)
im6.set_xlabel(r"$\phi_{23}$")
im6.set_ylabel(r"$\phi_{32}$")
im6.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im6.xaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im6.yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
im6.yaxis.set_major_locator(tck.MultipleLocator(base=0.5))
im6.grid()

plt.show()
