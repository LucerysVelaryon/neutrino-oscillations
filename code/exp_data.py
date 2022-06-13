import numpy as np

"""
Particle Data Group (2022)
--------------------------
bfp
"""
mtau = 1776.86e6 # ev

"""
Super-Kamiokande IV (2019)
--------------------------
bfp
"""
Delta_m_atm_square = 0.00253 # eV

"""
NuFIT 5.1 (2021) IO without SK-atm
----------------------------------
angle = array(bfp, 1sigma-, 1sigma+, 3sigma-, 3sigma+)
"""
theta_13_exp = np.array([8.60, 8.48, 8.72, 8.24, 8.98])*np.pi/180 # rad
theta_12_exp = np.array([33.45, 32.71, 34.22, 31.27, 35.87])*np.pi/180 # rad
theta_23_exp = np.array([49.5, 48.3, 50.5, 39.8, 52.1])*np.pi/180 # rad

"""
Simone Marcinano's model
------------------------
"""
llambda = 0.22 # angle of Cabibbo
m0 = 0.035 # eV
X = np.sqrt(Delta_m_atm_square/np.power(m0,2)-1)
