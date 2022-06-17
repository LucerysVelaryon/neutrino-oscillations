import numpy as np

"""
Particle Data Group (2022)
--------------------------
bfp
"""
MTAU = 1776.86e6  # ev

"""
Super-Kamiokande IV (2019)
--------------------------
bfp
"""
DELTA_M_ATM_SQUARE = 0.00253  # eV

"""
NuFIT 5.1 (2021) IO without SK-atm
----------------------------------
angle = array(bfp, 1sigma-, 1sigma+, 3sigma-, 3sigma+)
"""
THETA_12_EXP = np.array([33.45, 32.71, 34.22, 31.27, 35.87])*np.pi/180  # rad
THETA_13_EXP = np.array([8.60, 8.48, 8.72, 8.24, 8.98])*np.pi/180  # rad
THETA_23_EXP = np.array([49.5, 48.3, 50.5, 39.8, 52.1])*np.pi/180  # rad
ANGLES = {'tan_theta_12': [np.tan(theta) for theta in THETA_12_EXP],
          'sin_theta_13': [np.sin(theta) for theta in THETA_13_EXP],
          'tan_theta_23': [np.tan(theta) for theta in THETA_23_EXP]}

"""
Simone Marcinano's model
------------------------
"""
LAMBDA = 0.22  # angle of Cabibbo
M0 = 0.035  # eV
X = np.sqrt(DELTA_M_ATM_SQUARE/np.power(M0, 2) - 1)
DEFAULT_RANGES = {}
for e in [
        {'x_{}'.format(i): (-2*LAMBDA, 2*LAMBDA) for i in range(1, 5)},
        {'a_{}{}'.format(i, j): (LAMBDA, 10*LAMBDA) for i in range(1, 4)
        for j in range(1, 4)},
        {'phi_{}'.format(i): (0, 2*np.pi) for i in [22, 23, 32]},
        {'tan_theta_12': (0, 2), 'sin_theta_13': (0, 1), 'tan_theta_23': (0, 2)}
        ]:
    DEFAULT_RANGES.update(e)
PARAM_NAMES = {}
for e in [
        {'x_{}'.format(i): r'$x_{}$'.format(i) for i in range(1, 5)},
        {'a_{}{}'.format(i, j): r'$a_{{{}{}}}$'.format(i, j)
        for i in range(1, 4) for j in range(1, 4)},
        {'phi_{}'.format(i): r'$\phi_{{{}}}$'.format(i)
        for i in [22, 23, 32]},
        {'tan_theta_12': r'$\tan\theta_{12}$',
        'sin_theta_13': r'$\sin\theta_{13}$',
        'tan_theta_23': r'$\tan\theta_{23}$'}
        ]:
        PARAM_NAMES.update(e)
ALL_METHODS = ['neutrinosector', 'alltogether', 'numerical']
