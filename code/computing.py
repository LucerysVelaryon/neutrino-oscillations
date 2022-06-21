import random
import time
import csv
import numpy as np
import sympy as sp
import barre
from exp_data import *


# This algorithm is the core of the project, computing the neutrino mixing
# angles mostly using Simone Marciano's Master Thesis, but also articles on
# similar models by David Meloni and Giorgio Arcadi (see the Overleaf for the
# refs). All the numerical data is then saved in csv files, easy to open using
# Pandas.


def compute(methods=ALL_METHODS, n=100000):
    # A loading bar just to see the approximate remaining time.
    bar = barre.BarreDeProgression(titre='Computing...')
    files = []
    csv_writers = []
    for i in range(len(methods)):
        file = open('data/' + methods[i] + '.csv', 'a')
        csv_writers.append(csv.writer(file, delimiter=' '))
        files.append(file)

    for i in range(n):
        t = time.time()
        random.seed()
        count = 0

        # Creating additionnal 0-value in array and matrix so that the indices
        # in the code comply with those in the model.
        x = [0] + [random.uniform(DEFAULT_RANGES['x_{}'.format(i)][0],
                   DEFAULT_RANGES['x_{}'.format(i)][1]) for i in range(1, 5)]
        a = np.matrix([[0]*4] + [
            [0]+[random.uniform(DEFAULT_RANGES['a_{}{}'.format(i, j)][0],
            DEFAULT_RANGES['a_{}{}'.format(i, j)][1]) for i in range(1, 4)]
            for j in range(1, 4)
            ])
        # No need to compute all 9 possible phases for the computings, so
        # creating three variables, instead of the entire matrix, in order to
        # save computationnal time.
        phi_22 = random.uniform(DEFAULT_RANGES['phi_22'][0],
                                DEFAULT_RANGES['phi_22'][1])
        phi_23 = random.uniform(DEFAULT_RANGES['phi_23'][0],
                                DEFAULT_RANGES['phi_23'][1])
        phi_32 = random.uniform(DEFAULT_RANGES['phi_32'][0],
                                DEFAULT_RANGES['phi_32'][1])

        # The two following methods use the analytical formulae computed by
        # Simone Marciano in his thesis, without and with taking into account
        # the contribution from the charged lepton sector.

        if 'neutrinosector' in methods :
            tan_theta_12 = (1 - (-(1 + np.power(X, 2))*x[1] + x[2] + X*(2*x[3]
                            + X*x[4]))/(2*np.power(1 + np.power(X, 2), 3/2)))
            sin_theta_13 = LAMBDA/np.power(1 + np.power(X, 2), 3/2)*(X*x[2]
                           + np.power(X, 2)*x[3] - x[3] - X*x[4])
            csv_writers[count].writerow(x[1:] + [tan_theta_12, sin_theta_13])
            count += 1

        if 'alltogether' in methods :
            phiA = phi_32 + phi_23
            K = (1 + np.power(X, 2))/(-a[2,3]*a[3,2]*np.exp(1j*(2*phi_22
                 + phiA)) + a[2,2]*np.exp(1j*(phi_22 + 2*phiA)))
            tan_theta_12 = np.abs(1 - (-4*a[1,3]*K*(a[3,2]*np.exp(1j*phi_23)
                           - a[2,2]*X*np.exp(1j*phi_22))
                           - (1 + np.power(X, 2))*x[1] + x[2] + 2*X*x[3]
                           + np.power(X, 2)*x[4])/(2*np.power(1
                           + np.power(X, 2), 3/2))*LAMBDA)
            sin_theta_13 = np.abs((a[1,3]*K*(a[2,2]*np.exp(1j*phi_22)
                           + a[3,2]*np.exp(1j*phi_32)) + x[3] - X*(x[2]
                           + X*x[3] - x[4]))/np.power(1
                           + np.power(X, 2), 3/2)*LAMBDA)
            csv_writers[count].writerow(
                x[1:] + [a[1,3], a[2,2], a[2,3], a[3,2], phi_22, phi_23, phi_32,
                tan_theta_12, sin_theta_13]
                )
            count += 1

        # The following method takes into account the reasoning made by Simone
        # Marciano in his thesis, especially concerning the ranges of the
        # parameters of the model, but computes the mixing angles using in a
        # total numerical way, i.e. diagonalizing matrices at the end of the
        # reasoning.

        if 'numerical' in methods :
            m_nu = sp.Matrix(3, 3, [
                x[1]*LAMBDA, 1, X, 1, x[2]*LAMBDA, x[3]*LAMBDA, X, x[3]*LAMBDA,
                x[4]*LAMBDA
                ])
            m_l = sp.Matrix(3, 3, [
                a[1,1]*np.power(LAMBDA, 5),
                a[1,2]* np.power(LAMBDA, 3),
                a[1,3]*LAMBDA,
                a[2,1]*np.power(LAMBDA, 6),
                a[2,2]*np.power(LAMBDA, 2) * np.exp(1j*phi_22),
                a[2,3]*np.exp(1j*phi_23),
                a[3,1]*np.power(LAMBDA, 6),
                a[3,2]*np.power(LAMBDA, 2) * np.exp(1j*phi_32),
                1
                ])
            # The following method used to extract the eigenvectors gives the
            # transposed wanted matrix, so one must tranpose it back. This is
            # due to the package SymPy which does not directly give the
            # eigenvectors but a list of tuples of the form :
            #     [(eigenvalue, multiplicity, span), ...]
            # where the span is hence a list of eigenvectors associated with
            # the eigenvalue. For my analysis, the multiplicity is always 1.
            # For the following matrices, an order in the eigenvalues, hence in
            # the eigenvectors, must be chosen to define a unique U_nu (resp.
            # U_l) matrix.
            # For U_nu, IH is assumed : m1 > m2 >> m3.
            # For U_l, the known charged lepton masses and the definition of
            # m_l imposes ma = m_e < mb = m_mu < mc = m_tau.
            U_nu = sp.Matrix([
                list(tup[2][0]) for tup in sorted((m_nu * m_nu.H).eigenvects(),
                key=lambda e: complex(e[0]).real, reverse=True)
                ]).transpose()
            U_l = sp.Matrix([
                list(tup[2][0]) for tup in sorted((m_l * m_l.H).eigenvects(),
                key=lambda e: complex(e[0]).real)
                ]).transpose()
            U_PMNS = U_l.H * U_nu
            # Warning, due to Python, the indices of the U_PMNS terms used for
            # the computing of the angles are minus 1 the one used in
            # theoretical formulae in the literature (U_13 in an article will
            # be U_02 in this code).
            tan_theta_12 = np.abs(U_PMNS[0,1])/np.abs(U_PMNS[0,0])
            sin_theta_13 = np.abs(U_PMNS[0,2])
            tan_theta_23 = np.abs(U_PMNS[1,2])/np.abs(U_PMNS[2,2])
            csv_writers[count].writerow(
                x[1:] + [a[i,j] for i in range(1, 4) for j in range(1, 4)]
                + [phi_22, phi_23, phi_32, tan_theta_12, sin_theta_13,
                tan_theta_23]
                )
            count += 1

        bar.maj(i*100/n, time.time()-t)

    for file in files :
        file.close()
    bar.stop()
