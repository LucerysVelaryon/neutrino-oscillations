import random
import time
import csv
import numpy as np
import sympy as sp
import barre
from exp_data import *


"""
This algorithm is the core of the project, computing the neutrino mixing angles
mostly using Simone Marciano's Master Thesis, but also articles on similar
models by David Meloni and Giorgio Arcadi (see the Overleaf for the refs).
"""


def computing(methods = "all", n = 1000) :
    """
    A loading bar just to see the approximate remaining time.
    """
    bar = barre.BarreDeProgression(titre = 'Computing...')
    if methods == "all" :
        methods = ["neutrinosector", "alltogether", "numerical"]

    files = []
    csv_writers = []
    for i in range(len(methods)) :
        file = open("data/" + methods[i] + ".csv", 'a')
        csv_writers.append(csv.writer(file, delimiter = ' '))
        files.append(file)

    for i in range(n):
        t = time.time()
        random.seed()
        count = 0

        """
        Creating additionnal 0-value in array and matrix so that the
        indexes in the code comply with those in the model.
        """
        x = [0] + [random.uniform(-2*llambda, 2*llambda) for _ in range(4)]
        a = np.matrix([[0, 0, 0, 0]] + [[0] + [random.uniform(llambda, 10*llambda) for _ in range(3)] for _ in range(3)])
        """
        Setting a_33 to 1 to comply with the expression of the charged letpons
        mass matrix in the model, expressed in units of m_tau.
        """
        a[3,3] = 1
        """
        No need to compute all 9 possible phases for the computations, so
        creating three variables to save computationnal time.
        """
        phi_22 = random.uniform(0, 2*np.pi)
        phi_23 = random.uniform(0, 2*np.pi)
        phi_32 = random.uniform(0, 2*np.pi)

        """
        The two following methods use the analytical formulae computed by Simone
        Marciano in his thesis, without and with taking into account the
        contribution from the charged lepton sector.
        """

        if "neutrinosector" in methods :
            tan_theta_12 = 1 - (-(1 + np.power(X, 2))*x[1] + x[2] + X*(2*x[3] + X*x[4]))/(2*np.power(1 + np.power(X, 2), 3/2))
            sin_theta_13 = llambda/np.power(1 + np.power(X, 2), 3/2)*(X*x[2] + np.power(X, 2)*x[3] - x[3] - X*x[4])
            csv_writers[count].writerow(x[1:] + [tan_theta_12, sin_theta_13])
            count += 1

        if "alltogether" in methods :
            phiA = phi_32 + phi_23
            K = (1 + np.power(X, 2))/(-a[2,3]*a[3,2]*np.exp(1j*(2*phi_22 + phiA)) + a[2,2]*np.exp(1j*(phi_22 + 2*phiA)))
            tan_theta_12 = np.abs(1 - (-4*a[1,3]*K*(a[3,2]*np.exp(1j*phi_23) - a[2,2]*X*np.exp(1j*phi_22)) - (1 + np.power(X, 2))*x[1] + x[2] + 2*X*x[3] + np.power(X, 2)*x[4])/(2*np.power(1 + np.power(X, 2), 3/2))*llambda)
            sin_theta_13 = np.abs((a[1,3]*K*(a[2,2]*np.exp(1j*phi_22) + a[3,2]*np.exp(1j*phi_32)) + x[3] - X*(x[2] + X*x[3] - x[4]))/np.power(1 + np.power(X, 2), 3/2)*llambda)
            csv_writers[count].writerow(x[1:] + [a[1,3], a[2,2], a[2,3], a[3,2], phi_22, phi_23, phi_32, tan_theta_12, sin_theta_13])
            count += 1

        """
        The following method takes into account the resoning made by Simone
        Marciano in his thesis, especially concerning the ranges of the
        parameters of the model, but computes the mixing angles using in a total
        numerical way, i.e. diagonalizing matrices at the end of the reasoning.
        """

        if "numerical" in methods :
            m_nu = sp.Matrix(3, 3, [x[1]*llambda, 1, X,
                                    1, x[2]*llambda, x[3]*llambda,
                                    X, x[3]*llambda, x[4]*llambda])
            m_l = a[1:, 1:]*sp.Matrix(3, 3, [np.power(llambda, 5), np.power(llambda, 3), llambda,
                                             np.power(llambda, 6), np.power(llambda, 2)*np.exp(1j*phi_22), np.exp(1j*phi_23),
                                             np.power(llambda, 6), np.power(llambda, 2)*np.exp(1j*phi_32), 1])
            """
            The following method used to extract the eigenvectors gives the
            transposed wanted matrix, so one must tranpose it back.
            """
            U_nu = sp.Matrix([list(tup[2][0]) for tup in (m_nu*m_nu.H).eigenvects()]).transpose()
            U_l = sp.Matrix([list(tup[2][0]) for tup in (m_l*m_l.H).eigenvects()]).transpose()
            U_PMNS = U_l.H*U_nu
            tan_theta_12 = np.abs(U_PMNS[0,1])/np.abs(U_PMNS[0,0])
            sin_theta_13 = np.abs(U_PMNS[0,2])
            tan_theta_23 = np.abs(U_PMNS[1,2])/np.abs(U_PMNS[2,2])
            csv_writers[count].writerow(x[1:] + [a[i,j] for i in range(1,4) for j in range(1,4)] + [phi_22, phi_23, phi_32, tan_theta_12, sin_theta_13, tan_theta_23])
            count += 1

        bar.maj(i*100/n, time.time()-t)

    for file in files :
        file.close()
    bar.stop()
