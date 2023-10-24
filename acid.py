"""
Module for solving the cubical, quadratic, and linear equations for proton concetration.
"""

import matplotlib.pyplot as plt
import numpy as np
from math import pow

# define constants
K_w = 10**-14


class Acid:
    def __init__(self, pKa, C):
        # TODO: check if pKa is a list, etc.
        self.pKa = pKa
        self.Ka = [10**(-pKa[i]) for i in range(len(pKa))]
        self.C = C

    def alpha(self, i, pH):
        """
        Method returns the fraction of the acid that is dissociated as function of pH.

        Valid range for i: 0, 1, ..., n where n = number of dissociable protons.
        """
        H = 10**(-pH)
        n = len(self.Ka)

        # calculate denominator
        enum = pow(H, n-i)

        for j in range(1, i+1):
            enum *= self.Ka[j-1]

        # calculate value of enumerator
        denom = pow(H, n)

        for k in range(1, n+1):

            p = pow(H, n-k)

            for k in range(1, k+1):
                p *= self.Ka[k-1]

            denom += p

        return enum / denom

    def __repr__(self, i):
        n = len(self.pKa)
        if i == 0:
            return f"$H_{n}A$"
        if i == n:
            return f"$A^{{{-(n)}}}$"

        return f"$H_{{{n - i}}}A^{{{-i}}}$"

    def c(self, i, pH):
        """
        Method returns concentration of species with protonation state i for a given pH.
        Again: i = 0, 1, ..., n.
        """
        return self.C * self.alpha(i, pH)

    def logc(self, i, pH):
        """
        Method returns the log10 concentration of species with protonation state i for a given pH.
        Again: i = 0, 1, ..., n.
        """
        return np.log10(self.c(i, pH))
