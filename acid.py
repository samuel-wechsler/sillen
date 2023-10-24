
"""
This script defines an Acid class to model the behavior of acids in a solution. It includes methods to calculate
the fraction of acid dissociation, concentration of species with specific protonation states, and more.

For the formula used in the `alpha` method, please refer to the source:
"Insert Source Reference Here" (DOI 10.1007/s40828-016-0029-1)

Author: Samuel Wechsler
"""
import matplotlib.pyplot as plt
import numpy as np
from math import pow

# define constants
K_w = 10**-14


class Acid:
    def __init__(self, pKa, C):
        """
        Initialize an Acid object.

        Args:
            pKa (list of float): List of pKa values for the acid.
            C (float): Initial acid concentration.
        """
        # TODO: check if pKa is a list, etc.
        self.pKa = pKa
        self.Ka = [10**(-pKa[i]) for i in range(len(pKa))]
        self.C = C

    def alpha(self, i, pH):
        """
        Calculate fraction of acid dissociated at a given pH.

        Args:
            i (int): Protonation state (0, 1, ..., n).
            pH (float): pH value.

        Returns:
            float: Fraction of acid dissociated.

        Reference:
        Gambi, A., Toniolo, R. Acid–base logarithmic diagrams with computer algebra systems.
        ChemTexts 2, 9 (2016). https://doi.org/10.1007/s40828-016-0029-1
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
        """
        Return a representation of the acid species with a specific protonation state.

        Args:
            i (int): Protonation state (0, 1, ..., n).

        Returns:
            str: Representation of the acid species.
        """
        n = len(self.pKa)
        if i == 0:
            return f"$H_{n}A$"
        if i == n:
            return f"$A^{{{-(n)}}}$"

        return f"$H_{{{n - i}}}A^{{{-i}}}$"

    def c(self, i, pH):
        """
        Calculate concentration of species with a specific protonation state at a given pH.

        Args:
            i (int): Protonation state (0, 1, ..., n).
            pH (float): pH value.

        Returns:
            float: Concentration of the species.
        """
        return self.C * self.alpha(i, pH)

    def logc(self, i, pH):
        """
        Calculate log10 concentration of species with a specific protonation state at a given pH.

        Args:
            i (int): Protonation state (0, 1, ..., n).
            pH (float): pH value.

        Returns:
            float: Log10 concentration of the species.
        """
        return np.log10(self.c(i, pH))
