
"""
This script defines an Acid class to model the behavior of acids in a solution. It includes methods to calculate
the fraction of acid dissociation, concentration of species with specific protonation states, and more.

Author: Samuel Wechsler
"""
import matplotlib.pyplot as plt
import numpy as np
from math import pow
import sympy as sp

# define constants
K_w = 10**-14


class Acid:
    def __init__(self, pKa, C, name="A", charge=0):
        """
        Initialize an Acid object.

        Args:
            pKa (list of float): List of pKa values for the acid.
            C (float): Initial acid concentration.
            name (str): name of acid.
        """
        # TODO: check if pKa is a list, etc.
        self.pKa = pKa
        self.Ka = [10**(-pKa[i]) for i in range(len(pKa))]
        self.C = C

        # check if string input is valid
        if "$" in name:
            raise Exception(
                "Invalid character in acid name. Please avoid using '$' as it may not render properly.")
        self.name = name
        self.charge = charge

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
        protons = n - i  # number of protons
        charge = self.charge + protons  # charge of species

        # define charge symbols
        charge_symbols = {
            0: "",
            1: "+",
            -1: "-",
        }

        # format charge
        charge_str = charge_symbols.get(
            charge, str(abs(charge)) + ("+" if charge > 0 else "-")
        )

        if protons == 0:
            return f"${self.name}^{{{charge_str}}}$"
        elif protons == 1:
            return f"$H{sp.latex(sp.sympify(self.name))}^{{{charge_str}}}$"
        else:
            return f"$H_{protons}{sp.latex(sp.sympify(self.name))}^{{{charge_str}}}$"

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
