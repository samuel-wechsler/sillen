"""
This script defines a function, `plot_sillen`, which generates a Sillén diagram for acids with n dissociable protons. 
The Sillén diagram is based on their respective pKa values and initial concentration C. It uses the `Acid` class 
from the 'acid' module for calculations.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from acid import Acid


def plot_sillen(acids, pH_range=np.linspace(0, 14, 100)):
    """
    Plots Sillén diagram for acids with n dissociable protons given their respective pKa values
    and initial concentration C.

    Args:
        pKa (list of float): List of pKa values for the acid.
        C (float): Initial acid concentration.
        pH_range (numpy.ndarray, optional): pH values for the x-axis (default is a range from 0 to 14).

    The function generates a Sillén diagram showing the pH-dependent concentration of species with different protonation states.
    """
    if not isinstance(acids, list):
        acids = [acids]

    # create plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.set_xlabel('pH', fontsize=16)
    ax.set_ylabel('log[Spezies]', fontsize=16)

    # Set the x and y axis limits
    ax.set_xlim(0, 14)
    ax.set_ylim(-14, 0)
    ax.set_aspect('equal')

    # Set the grid properties
    ax.grid(True, which='major', linestyle='-',
            linewidth=0.5, color='black')
    ax.grid(True, which='minor', linestyle='--',
            linewidth=0.2, color='gray')
    ax.minorticks_on()

    # Major grid for integer values on the y-axis
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))

    # Set a light green background
    ax.set_facecolor('#f0f0f0')

    # Choose a color palette from Seaborn
    colors = sns.color_palette("tab10", sum(
        [len(acid.pKa) for acid in acids]) + 1)

    # # Choose the "viridis" color map from Matplotlib
    # colors = plt.cm.viridis(np.linspace(
    #     0, 1, sum([len(acid.pKa) for acid in acids]) + 1))

    for acid in acids:
        # create list of concentrations for each protonation state
        for i in range(len(acid.pKa)+1):
            logc = [acid.logc(i, pH) for pH in pH_range]
            ax.plot(pH_range, logc, label=acid.__repr__(i), color=colors[i])

    # plot H, and OH
    ax.plot(pH_range, [-pH for pH in pH_range],
            linewidth=0.5, label='$H^+$', color='black')

    ax.plot(pH_range, [-14 + pH for pH in pH_range], linewidth=0.5,
            label='$OH^-$', color='black')

    # ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1), borderaxespad=0.0)
    ax.legend(bbox_to_anchor=(1, 1), loc=2, frameon=False, fontsize=16)

    plt.show()


# a = Acid([4.75, 8.4, 12.3, 13.5], 0.1, name="A")
b = Acid([3.13, 4.76, 6.4], 0.1, name=f"PO4", charge=-3)
plot_sillen([b])
