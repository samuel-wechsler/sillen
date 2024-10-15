"""
This script defines a function, `plot_sillen`, which generates a Sillén diagram for acids with n dissociable protons. 
The Sillén diagram is based on their respective pKa values and initial concentration C. It uses the `Acid` class 
from the 'acid' module for calculations.
"""

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
import interplot as ip

from acid import Acid


PHG_LINEWIDTH = 10
PHG_OPACITY = 0.4


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


def plot_sillen_interactive(
    acids,
    deprot_levels=None,
    pH_range=np.linspace(0, 14, 100),
    fig=None,
    **kwargs,
):
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

    def mpl_grid(fig, ax):
        ax[0, 0].grid(True, which='major', linestyle='-',
                linewidth=0.5, color='black')
        ax[0, 0].grid(True, which='minor', linestyle='--',
                linewidth=0.2, color='gray')
        ax[0, 0].xaxis.set_major_locator(plt.MultipleLocator(1))
        ax[0, 0].yaxis.set_major_locator(plt.MultipleLocator(1))
        ax[0, 0].xaxis.set_minor_locator(plt.MultipleLocator(0.2))
        ax[0, 0].yaxis.set_minor_locator(plt.MultipleLocator(0.2))
        return fig, ax

    # create plot
    args = dict(
        fig=fig,
        interactive=False,
        title="Sillén diagram",
        xlabel="pH",
        ylabel="log[Spezies]",
        xlim=(0, 14),
        ylim=(-14, 0),
        mpl_custom_func=mpl_grid
    )
    args.update(kwargs)
    fig = ip.Plot.init(
        **args
    )

    phg_left = [10**-pH for pH in pH_range]
    phg_right = [10**(-14 + pH) for pH in pH_range]

    for acid, deprot_level in ip.zip_smart(acids, deprot_levels):
        # create list of concentrations for each protonation state
        for i in range(len(acid.pKa)+1):
            logc = [acid.logc(i, pH) for pH in pH_range]
            fig.add_line(pH_range, logc, label=re.sub(r'[{}_^$]?(:?\\mathrm)?', r'', acid.__repr__(i)))
            if deprot_level is not None and i < deprot_level:
                phg_left += 10**np.array(logc)
            if deprot_level is not None and i > deprot_level:
                phg_right += 10**np.array(logc)

    if deprot_levels is not None:
        fig.add_line(pH_range, np.log10(phg_left), label='PHG left side', linewidth=PHG_LINEWIDTH, opacity=PHG_OPACITY, color="blue")
        fig.add_line(pH_range, np.log10(phg_right), label='PHG right side', linewidth=PHG_LINEWIDTH, opacity=PHG_OPACITY, color="red")

    # plot H, and OH
    fig.add_line(pH_range, [-pH for pH in pH_range],
            linewidth=0.5, label='H+', color='black')
    fig.add_line(pH_range, [-14 + pH for pH in pH_range], linewidth=0.5,
            label='OH-', color='black')

    fig.post_process()

    return fig

# a = Acid([4.75, 8.4, 12.3, 13.5], 0.1, name="A")
b = Acid([3.13, 4.76, 6.4], 0.1, name="PO4", charge=-3)

if __name__ == "__main__":
    plot_sillen([b])
