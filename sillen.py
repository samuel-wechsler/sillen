import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

from acid import Acid


def plot_sillen(pKa, C, pH_range=np.linspace(0, 14, 100)):
    """
    Plots sillén diagram for acid n dissociable protons given their respective pKa values
    and initial concentration C.
    """

    # create acid object
    acid = Acid(pKa, C)

    # create plot
    fig, ax = plt.subplots(figsize=(8, 6))

    ax.set_xlabel('pH', fontsize=16)
    ax.set_ylabel('log[Spezies]', fontsize=16)

    # Set the x and y axis limits
    ax.set_xlim(0, 14)
    ax.set_ylim(-14, 0)
    ax.set_aspect('equal')

    # Set the grid properties
    ax.grid(True, which='major', linestyle='-', linewidth=0.5, color='black')
    ax.grid(True, which='minor', linestyle='--', linewidth=0.2, color='gray')
    ax.minorticks_on()

    # Major grid for integer values on the y-axis
    ax.yaxis.set_major_locator(plt.MultipleLocator(1))

    # Set a light green background
    ax.set_facecolor('#f0f0f0')

    # Choose a color palette from Seaborn
    colors = sns.color_palette("tab10", len(acid.Ka)+1)

    # create list of concentrations for each protonation state
    for i in range(len(pKa)+1):
        logc = [acid.logc(i, pH) for pH in pH_range]
        ax.plot(pH_range, logc, label=acid.__repr__(i), color=colors[i])

    # plot H, and OH
    ax.plot(pH_range, [-pH for pH in pH_range],
            linewidth=0.5, label='$H^+$', color='black')

    ax.plot(pH_range, [-14 + pH for pH in pH_range], linewidth=0.5,
            label='$OH^-$', color='black')

    #ax.legend(loc='upper right', bbox_to_anchor=(1.1, 1), borderaxespad=0.0)
    ax.legend(bbox_to_anchor=(1, 1), loc=2, frameon=False, fontsize=16)

    plt.show()


plot_sillen([3.13, 4.76, 6.4], 1)
