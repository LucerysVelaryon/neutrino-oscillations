import matplotlib.pyplot as plt
from math import *
import numpy as np
import csv
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

df = pd.read_csv("data/alltogether.csv",
                 names = ['phi_22', 'phi_23', 'phi_32', 'x_1', 'x_2', 'x_3',
                 'x_4', 'a_13', 'a_22', 'a_23', 'a_32', 'phiA', 'tan_theta_12',
                 'sin_theta_13', 'precision'])

mask_1sigma = (df['precision'] == 2)
mask_3sigma = (df['precision'] == 1)
mask_tot = (df['precision'] >= 1)

def plotting_phi(phi1, phi2):
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10, 10),
                            gridspec_kw={'height_ratios': [1, 3],
                            'width_ratios': [3, 1]})
    fig.subplots_adjust(wspace=0, hspace=0)

    axs[0, 1].remove()

    axs[1, 0].set_xlabel(r"$\phi_{"+phi1+"}$")
    axs[1, 0].set_ylabel(r"$\phi_{"+phi2+"}$")
    axs[1, 0].grid()
    axs[1, 0].xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
    axs[1, 0].xaxis.set_major_locator(tck.MultipleLocator(base=0.5))
    axs[1, 0].yaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
    axs[1, 0].yaxis.set_major_locator(tck.MultipleLocator(base=0.5))
    sns.kdeplot(x=np.divide(df[mask_tot]['phi_'+phi1],np.pi),
                y=np.divide(df[mask_tot]['phi_'+phi2],np.pi),
                color='green', fill=True, cut=0,
                cmap="Greens", ax=axs[1, 0], cbar=True,
                cbar_ax=inset_axes(axs[1, 0], width="3%", height="30%", loc=3,
                borderpad=1), cbar_kws={'orientation':'vertical', 'format':'%.2f',
                'label':r'KDE 3$\sigma$ range'})
    axs[1, 0].scatter(np.divide(df[mask_1sigma]['phi_'+phi1],np.pi),
                      np.divide(df[mask_1sigma]['phi_'+phi2],np.pi),
                      label=r'$1\sigma$ range', s=1, c='red', alpha=0.8)

    axs[0, 0].set_ylabel(r"Density of occurences")
    axs[0, 0].set_yticks(np.arange(0.1, 0.7, 0.2))
    axs[0, 0].set_xlim(-0.1,2.1)
    axs[0, 0].grid()
    axs[0, 0].hist([np.divide(df[mask_1sigma]['phi_'+phi1],np.pi),
                    np.divide(df[mask_3sigma]['phi_'+phi1],np.pi)],
                    bins=50, range=(0,2), density=True, stacked=True,
                    label=[r'$1\sigma$ range','$3\sigma$ range'],
                    color=['red','green'], alpha=0.3)
    sns.kdeplot(x=np.divide(df[mask_tot]['phi_'+phi1],np.pi), color = 'green', cut=0,
                ax=axs[0, 0], label=r'KDE $3\sigma$ range')

    axs[1, 1].set_xlabel(r"Density of occurences")
    axs[1, 1].set_xticks(np.arange(0.1, 0.7, 0.2))
    axs[1, 1].set_ylim(-0.1,2.1)
    axs[1, 1].grid()
    axs[1, 1].hist([np.divide(df[mask_1sigma]['phi_'+phi2],np.pi),
                    np.divide(df[mask_3sigma]['phi_'+phi2],np.pi)],
                    bins=50, range=(0,2), density=True, stacked=True,
                    color=['red','green'], alpha=0.3, orientation='horizontal')
    sns.kdeplot(y=np.divide(df[mask_tot]['phi_'+phi2],np.pi), color='green', cut=0,
                ax=axs[1, 1], clip=(0,2))

    fig.legend(loc='lower left', bbox_to_anchor=(0.717,0.7))

    plt.savefig('alltogether_phi'+phi1+'_phi'+phi2+'.svg', format='svg', bbox_inches='tight')

plotting_phi('22', '23')
plotting_phi('22', '32')
plotting_phi('23', '32')