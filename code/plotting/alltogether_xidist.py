import matplotlib.pyplot as plt
from math import *
import numpy as np
import csv
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# NuFIT 5.1 (2021) IO without SKM
# bfp, 1sigma range, 3sigma range
theta_13_exp = np.array([8.60, 8.48, 8.72, 8.24, 8.98])*np.pi/180
theta_12_exp = np.array([33.45, 32.71, 34.22, 31.27, 35.87])*np.pi/180

df = pd.read_csv("data/alltogether.csv",
                 names = ['phi_22', 'phi_23', 'phi_32', 'x_1', 'x_2', 'x_3',
                 'x_4', 'a_13', 'a_22', 'a_23', 'a_32', 'phiA', 'tan_theta_12',
                 'sin_theta_13', 'precision'])

df['precision'] = 0

df.loc[lambda df: (df['sin_theta_13'] >= np.sin(theta_13_exp[3]))
                        & (df['sin_theta_13'] <= np.sin(theta_13_exp[4]))
                        & (df['tan_theta_12'] >= np.tan(theta_12_exp[3]))
                        & (df['tan_theta_12'] <= np.tan(theta_12_exp[4])),
          'precision'] = 1

df.loc[lambda df: (df['sin_theta_13'] >= np.sin(theta_13_exp[1]))
                        & (df['sin_theta_13'] <= np.sin(theta_13_exp[2]))
                        & (df['tan_theta_12'] >= np.tan(theta_12_exp[1]))
                        & (df['tan_theta_12'] <= np.tan(theta_12_exp[2])),
          'precision'] = 2

mask_1sigma = (df['precision'] == 2)
mask_3sigma = (df['precision'] == 1)
mask_tot = (df['precision'] >= 1)

def plotting_xi(xa, xb):
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10, 10),
                            gridspec_kw={'height_ratios': [1, 3],
                            'width_ratios': [3, 1]})
    fig.subplots_adjust(wspace=0, hspace=0)

    axs[0, 1].remove()

    axs[1, 0].set_xlabel(r"$x_{"+xa+"}$")
    axs[1, 0].set_ylabel(r"$x_{"+xb+"}$")
    axs[1, 0].grid()
    sns.kdeplot(x=df[mask_tot]['x_'+xa],
                y=df[mask_tot]['x_'+xb],
                color='green', fill=True, cut=0, clip=[(-0.22,0.22),(-0.22,0.22)],
                cmap="Greens", ax=axs[1, 0], cbar=True,
                cbar_ax=inset_axes(axs[1, 0], width="3%", height="30%", loc=3,
                borderpad=1), cbar_kws={'orientation':'vertical', 'format':'%.2f',
                'label':r'KDE 3$\sigma$ range'})
    axs[1, 0].scatter(df[mask_1sigma]['x_'+xa],
                      df[mask_1sigma]['x_'+xb],
                      label=r'$1\sigma$ range', s=1, c='red', alpha=0.8)

    axs[0, 0].set_ylabel(r"Density of occurences")
    axs[0, 0].set_yticks(np.arange(0.5, 3, 0.5))
    axs[0, 0].set_xlim(-0.25,0.25)
    axs[0, 0].grid()
    axs[0, 0].hist([df[mask_1sigma]['x_'+xa],
                    df[mask_3sigma]['x_'+xa]],
                    bins=50, density=True, stacked=True,
                    label=[r'$1\sigma$ range','$3\sigma$ range'],
                    color=['red','green'], alpha=0.3)
    sns.kdeplot(x=df[mask_tot]['x_'+xa], color = 'green', cut=0, clip=(-0.22,0.22),
                ax=axs[0, 0], label=r'KDE $3\sigma$ range')

    axs[1, 1].set_xlabel(r"Density of occurences")
    axs[1, 1].set_xticks(np.arange(1, 3, 1))
    axs[1, 1].set_ylim(-0.25,0.25)
    axs[1, 1].grid()
    axs[1, 1].hist([df[mask_1sigma]['x_'+xb],
                    df[mask_3sigma]['x_'+xb]],
                    bins=50, density=True, stacked=True,
                    color=['red','green'], alpha=0.3, orientation='horizontal')
    sns.kdeplot(y=df[mask_tot]['x_'+xb], color='green', cut=0, clip=(-0.22,0.22),
                ax=axs[1, 1])

    fig.legend(loc='lower left', bbox_to_anchor=(0.717,0.7))

    plt.savefig('alltogether_x'+xa+'_x'+xb+'.svg', format='svg', bbox_inches='tight')

plotting_xi('1', '2')
plotting_xi('1', '3')
plotting_xi('1', '4')
plotting_xi('2', '3')
plotting_xi('2', '4')
plotting_xi('3', '4')
