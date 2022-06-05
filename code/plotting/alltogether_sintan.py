import matplotlib.pyplot as plt
from math import *
import numpy as np
import csv
import pandas as pd
import seaborn as sns
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

alldf = pd.read_csv("data/alltogether.csv",
                 names = ['phi_22', 'phi_23', 'phi_32', 'x_1', 'x_2', 'x_3',
                 'x_4', 'a_13', 'a_22', 'a_23', 'a_32', 'phiA', 'tan_theta_12',
                 'sin_theta_13', 'precision'])

mask_1sigma = (alldf['precision'] == 2)
mask_3sigma = (alldf['precision'] == 1)
mask_tot = (alldf['precision'] >= 1)
mask_rest = (alldf['precision'] == 0)
df = alldf[(alldf['sin_theta_13'] <= 1)]

theta_13_exp = np.array([8.62, 8.5, 8.74, 8.25, 8.98])
theta_12_exp = np.array([33.45, 32.7, 34.22, 31.27, 35.87])
theta_13_exp = theta_13_exp*np.pi/180
theta_12_exp = theta_12_exp*np.pi/180

fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', figsize=(10, 10),
                        gridspec_kw={'height_ratios': [1, 3],
                        'width_ratios': [3, 1]})
fig.subplots_adjust(wspace=0, hspace=0)

axs[0, 1].remove()

axs[1, 0].set_xlabel(r"$\sin\theta_{13}$")
#axs[1, 0].set_xscale('log')
axs[1, 0].set_ylabel(r"$\tan\theta_{12}$")
#axs[1, 0].set_xlim(1E-5,2)
axs[1, 0].set_xlim(-0.1,1.1)
axs[1, 0].set_ylim(-0.1,2.2)
axs[1, 0].grid()
the_bar = inset_axes(axs[1, 0], width="3%", height="30%", loc=1,
borderpad=1)
sns.kdeplot(x=np.power((alldf['sin_theta_13'].iloc[0:100000]),1),
            y=alldf['tan_theta_12'],
            color='blue', fill=True, cut=0, clip=[(0,1),(0,3)],
            cmap="Blues", ax=axs[1, 0], cbar=True,
            cbar_ax=the_bar, cbar_kws={'orientation':'vertical', 'format':'%.2f',
            'label':r'KDE'})
the_bar.yaxis.set_ticks_position('left')
the_bar.yaxis.set_label_position('left')
axs[1, 0].axhline(y=np.tan(theta_12_exp[0]), c="darkred", label=r"bfp")
axs[1, 0].axhspan(ymin=np.tan(theta_12_exp[3]), ymax=np.tan(theta_12_exp[4]), color='darkred', alpha=0.3, label=r"$3\sigma$ range")
axs[1, 0].axvline(x=np.power(np.sin(theta_13_exp[0]),1), c="darkred")
axs[1, 0].axvspan(xmin=np.power(np.sin(theta_13_exp[3]),1), xmax=np.power(np.sin(theta_13_exp[4]),1), color='darkred', alpha=0.3)

axs[0, 0].set_ylabel(r"Density of occurences")
axs[0, 0].set_yticks(np.arange(1, 4.1, 1))
axs[0, 0].grid()
"""
counts_tot, bins_tot = np.histogram(np.power(df['sin_theta_13'][mask_tot],1), bins=np.logspace(-5,0,50))
counts_rest, bins_rest = np.histogram(np.power(df['sin_theta_13'][mask_rest],1), bins=np.logspace(-5,0,50))
axs[0, 0].hist([bins_tot[:-1], bins_rest[:-1]],
                bins=np.logspace(-5, 0, 50), density=False, stacked=True,
                weights=[counts_tot/(np.sum(counts_tot)+np.sum(counts_rest))*9.8,
                counts_rest/(np.sum(counts_tot)+np.sum(counts_rest))*9.8],
                label=['','Possible values'],
                color=['darkred','#4080bf'], alpha=0.3)
"""
axs[0, 0].hist([df['sin_theta_13'][mask_tot],
                df['sin_theta_13'][mask_rest]],
                bins=50, range=(0,1), density=True, stacked=True,
                label=['','Possible values'],
                color=['darkred','#4080bf'], alpha=0.3)
sns.kdeplot(x=np.power(df['sin_theta_13'][df['sin_theta_13'] >= 0],1), color='#4080bf', cut=0,
            ax=axs[0, 0], label=r'KDE')

axs[1, 1].set_xlabel(r"Density of occurences")
axs[1, 1].set_xticks(np.arange(0.2, 2, 0.4))
axs[1, 1].grid()
axs[1, 1].hist([df['tan_theta_12'][mask_tot],
                df['tan_theta_12'][mask_rest]],
                bins=75, range=(0,3), density=True, stacked=True,
                color=['darkred','#4080bf'], alpha=0.3, orientation='horizontal')
sns.kdeplot(y=df['tan_theta_12'][df['tan_theta_12'] <= 3], color='#4080bf',
            ax=axs[1, 1])

fig.legend(loc='lower left', bbox_to_anchor=(0.717,0.7))

plt.savefig('alltogether_sintan.svg', format='svg', bbox_inches='tight')
