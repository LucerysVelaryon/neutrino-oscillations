import csv
import numpy as np
import pandas as pd
import seaborn as sb
import smoothing as sm
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import mpl_toolkits.axes_grid1.inset_locator as inloc
from exp_data import *


def load_distributions(methods=ALL_METHODS):
    dataframes = {
        m: pd.read_csv('data/' + m + '.csv', delimiter=' ') for m in methods
        }
    for m in methods :
        dataframes[m] = dataframes[m][dataframes[m]['sin_theta_13'] <= 1]
        if m != 'numerical' :
            dataframes[m].loc[:, 'tan_theta_23'] = np.tan(THETA_23_EXP[0])
        dataframes[m]['precision'] = 0
        dataframes[m].loc[
            lambda df: (df['tan_theta_12'] >= ANGLES['tan_theta_12'][3])
            & (df['tan_theta_12'] <= ANGLES['tan_theta_12'][4])
            & (df['sin_theta_13'] >= ANGLES['sin_theta_13'][3])
            & (df['sin_theta_13'] <= ANGLES['sin_theta_13'][4])
            & (df['tan_theta_23'] >= ANGLES['tan_theta_23'][3])
            & (df['tan_theta_23'] <= ANGLES['tan_theta_23'][4]),
            'precision'
            ] = 1
        dataframes[m].loc[
            lambda df: (df['tan_theta_12'] >= ANGLES['tan_theta_12'][1])
            & (df['tan_theta_12'] <= ANGLES['tan_theta_12'][2])
            & (df['sin_theta_13'] >= ANGLES['sin_theta_13'][1])
            & (df['sin_theta_13'] <= ANGLES['sin_theta_13'][2])
            & (df['tan_theta_23'] >= ANGLES['tan_theta_23'][1])
            & (df['tan_theta_23'] <= ANGLES['tan_theta_23'][2]),
            'precision'
            ] = 2
    return dataframes


def refine_ranges(dataframe, parameter, upperlim=None, bottomlim=None):
    if upperlim != None:
        dataframe = dataframe[dataframe[parameter] <= upperlim]
    if bottomlim != None:
        dataframe = dataframe[dataframe[parameter] >= bottomlim]
    return dataframe


def extend(some_range, log_scale):
    some_range = list(some_range)
    if log_scale:
        margin = np.abs(some_range[1]-some_range[0])
        if some_range[0] == 0 :
            some_range[0] = 1E-5 + margin
    else:
        margin = np.abs(some_range[1]-some_range[0]) * 0.05
    return (some_range[0] - margin, some_range[1] + margin)


def plot_two_parameters(method, dataframe, parameters, logs=(False, False),
                        powers=(1, 1), ranges=[(), ()], zoom=[],
                        zoomloc=[1, 2, 4], barloc=4):

    if (parameters[0] and parameters[1]) in list(ANGLES.keys()):
        color = ['#4080bf', 'Blues']
        p = 0
    else:
        color = ['green', 'Greens']
        p = 1

    before = ['', '']
    after = ['', '']
    for i in [0, 1]:
        if not ranges[i]:
            ranges[i] = DEFAULT_RANGES[parameters[i]]
        if powers[i] != 1:
            after[i] = r'$)^' + str(powers[i]) + '$'
            before[i] = r'$($'

    fig, axs = plt.subplots(
        2, 2, sharex='col', sharey='row', figsize=(10, 10),
        gridspec_kw={'height_ratios': [1, 3], 'width_ratios': [3, 1]}
        )
    fig.subplots_adjust(wspace=0, hspace=0)

    axs[0, 1].remove()

    axs[1, 0].set_xlabel(before[0] + PARAM_NAMES[parameters[0]] + after[0])
    if logs[0]:
        axs[1, 0].set_xscale('log')
    axs[1, 0].set_xlim(
        extend(ranges[0], logs[0])[0], extend(ranges[0], logs[0])[1]
        )
    axs[1, 0].set_ylabel(before[1] + PARAM_NAMES[parameters[1]] + after[1])
    if logs[1]:
        axs[1, 0].set_yscale('log')
    axs[1, 0].set_ylim(
        extend(ranges[1], logs[1])[0], extend(ranges[1], logs[1])[1]
        )
    axs[1, 0].grid()
    the_bar = inloc.inset_axes(
        axs[1, 0], width='3%', height='30%', loc=barloc, borderpad=1
        )
    sb.kdeplot(
        x=np.power(dataframe[parameters[0]].iloc[0:100000], powers[0]),
        y=np.power(dataframe[parameters[1]].iloc[0:100000], powers[1]),
        fill=True, cut=0, ax=axs[1, 0], cbar=True, cbar_ax=the_bar,
        cmap=color[1], cbar_kws={'orientation': 'vertical',
        'format': '%.2f', 'label': r'KDE'}
        )
    if barloc in [1, 4]:
        the_bar.yaxis.set_ticks_position('left')
        the_bar.yaxis.set_label_position('left')
    if parameters[0] in list(ANGLES.keys()):
        axs[1, 0].axvline(
            x=np.power(ANGLES[parameters[0]][0], powers[0]),
            color='darkred', label=r'bfp'
            )
        axs[1, 0].axvspan(
            xmin=np.power(ANGLES[parameters[0]][3], powers[0]),
            xmax=np.power(ANGLES[parameters[0]][4], powers[0]),
            color='darkred', alpha=0.3, label=r'$3\sigma$ range'
            )
    if parameters[1] in list(ANGLES.keys()):
        axs[1, 0].axhline(
            y=np.power(ANGLES[parameters[1]][0], powers[1]),
            color='darkred'
            )
        axs[1, 0].axhspan(
            ymin=np.power(ANGLES[parameters[1]][3], powers[1]),
            ymax=np.power(ANGLES[parameters[1]][4], powers[1]),
            color='darkred', alpha=0.3
            )
    if len(zoom) == 4:
        axins = inloc.inset_axes(
            axs[1, 0], width='45%', height='45%', loc=zoomloc[0], borderpad=1.5
            )
        axins.set_xlim(zoom[0], zoom[1])
        axins.set_ylim(zoom[2], zoom[3])
        sb.kdeplot(
            x=np.power(dataframe[parameters[0]].iloc[0:100000], powers[0]),
            y=np.power(dataframe[parameters[1]].iloc[0:100000], powers[1]),
            fill=True, cut=0, cmap=color[1], ax=axins
            )
        axins.set_xlabel('')
        axins.set_ylabel('')
        axins.grid()
        inloc.mark_inset(
            axs[1, 0], axins, loc1=zoomloc[1], loc2=zoomloc[2],
            fc="none", ec="0.5"
            )

    axs[0, 0].set_ylabel(r'Density of occurences')
    axs[0, 0].grid()
    if logs[0]:
        d, u = ranges[0]
        d = d//10
        u = (u+10)//10
        if d == 0:
            d = -5
        counts_tot, bins_tot = np.histogram(
            np.power(dataframe[parameters[0]][dataframe['precision'] >= p+1],
                     powers[0]),
            bins=np.logspace(d, u, 50)
            )
        counts_rest, bins_rest = np.histogram(
            np.power(dataframe[parameters[0]][dataframe['precision'] == p],
                     powers[0]),
            bins=np.logspace(d, u, 50)
            )
        axs[0, 0].hist(
            [bins_tot[:-1], bins_rest[:-1]],
            bins=np.logspace(d, u, 50), density=False, stacked=True,
            weights=[counts_tot/(np.sum(counts_tot)+np.sum(counts_rest))*9.8,
                     counts_rest/(np.sum(counts_tot)+np.sum(counts_rest))*9.8],
            color=['darkred', color[0]], alpha=0.3,
            label=['','Possible values']
            )
    else:
        axs[0, 0].hist(
            [np.power(dataframe[parameters[0]][dataframe['precision'] >= p+1],
                      powers[0]),
            np.power(dataframe[parameters[0]][dataframe['precision'] == p],
                     powers[0])],
            bins=50, density=True, stacked=True, color=['darkred', color[0]],
            alpha=0.3, label=['', 'Possible values']
            )
    sb.kdeplot(
        x=np.power(dataframe[parameters[0]], powers[0]), color=color[0], cut=0,
        ax=axs[0, 0], label=r'KDE'
        )
    axs[0, 0].set_ylim(bottom=1E-9)

    axs[1, 1].set_xlabel(r'Density of occurences')
    axs[1, 1].grid()
    if logs[1]:
        d, u = ranges[1]
        d = d//10
        u = (u+10)//10
        if d == 0:
            d = -5
        counts_tot, bins_tot = np.histogram(
            np.power(dataframe[parameters[1]][dataframe['precision'] >= p+1],
                     powers[1]),
            bins=np.logspace(d, u, 50)
            )
        counts_rest, bins_rest = np.histogram(
            np.power(dataframe[parameters[1]][dataframe['precision'] == p],
                     powers[1]),
            bins=np.logspace(d, u, 50)
            )
        axs[1, 1].hist(
            [bins_tot[:-1], bins_rest[:-1]],
            bins=np.logspace(d, u, 50), density=False, stacked=True,
            weights=[counts_tot/(np.sum(counts_tot)+np.sum(counts_rest))*9.8,
                     counts_rest/(np.sum(counts_tot)+np.sum(counts_rest))*9.8],
            color=['darkred', color[0]], alpha=0.3, orientation='horizontal'
            )
    else:
        axs[1, 1].hist(
            [np.power(dataframe[parameters[1]][dataframe['precision'] >= p+1],
                      powers[1]),
            np.power(dataframe[parameters[1]][dataframe['precision'] == p],
                     powers[1])],
            bins=50, density=True, stacked=True, color=['darkred', color[0]],
            alpha=0.3, orientation='horizontal'
            )
    sb.kdeplot(
        y=np.power(dataframe[parameters[1]], powers[1]), color=color[0], cut=0,
        ax=axs[1, 1]
        )
    axs[1, 1].set_xlim(left=1E-9)

    fig.legend(loc='lower left', bbox_to_anchor=(0.717, 0.7))

    for i in [0, 1]:
        if powers[i] != 1:
            after[i] = 'power' + str(powers[i])
            before[i] = ''
        if logs[i]:
            before[i] = 'log'

    plt.savefig(
        'figures/' + method + '_' + before[0] + parameters[0].replace('_', '')
        + after[0] + '_' + before[1] + parameters[1].replace('_', '') + after[1]
        + '.svg', format='svg', bbox_inches='tight'
        )
