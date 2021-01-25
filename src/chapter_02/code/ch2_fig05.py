import os
import glob
import pickle
import re

# Our numerical workhorses
import numpy as np
import pandas as pd

# Import the project utils
import sys
sys.path.insert(0, '../analysis/')
import mwc_induction_utils as mwc

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

# Seaborn, useful for graphics
import seaborn as sns

mwc.set_plotting_style()

# Load the master data file
datadir = '../../data/'
df = pd.read_csv(datadir + 'flow_master.csv', comment='#')

# Now we remove the autofluorescence and delta values
df = df[(df.rbs != 'auto') & (df.rbs != 'delta')]


# Load the flat-chain  used for parameter estimation
with open('../../data/mcmc/main_text_KaKi.pkl', 'rb') as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()

# map value of the parameters
max_idx = np.argmax(gauss_flatlnprobability, axis=0)
ea, ei, sigma = gauss_flatchain[max_idx]
ka, ki = np.exp(-ea), np.exp(-ei)

# Convert the flatchains to units of concentration.
ka_fc = np.exp(-gauss_flatchain[:, 0])
ki_fc = np.exp(-gauss_flatchain[:, 1])

# Plot the theory vs data for the titration curves
# Define the IPTG concentrations to evaluate
IPTG = np.logspace(-7, -2, 100)
IPTG_lin = np.array([0, 1E-7])

# Set the colors for the strains
colors = sns.color_palette('colorblind', n_colors=7)
colors[4] = sns.xkcd_palette(['dusty purple'])[0]
reps = [1740, 1220, 260, 124, 60, 22]
for i, r in enumerate(reps):
    if i == 0:
        rep_colors = {r: colors[i]}
    else:
        rep_colors[r] = colors[i]
# Define the operators and their respective energies
operators = ['O1', 'O2', 'O3']
energies = {'O1': -15.3, 'O2': -13.9, 'O3': -9.7, 'Oid': -17}


# Set up the gridspec figure axis.
fig = plt.figure(figsize=(8.5, 6))
gs = gridspec.GridSpec(4, 4)
ax1 = plt.subplot(gs[0:2, 0:2])
ax2 = plt.subplot(gs[0:2, 2:])
ax3 = plt.subplot(gs[2:4, 0:2])
ax4 = plt.subplot(gs[2, 2:4])
ax5 = plt.subplot(gs[3, 2:4])

# Lump  together the axes into a useable form.
ax = [ax1, ax2, ax3, ax4, ax5]

# fig, ax = plt.subplots(3, 3, figsize=(12, 10))
# ax = ax.ravel()
# Loop through operators
for i, op in enumerate(operators):
    print(op)
    data = df[df.operator == op]
    # loop through RBS mutants
    for j, rbs in enumerate(df.rbs.unique()):
        # plot the theory using the parameters from the fit.
        # Log-scale
        ax[i].plot(IPTG, mwc.fold_change_log(IPTG * 1E6,
                                             ea=ea, ei=ei, epsilon=4.5,
                                             R=np.array(
                                                 df[(df.rbs == rbs)].repressors.unique()),
                                             epsilon_r=energies[op]),
                   color=colors[j], label=None, zorder=1)
        # Linear scale
        ax[i].plot(IPTG_lin, mwc.fold_change_log(IPTG_lin * 1E6,
                                                 ea=ea, ei=ei, epsilon=4.5,
                                                 R=np.array(
                                                     df[(df.rbs == rbs)].repressors.unique()),
                                                 epsilon_r=energies[op]),
                   color=colors[j], label=None, zorder=1, linestyle=':')
        # plot 95% HPD region using the variability in the MWC parameters
        # Log scale
        cred_region = mwc.mcmc_cred_region(IPTG * 1e6,
                                           gauss_flatchain, epsilon=4.5,
                                           R=df[(df.rbs == rbs)
                                                ].repressors.unique(),
                                           epsilon_r=energies[op])
        ax[i].fill_between(IPTG, cred_region[0, :], cred_region[1, :],
                           alpha=0.3, color=colors[j])
        # Compute the mean value for each concentration
        fc_mean = data[data.rbs == rbs].groupby('IPTG_uM').fold_change_A.mean()
        # compute the standard error of the mean
        fc_err = data[data.rbs == rbs].groupby('IPTG_uM').fold_change_A.std() / \
            np.sqrt(data[data.rbs == rbs].groupby('IPTG_uM').size())

        # plot the experimental data
        # Distinguish between the fit data and the predictions
        if (op == 'O2') & (rbs == 'RBS1027'):
            ax[i].errorbar(np.sort(data[data.rbs == rbs].IPTG_uM.unique()) / 1E6,
                           fc_mean, yerr=fc_err, linestyle='none', color=colors[j])
            ax[i].plot(np.sort(data[data.rbs == rbs].IPTG_uM.unique()) / 1E6,
                       fc_mean, marker='o', linestyle='none',
                       markeredgewidth=1, markersize=6, markeredgecolor=colors[j],
                       markerfacecolor='w',
                       label=df[df.rbs == 'RBS1027'].repressors.unique()[
                0] * 2,
                zorder=100)
        else:
            ax[i].errorbar(np.sort(data[data.rbs == rbs].IPTG_uM.unique()) / 1E6,
                           fc_mean, yerr=fc_err,
                           fmt='o', label=df[df.rbs == rbs].repressors.unique()[0] * 2,
                           color=colors[j], zorder=100, markersize=6)

    # Add operator and binding energy labels.
    ax[i].set_title(r'%s  $\Delta\varepsilon_{RA} = %s\, k_BT$'
                    % (op, energies[op]), backgroundcolor='#ffedce',
                    fontsize=12, y=1.03)


# Set the sclae and labels.
    ax[i].set_xscale('symlog', linthreshx=1E-7, linscalex=0.5)
    ax[i].set_xlabel('IPTG (M)', fontsize=12)
    ax[i].set_ylabel('fold-change', fontsize=12)
    ax[i].set_ylim([-0.01, 1.1])
    ax[i].set_xlim([-5E-9, 1E-2])
    ax[i].tick_params(labelsize=12)

ax[0].legend(title='rep. / cell', loc='upper left')


# Load all of the flatchains and plot the Ka Ki values.
chains = np.sort(glob.glob('../../data/mcmc/SI_I_*.pkl'))
dfs = []
for i, c in enumerate(chains):
    with open(c, 'rb') as file:
        unpickler = pickle.Unpickler(file)
        flatchain = unpickler.load()
        flatlnprob = unpickler.load()
    max_idx = np.argmax(flatlnprob, axis=0)
    ea, ei, sigma = flatchain[max_idx]
    ka_mode = np.exp(-ea)
    ki_mode = np.exp(-ei)
    ea = flatchain[:, 0]
    ei = flatchain[:, 1]
    ka_hpd = mwc.hpd(ea, mass_frac=0.95)
    ka_hpd = np.exp(-ka_hpd)
    ki_hpd = mwc.hpd(ei, mass_frac=0.95)
    ki_hpd = np.exp(-ki_hpd)
    # Parse the file name for operator and repressor copy number.
    split = c.split('_')
    op = split[2]
    print(op)
    R = int(split[3].rstrip('.pkl')[1:])
    print(R)

    # Make a DataFrame
    df = pd.DataFrame([op, R, ka_mode, ki_mode, ka_hpd[0], ka_hpd[1],
                       ki_hpd[0], ki_hpd[1]]).T
    df.columns = ['operator', 'repressors', 'Ka_uM', 'Ki_uM', 'Ka_low',
                  'Ka_high', 'Ki_low', 'Ki_high']
    dfs.append(df)
df = pd.concat(dfs, ignore_index=True)
df = df.sort_values('repressors')
# Set up the colors and figure axis.
ops = df['operator'].unique()
en_colors = sns.color_palette('viridis', n_colors=len(ops))

color_dict = {'O1': colors[0], 'O2': colors[1], 'O3': colors[2]}
# Group the dataframe.
grouped = df.groupby(['repressors', 'operator'])

i = 1
j = 0
ff = {'O1': -0.15, 'O2': 0, 'O3': 0.15}
glyph = {'O1': 'o', 'O2': 'D', 'O3': 's'}
for o in glyph.keys():
    ax[4].plot([], [], glyph[o], markerfacecolor='w', markeredgecolor='k',
               markeredgewidth=2, ms=5.5, label=o)
ax[4].legend(loc='lower center', ncol=3)
for g, d in grouped:
    ax[3].plot(i + ff[g[1]], d['Ka_uM'], glyph[g[1]], markerfacecolor='w',
               markeredgecolor=rep_colors[g[0]], markeredgewidth=2, ms=5.5)
    ax[3].vlines(i + ff[g[1]], d['Ka_low'], d['Ka_high'], color=rep_colors[g[0]],
                 lw=1.5)

    ax[4].plot(i + ff[g[1]], d['Ki_uM'], glyph[g[1]], markerfacecolor='w',
               markeredgecolor=rep_colors[g[0]], markeredgewidth=2, ms=5.5)
    ax[4].vlines(i + ff[g[1]], d['Ki_low'], d['Ki_high'], color=rep_colors[g[0]],
                 lw=1.5)

    j += 1
    if j % len(ops) == 0:
        i += 1
repressors = df['repressors'].unique()
sel_dat = df[(df['repressors'] == 260) & (df['operator'] == 'O2')]
ka_paper = sel_dat['Ka_uM']
ka_hpd = np.linspace(sel_dat['Ka_low'], sel_dat['Ka_high'], 1000)
ki_paper = sel_dat['Ki_uM']
ki_hpd = np.linspace(sel_dat['Ki_low'], sel_dat['Ki_high'], 1000)
x_vals = np.linspace(0, i, 1000)
ax[3].hlines(ka_paper, 0.8, len(repressors) + .2,
             linestyle='--', color='k', zorder=1, alpha=0.6)
ax[4].hlines(ki_paper, 0.8, len(repressors) + .2,
             linestyle='--', color='k', zorder=1, alpha=.6)
ax[3].fill_between(np.linspace(0, i, 1000), sel_dat['Ka_low'].values[0],
                   sel_dat['Ka_high'].values[0], color='k', alpha=0.4)
ax[4].fill_between(np.linspace(0, i, 1000), sel_dat['Ki_low'].values[0],
                   sel_dat['Ki_high'].values[0], color='k', alpha=0.4)

for j, a in enumerate([ax[3], ax[4]]):
    a.set_xlim([1 - .2, (i - 1) + .2])
    a.set_xticks(np.arange(1, i, 1))
    a.set_yscale('log')
    a.set_xticks([])

# ax[3].set_xticklabels([])
# ax[4].set_xticklabels(np.sort(df['repressors'].unique()), fontsize=13)
ax[3].set_ylim([1E-3, 1E5])
ax[4].set_ylim([0.5E-3, 1E2])
ax[3].set_yticks([1E-2, 1E0, 1E2, 1E4])
ax[3].tick_params(labelsize=12)
ax[3].tick_params(labelsize=12)
ax[4].set_yticks([1E-2, 1E-1, 1E0, 1E1])
ax[3].set_ylabel('$K_A\,\,(\mu\mathrm{M})$', fontsize=12)
ax[4].set_ylabel('$K_I\,\,(\mu\mathrm{M})$', fontsize=12)
"""


## Separate the data for calculation of other properties.
#grouped = pd.groupby(df, 'operator')

##  Define functions for each property calculation.
#def dyn_range(num_rep, ep_r, ka_ki, ep_ai=4.5, n_sites=2, n_ns=4.6E6):
    #pact_leak = 1 / (1 + np.exp(-ep_ai))
    #pact_sat = 1 / (1 + np.exp(-ep_ai) * (ka_ki)**n_sites)
    #leak = (1 + pact_leak * (num_rep / n_ns) * np.exp(-ep_r))**-1
    #sat = (1 + pact_sat * (num_rep / n_ns) * np.exp(-ep_r))**-1
    #return sat - leak

#
#def dyn_cred_region(num_rep, ka_flatchain, ki_flatchain,
                    #ep_r, mass_frac=0.95, epsilon=4.5):
    #cred_region = np.zeros([2, len(num_rep)])
    ## Loop through each repressor copy number and compute the fold-changes
    ## for each concentration.
    #ka_ki = ka_flatchain / ki_flatchain
    #for i, R in enumerate(num_rep):
        #drng = dyn_range(R, ep_r, ka_ki, ep_ai=epsilon)
        #cred_region[:, i] = mwc.hpd(drng, mass_frac)
    #return cred_region

#def leakiness(num_rep, ep_r, ep_ai, n_ns=4.6E6):
    #pact = 1 / (1 + np.exp(-ep_ai))
    #return (1 + pact * (num_rep / n_ns) * np.exp(-ep_r))**-1

#def saturation(num_rep, ep_r, ep_ai, ka_ki, n_sites=2, n_ns=4.6E6):
    #pact = 1 / (1 + np.exp(-ep_ai) * ka_ki**n_sites)
    #return (1 + pact * (num_rep/n_ns)*np.exp(-ep_r))**-1

#
#def saturation_cred_region(num_rep, ep_r, ep_ai, ka_flatchain, ki_flatchain,
                           #n_sites=2, n_ns=4.6E6, mass_frac=0.95):
    #pact = 1 / (1 + np.exp(-ep_ai) * (ka_flatchain / ki_flatchain)**n_sites)
    #cred_region = np.zeros([2, len(num_rep)])
    #for i, R in enumerate(num_rep):
        #fc = (1 + pact * (R / n_ns) * np.exp(-ep_r))**-1
        #cred_region[:, i] = mwc.hpd(fc, mass_frac)
    #return cred_region

#"""
# The following equations are borrowed from Stephanie Barnes.
#"""
# def pact(IPTG, K_A, K_I, e_AI):
#'''
# Computes the probability that a repressor is active
# Parameters
#----------
#IPTG : array-like
# Array of IPTG concentrations in uM
#K_A : float
# Dissociation constant for active repressor
#K_I : float
# Dissociation constant for inactive repressor
#e_AI : float
# Energetic difference between the active and inactive state
# Returns
#-------
# probability that repressor is active
#'''
# pact = (1 + IPTG * 1 / K_A)**2 / \
#(((1 + IPTG * 1 / K_A))**2 + np.exp(-e_AI) * (1 + IPTG * 1 / K_I)**2)
# return pact

#
# def fold_change(IPTG, K_A, K_I, e_AI, R, Op):
#'''
# Computes fold-change for simple repression
# Parameters
#----------
#IPTG : array-like
# Array of IPTG concentrations in uM
#K_A : float
# Dissociation constant for active repressor
#K_I : float
# Dissociation constant for inactive repressor
#e_AI : float
# Energetic difference between the active and inactive state
#R : float
# Number of repressors per cell
#Op : float
# Operator binding energy
# Returns
#-------
# probability that repressor is active
#'''
# return 1 / (1 + R / 4.6E6 * pact(IPTG, K_A, K_I, e_AI) * np.exp(-Op))

# def EC50(K_A, K_I, e_AI, R, Op):
#'''
# Computes the concentration at which half of the repressors are in the active state
# Parameters
#----------
#K_A : float
# Dissociation constant for active repressor
#K_I : float
# Dissociation constant for inactive repressor
#e_AI : float
# Energetic difference between the active and inactive state

# Returns
#-------
# Concentration at which half of repressors are active (EC50)
#'''
#t = 1 + (R / 4.6E6) * np.exp(-Op) + (K_A / K_I)**2 * (2 * np.exp(-e_AI) + 1 + (R / 4.6E6) * np.exp(-Op))
#b = 2 * (1 + (R / 4.6E6) * np.exp(-Op)) + np.exp(-e_AI) + (K_A / K_I)**2 * np.exp(-e_AI)
# return K_A * ((K_A / K_I - 1)/(K_A / K_I - (t/b)**(1/2)) -1)

#
# def ec50_cred_region(num_rep, Op, e_AI, K_A, K_I,
# mass_frac=0.95):
#cred_region = np.zeros([2, len(num_rep)])
# for i, R in enumerate(num_rep):
#ec50_rng = EC50(K_A, K_I, e_AI, R, Op)
#cred_region[:, i] = mwc.hpd(ec50_rng, mass_frac)
# return cred_region

# def effective_Hill(K_A, K_I, e_AI, R, Op):
#'''
# Computes the effective Hill coefficient
# Parameters
#----------
#K_A : float
# Dissociation constant for active repressor
#K_I : float
# Dissociation constant for inactive repressor
#e_AI : float
# Energetic difference between the active and inactive state
# Returns
#-------
# effective Hill coefficient
#'''
#c = EC50(K_A, K_I, e_AI, R, Op)
# return 2/(fold_change(c, K_A, K_I, e_AI, R, Op) - fold_change(0, K_A, K_I, e_AI, R, Op))*\
#(-(fold_change(c, K_A, K_I, e_AI, R, Op))**2 * R / 4.6E6 * np.exp(-Op) * \
# 2 * c * np.exp(-e_AI) * (1/K_A * (1 + c/K_A) * (1 + c/K_I)**2 - 1/K_I *\
#(1 + c/K_A)**2 * (1 + c/K_I)) / ((1 + c/K_A)**2 + np.exp(-e_AI) * (1 + c/K_I)**2)**2)

# def effective_hill_cred(num_rep, Op, e_AI, K_A, K_I,
# mass_frac=0.95):
#cred_region = np.zeros([2, len(num_rep)])
# for i, R in enumerate(num_rep):
# Compute the EC50
#c = EC50(K_A, K_I, e_AI, R, Op)
# Compute the hill
# e_hill = 2/(fold_change(c, K_A, K_I, e_AI, R, Op) - fold_change(0, K_A, K_I, e_AI, R, Op))*\
#(-(fold_change(c, K_A, K_I, e_AI, R, Op))**2 * R / 4.6E6 * np.exp(-Op) * \
# 2 * c * np.exp(-e_AI) * (1/K_A * (1 + c/K_A) * (1 + c/K_I)**2 - 1/K_I *\
#(1 + c/K_A)**2 * (1 + c/K_I)) / ((1 + c/K_A)**2 + np.exp(-e_AI) * (1 + c/K_I)**2)**2)
#cred_region[:, i] = mwc.hpd(e_hill, mass_frac)

# return cred_region

# Compute the dynamic range
#drs = []
# for g, d in grouped:
#unique_IPTG = d.IPTG_uM.unique()
#min_IPTG = np.min(unique_IPTG)
#max_IPTG = np.max(unique_IPTG)
# Group the new data by repressors.
#grouped_rep = pd.groupby(d, ['rbs', 'date', 'username'])
# rbs_ind = {'HG104' : 0, 'RBS1147': 1, 'RBS446' : 2, 'RBS1027': 3,
#'RBS1': 4, 'RBS1L': 5}
#rep_dr = [[], [], [], [], [], []]
#rep_std = []
# for g_rep, d_rep in grouped_rep:
# if g_rep[2] != 'sloosbarnes':
#dr = d_rep[d_rep.IPTG_uM==max_IPTG].fold_change_A.values - d_rep[d_rep.IPTG_uM==min_IPTG].fold_change_A.values
# rep_dr[rbs_ind[g_rep[0]]].append(dr[0])

# Compute the means.
# for i, dr in enumerate(rep_dr):
#rep_dr[i] = np.mean(dr)
#rep_std.append(np.std(dr) / np.sqrt(len(dr)))

#reps = np.sort(df.repressors.unique())
#dr_df = pd.DataFrame([reps, rep_dr, rep_std]).T
#dr_df.columns = ['repressors', 'dynamic_range', 'err']
#dr_df.insert(0, 'operator', g)
# drs.append(dr_df)
#drng = pd.concat(drs, axis=0)

#
# Load in the flatchains for the calculation of the effective hill and EC50
#repressors = ['R22', 'R60', 'R124', 'R260', 'R1220', 'R1740']
#flatchains = [[], [], []]
#kas = [[], [], []]
#kis = [[], [], []]
# for i, op in enumerate(operators):
# for j, R in enumerate(repressors):
# with open('../../data/mcmc/SI_I_' + op + '_' + R + '.pkl', 'rb') as file:
# print(j)
#unpickler = pickle.Unpickler(file)
#gauss_flatchain = unpickler.load()
# flatchains[i].append(gauss_flatchain)
#gauss_flatlnprobability = unpickler.load()
#ind = np.argmax(gauss_flatlnprobability)
#kas[i].append(np.exp(-gauss_flatchain[ind, 0]))
#kis[i].append(np.exp(-gauss_flatchain[ind, 1]))

#
# Plot the leakiness, saturation, dynamic range, effective hill, and ec50.
#rep_range = np.logspace(0, 4, 200)
#ka, ki = np.exp(-ea) , np.exp(-ei)
#ka_ki = ka / ki

#en_colors = sns.color_palette('viridis', n_colors=len(operators))
#repressor_numbers = [22, 60, 124, 260, 1220, 1740]
# for i, op in enumerate(operators):
# Compute the dynamic range.
#sat = saturation(rep_range, energies[op], 4.5, np.exp(-ea)/np.exp(-ei))
#leak = leakiness(rep_range, energies[op], 4.5)
#dyn_rng = dyn_range(rep_range, energies[op], ka_ki)
#ec50 = EC50(ka/1E6, ki/1E6, 4.5, rep_range, energies[op])
#hill = effective_Hill(ka, ki, 4.5, rep_range, energies[op])
# ax[3].plot(rep_range, leak, color=en_colors[i], label='__nolegend__',
# markersize=6)
# ax[4].plot(rep_range, sat, color=en_colors[i], label='__nolegend__',
# markersize=6)
# ax[5].plot(rep_range, dyn_rng, color=en_colors[i], label='__nolegend__',
# markersize=6)
#ax[6].plot(rep_range, ec50, color=en_colors[i], label='__nolegend__')
#ax[7].plot(rep_range, hill, color=en_colors[i], label='__nolegend__')
# Compute the credible regions.
# sat_cred = saturation_cred_region(rep_range, energies[op], 4.5, ka_fc,
# ki_fc)
# dyn_cred = dyn_cred_region(rep_range,
#ka_fc, ki_fc, epsilon=4.5,
# ep_r=energies[op])
#ec50_cred = ec50_cred_region(rep_range, energies[op], 4.5, ka_fc/1E6, ki_fc/1E6)
#hill_cred = effective_hill_cred(rep_range, energies[op], 4.5, ka_fc, ki_fc)
# ax[4].fill_between(rep_range, sat_cred[0,:], sat_cred[1,:],
# alpha=0.3, color=en_colors[i])
# ax[5].fill_between(rep_range, dyn_cred[0,:], dyn_cred[1,:],
# alpha=0.3, color=en_colors[i])
# ax[6].fill_between(rep_range, ec50_cred[0,:], ec50_cred[1,:],
# alpha=0.3, color=en_colors[i])
# ax[7].fill_between(rep_range, hill_cred[0,:], hill_cred[1,:],
# alpha=0.3, color=en_colors[i])

# Plot the inferred parameter values for EC50 and effective hill
# for j, R in enumerate(repressor_numbers):
#ec50_inf = EC50(kas[i][j], kis[i][j], 4.5, R, energies[op])
# hill_inf = effective_Hill(kas[i][j], kis[i][j], 4.5, R,
# energies[op])
# convert the flatchains to units of concentration
#_ka_fc = np.exp(-flatchains[i][j][:,0])
#_ki_fc = np.exp(-flatchains[i][j][:,1])
#ec50_cred = EC50(_ka_fc, _ki_fc, 4.5, R, energies[op])
#ec50_cred = mwc.hpd(ec50_cred, mass_frac=0.95)
#hill_cred = effective_Hill(_ka_fc, _ki_fc, 4.5, R, energies[op])
#hill_cred = mwc.hpd(hill_cred, 0.95)

#
# if j == 0:
#label = energies[op]
# else:
#label = '__nolegend__'

#ax[6].vlines(R, ec50_cred[0]/1E6, ec50_cred[1]/1E6, color=en_colors[i], zorder=4-i, label='__nolegend__')
#ax[7].vlines(R, hill_cred[0], hill_cred[1], color=en_colors[i], zorder=4-i, label='__nolegend__')
#ax[6].plot(R, ec50_inf/1E6, 's', markerfacecolor='w', markeredgecolor=en_colors[i], ms=6, markeredgewidth=1.5, zorder=4-i)
#ax[7].plot(R, hill_inf, 's', ms=6,markerfacecolor='w', markeredgecolor=en_colors[i], markeredgewidth=1.5, zorder=4-i)

#
#
#
#
# Get the dynamic range data and plot.
# for i, op in enumerate(operators):
#dyn_rng = drng[drng.operator==op]
# ax[5].errorbar(2 * dyn_rng.repressors, dyn_rng.dynamic_range, yerr=dyn_rng.err, color=en_colors[i], fmt='o', linestyle='none',
# label=energies[op], markersize=5)

# Plot the leakiness and saturation data.
#grouped = pd.groupby(df, ['operator', 'repressors'])
# Define the colors so I don't have to make another data frame.
#op_colors = {'O1': en_colors[0], 'O2': en_colors[1], 'O3': en_colors[2]}
# for g, d in grouped:
# Extract the unique IPTG values.
#unique_IPTG = d['IPTG_uM'].unique()

# Slice the min and max IPTG values.
#leak_vals = d[d['IPTG_uM'] == np.min(unique_IPTG)].fold_change_A
#sat_vals = d[d['IPTG_uM'] == np.max(unique_IPTG)].fold_change_A

# Compute the mean and standard errors of reach.
#mean_leak = np.mean(leak_vals)
#sem_leak = np.std(leak_vals) / np.sqrt(len(leak_vals))
#mean_sat = np.mean(sat_vals)
#sem_sat =  np.std(sat_vals) / np.sqrt(len(sat_vals))

# Plot the data with the appropriate legends..
# if g[1] == 11:
#legend = energies[g[0]]
# else:
#legend = '__nolegend__'
# ax[3].plot(2 * g[1], mean_leak, 'o', color=op_colors[g[0]],
# markersize=5, label='__nolegend__')
# ax[4].plot(2 * g[1], mean_sat, 'o', color=op_colors[g[0]],
# markersize=5, label='__nolegend__')
# ax[3].errorbar(2 * g[1], mean_leak, sem_leak, linestyle='none',
# color=op_colors[g[0]], fmt='o', markersize=6, label=legend)
# ax[4].errorbar(2 * g[1], mean_sat, sem_sat, linestyle='none',
# color=op_colors[g[0]], fmt='o', markersize=6, label=legend)

# Add labels and format axes.
#ylabels = ['leakiness', 'saturation', 'dynamic range', '$[EC_{50}]\, \,(\mathrm{M})$', 'effective Hill coefficient']
# for i in range(len(ylabels)):
#ax[i+3].set_xlabel('repressors per cell', fontsize=14)
# ax[i+3].set_xscale('log')
#ax[i+3].set_ylabel(ylabels[i], fontsize=14)
# ax[3].set_yscale('log')
# ax[6].set_yscale('log')
# ax[8].set_axis_off()

# for i, a in enumerate(ax):
# if i < 3:
#a.set_xlim([-5E-9, 1E-2])
#a.set_xticks([0, 1E-6, 1E-4, 1E-2])
# else:
#a.set_xlim([1, 1E4])
#a.set_xticks([1, 10, 100, 1000, 10000])
# a.tick_params(labelsize=14)
# if i < 5:
#a.set_ylim([-0.01, 1.1])
#l = ax[3].legend(loc='lower left', title='binding\n energy ($k_BT$)')
#plt.setp(l.get_title(), multialignment='center')
#ax[0].legend(loc='upper left', title='rep. / cell')
# plt.tight_layout()
# plt.show()

# add plot letter labels
#plt.figtext(0., .96, 'A', fontsize=20)
#plt.figtext(0.33, .96, 'B', fontsize=20)
#plt.figtext(0.65, .96, 'C', fontsize=20)
#plt.figtext(0.0, .63, 'D', fontsize=20)
#plt.figtext(0.33, .63, 'E', fontsize=20)
#plt.figtext(0.65, .63, 'F', fontsize=20)
#plt.figtext(0.0, .32, 'G', fontsize=20)
#plt.figtext(0.33, .32, 'H', fontsize=20)

# plt.tight_layout()
# plt.savefig(output + '/fig5.pdf',
# bbox_inches='tight')
plt.savefig('../../figures/main_figs/fig6.svg', bbox_inches='tight')
