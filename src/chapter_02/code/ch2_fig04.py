# %%
import os
import glob
import pickle
import re

# Our numerical workhorses
import numpy as np
import pandas as pd
import scipy.stats

# Import the project utils
import phd

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
import matplotlib.colors as plc
import matplotlib as mpl

# Seaborn, useful for graphics
import seaborn as sns
import corner

phd.viz.pboc_style_mpl()
mpl.rcParams["figure.dpi"] = 110

# Find home directory for repo
import git

repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Find data directory
datadir = f"{homedir}/src/data/chapter_02/"
# %%

# Read the data
df = pd.read_csv(f"{datadir}flow_master.csv", comment="#")

# Remove the autofluorescence and delta values
df = df[(df.rbs != "auto") & (df.rbs != "delta")]

# Read MCMC inference of Ka Ki for O2 RBS1027
with open(f"{datadir}KaKi_mcmc.pkl", "rb") as file:
    unpickler = pickle.Unpickler(file)
    gauss_flatchain = unpickler.load()
    gauss_flatlnprobability = unpickler.load()

# map value of the parameters
max_idx = np.argmax(gauss_flatlnprobability, axis=0)
ea, ei, sigma = gauss_flatchain[max_idx]

ka_fc = np.exp(-gauss_flatchain[:, 0])
ki_fc = np.exp(-gauss_flatchain[:, 1])

# %%

# Plot the theory vs data for all 4 operators with the credible region

# Define the IPTG concentrations to evaluate
IPTG = np.logspace(-7, -2, 100)
IPTG_lin = np.array([0, 1e-7])


# Set the colors for the strains
colors = sns.color_palette("colorblind", n_colors=7)
colors[4] = sns.xkcd_palette(["dusty purple"])[0]


# Define the operators and their respective energies
operators = ["O1", "O2", "O3"]
energies = {"O1": -15.3, "O2": -13.9, "O3": -9.7, "Oid": -17}

# Initialize the figure.
fig, ax = plt.subplots(3, 3, figsize=(8.5, 7.5))
ax = ax.ravel()


# Plot the predictions.
for i, op in enumerate(operators):
    print(op)
    data = df[df.operator == op]
    # loop through RBS mutants
    for j, rbs in enumerate(df.rbs.unique()):
        # plot the theory using the parameters from the fit.
        if (op == "O2") & (rbs == "RBS1027"):
            label = None
        else:
            label = df[df.rbs == rbs].repressors.unique()[0] * 2
        # Log scale
        ax[i].plot(
            IPTG,
            phd.thermo.fold_change_log(
                IPTG * 1e6,
                ea=ea,
                ei=ei,
                epsilon=4.5,
                R=df[(df.rbs == rbs)].repressors.unique(),
                epsilon_r=energies[op],
            ),
            color=colors[j],
            label=label,
        )
        # Linear scale
        ax[i].plot(
            IPTG_lin,
            phd.thermo.fold_change_log(
                IPTG_lin * 1e6,
                ea=ea,
                ei=ei,
                epsilon=4.5,
                R=df[(df.rbs == rbs)].repressors.unique(),
                epsilon_r=energies[op],
            ),
            color=colors[j],
            linestyle=":",
            label=None,
        )

        # plot 95% HPD region using the variability in the MWC parameters
        cred_region = phd.thermo.mcmc_cred_region(
            IPTG * 1e6,
            gauss_flatchain,
            epsilon=4.5,
            R=df[(df.rbs == rbs)].repressors.unique(),
            epsilon_r=energies[op],
        )
        ax[i].fill_between(
            IPTG,
            cred_region[0, :],
            cred_region[1, :],
            alpha=0.3,
            color=colors[j],
        )  # compute the mean value for each concentration
        fc_mean = data[data.rbs == rbs].groupby("IPTG_uM").fold_change_A.mean()
        # # compute the standard error of the mean
        fc_err = data[data.rbs == rbs].groupby("IPTG_uM").fold_change_A.std() / np.sqrt(
            data[data.rbs == rbs].groupby("IPTG_uM").size()
        )

        # plot the experimental data
        # Distinguish between the fit data and the predictions
        if (op == "O2") & (rbs == "RBS1027"):
            ax[i].errorbar(
                np.sort(data[data.rbs == rbs].IPTG_uM.unique()) / 1e6,
                fc_mean,
                yerr=fc_err,
                linestyle="none",
                color=colors[j],
                label=None,
            )
            ax[i].plot(
                np.sort(data[data.rbs == rbs].IPTG_uM.unique()) / 1e6,
                fc_mean,
                marker="o",
                markersize=6,
                linestyle="none",
                markeredgewidth=1,
                markeredgecolor=colors[j],
                markerfacecolor="w",
                label=df[df.rbs == "RBS1027"].repressors.unique()[0] * 2,
            )

    # Add operator and binding energy labels.
    ax[i].set_title(
        r"%s  $\Delta\varepsilon_{RA} = %s\, k_BT$" % (op, energies[op]),
        backgroundcolor="#ffedce",
        fontsize=12,
        y=1.03,
    )
    ax[i].set_xscale("symlog", linthresh=1e-7, linscale=0.5)
    ax[i].set_xlabel("IPTG (M)", fontsize=12)
    ax[i].set_ylabel("fold-change", fontsize=12)
    ax[i].set_ylim([-0.01, 1.1])
    ax[i].set_xlim([-5e-9, 1e-2])
    ax[i].set_xticks([0, 1e-6, 1e-4, 1e-2])
    ax[i].tick_params(labelsize=10)


# Plot the properties

rep_range = np.logspace(0, 4, 200)
ka_ki = np.exp(-ea) / np.exp(-ei)
en_colors = sns.color_palette("viridis", n_colors=len(operators))
titles = [
    "leakiness",
    "saturation",
    "dynamic range",
    "EC50 ($\mu$M)",
    "effective Hill coefficient",
]
for i, op in enumerate(operators):
    # Compute the properties
    leak = phd.thermo.leakiness(rep_range, energies[op], ep_ai=4.5)
    sat = phd.thermo.saturation(rep_range, energies[op], 4.5, np.exp(-ea) / np.exp(-ei))
    dyn_rng = phd.thermo.dyn_range(rep_range, energies[op], ka_ki)
    ec50 = phd.thermo.EC50(np.exp(-ea), np.exp(-ei), 4.5, rep_range, energies[op])
    e_hill = phd.thermo.effective_Hill(np.exp(-ea), np.exp(-ei), 4.5, rep_range, energies[op])

    ax[3].plot(rep_range, leak, color=en_colors[i], label=energies[op])
    ax[4].plot(rep_range, sat, color=en_colors[i], label=energies[op])
    ax[5].plot(rep_range, dyn_rng, color=en_colors[i], label=energies[op])
    ax[6].plot(rep_range, ec50 / 1e6, color=en_colors[i])
    ax[7].plot(rep_range, e_hill, color=en_colors[i])
    ax[i + 3].set_xlabel("repressors per cell", fontsize=12)
    ax[i + 3].set_ylabel(titles[i], fontsize=12)

    # Plot the credible regions
    sat_cred = phd.thermo.saturation_cred_region(rep_range, energies[op], 4.5, ka_fc, ki_fc)
    dyn_cred = phd.thermo.dyn_cred_region(rep_range, ka_fc, ki_fc, epsilon=4.5, ep_r=energies[op])
    ec50_cred = phd.thermo.ec50_cred_region(
        rep_range, energies[op], 4.5, ka_fc, ki_fc, mass_frac=0.95
    )
    hill_cred = phd.thermo.effective_hill_cred(
        rep_range, energies[op], 4.5, ka_fc, ki_fc, mass_frac=0.95
    )
    ax[5].fill_between(
        rep_range, dyn_cred[0, :], dyn_cred[1, :], alpha=0.3, color=en_colors[i]
    )
    ax[4].fill_between(
        rep_range, sat_cred[0, :], sat_cred[1, :], alpha=0.3, color=en_colors[i]
    )
    ax[6].fill_between(
        rep_range,
        ec50_cred[0, :] / 1e6,
        ec50_cred[1, :] / 1e6,
        alpha=0.3,
        color=en_colors[i],
    )
    ax[7].fill_between(
        rep_range, hill_cred[0, :], hill_cred[1, :], alpha=0.3, color=en_colors[i]
    )

    ax[i + 3].set_xlim([1, 1e4])
    #

ax[6].set_xlim([1, 1e4])
ax[7].set_xlim([1, 1e4])
ax[6].set_ylabel("$[EC_{50}]\,\,$(M)", fontsize=12)
ax[7].set_ylabel("effective Hill coefficient", fontsize=12)
leg_1 = ax[0].legend(loc="upper left", title="rep. / cell", fontsize=8, handlelength=1)
leg_2 = ax[3].legend(
    title="   binding\n energy ($k_BT$)", loc="lower left", fontsize=8, handlelength=1
)
leg_1.get_title().set_fontsize(8)
leg_2.get_title().set_fontsize(8)
ax[3].set_yscale("log")
ax[6].set_yscale("log")
ax[6].set_yticks([1e-6, 1e-5, 1e-4])
ax[7].set_yticks([1.2, 1.4, 1.6, 1.8])
ax[8].set_axis_off()

for i in range(3, len(ax)):
    ax[i].set_xscale("log")
    ax[i].set_xticks([1, 10, 100, 1000, 1e4])
    ax[i].set_xlabel("repressors per cell", fontsize=12)

# Add plot letter label
plt.figtext(0.01, 0.96, "(C)", fontsize=12)
plt.figtext(0.33, 0.96, "(D)", fontsize=12)
plt.figtext(0.64, 0.96, "(E)", fontsize=12)
plt.figtext(0.01, 0.65, "(F)", fontsize=12)
plt.figtext(0.34, 0.65, "(G)", fontsize=12)
plt.figtext(0.64, 0.65, "(H)", fontsize=12)
plt.figtext(0.01, 0.32, "(I)", fontsize=12)
plt.figtext(0.34, 0.32, "(J)", fontsize=12)


plt.tight_layout()
plt.savefig("../figs/ch2_fig04_C-J.pdf", bbox_inches="tight")

# %%

# Generate the jointplot to insert into the figure via illustrator.
lab = ["$K_A\,\,(\mu\mathrm{M})$", "$K_I\,\,(\mu\mathrm{M})$"]
ka_ki_df = pd.DataFrame(np.array([ka_fc, ki_fc]).T, columns=lab)
inds = np.arange(0, len(ka_fc), 1)
np.random.seed(666)

# Calculate the point density

plt.close("all")
g = sns.JointGrid(
    lab[0], lab[1], ka_ki_df, xlim=(100, 200), ylim=(0.45, 0.625), space=0.05, height=2
)

g.ax_joint.plot(
    ka_fc, ki_fc, ".", color="#937D69", ms=2, alpha=0.05, rasterized=True, zorder=1
)
g.plot_joint(
    sns.kdeplot,
    cmap=sns.cubehelix_palette(n_colors=10, as_cmap=True, reverse=True),
    zorder=10,
    linewidth=1,
    n_levels=5,
    shade=True,
    alpha=0.5,
    shade_lowest=False,
)


# Plot the mode and HPD on the marginals.
ind = np.where(gauss_flatlnprobability == gauss_flatlnprobability.max())[0]
ka_mode = ka_fc[ind][0]
ki_mode = ki_fc[ind][0]
ka_cred = phd.stats.hpd(ka_fc, mass_frac=0.95)
ki_cred = phd.stats.hpd(ki_fc, mass_frac=0.95)

g.ax_marg_y.plot(6, ki_mode, "o", color=colors[4])
g.ax_marg_y.vlines(6, ki_cred[0], ki_cred[1], color=colors[4])
g.ax_marg_x.plot(ka_mode, 0.015, "o", color=colors[4])
g.ax_marg_x.hlines(0.015, ka_cred[0], ka_cred[1], color=colors[4])

g.ax_joint.set_xlim([100, 230])
g.ax_joint.set_ylim([0.45, 0.65])
g.plot_marginals(sns.kdeplot, shade=True, color=colors[4], zorder=1, linewidth=1)

# # Plot the mode and HPD for each marginal distribution
g.fig.set_figwidth(5.75)
g.fig.set_figheight(3.25)

# Save it.

# plt.tight_layout()
plt.savefig("../figs/ch2_fig04_B.pdf", bbox_inches="tight")

# %%
