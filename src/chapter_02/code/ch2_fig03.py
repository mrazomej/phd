#%%
# For interaction with the operating system.
import glob

# For scientific computing
import numpy as np
import pandas as pd
import scipy.stats

# Import the phd utils
import phd

# Import matplotlib stuff for plotting
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker
import matplotlib as mpl

# Seaborn, useful for graphics
import seaborn as sns

phd.viz.pboc_style_mpl()
mpl.rcParams["figure.dpi"] = 110

# Find home directory for repo
import git

repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir

# Find data directory
datadir = f"{homedir}/src/data/chapter_02/"


#%%

# Purpose is to generate figures of distributions from flow cytometry for
# a flowchart diagram.
flow_glob = glob.glob(f"{datadir}/flow_csv/201608*O2*RBS1027*.csv")
cell_cloud = pd.read_csv(flow_glob[3], comment="#")

# Fit a two-dimensional gaussian to the data.
mu, cov = phd.flow.fit_2D_gaussian(cell_cloud, log=True)

# Compute the statistic for each of the pair of log scattering data
interval_array = phd.flow.gauss_interval(cell_cloud, mu, cov, log=True)
alpha = 0.40

# Find which data points fall inside the interval
idx = interval_array <= scipy.stats.chi2.ppf(alpha, 2)
disc_idx = interval_array >= scipy.stats.chi2.ppf(alpha, 2)
# Select those data which lay within the 40th percentile.
selected_cells = cell_cloud[idx]
discarded_cells = cell_cloud[disc_idx]

#%%

# Now just generate the plot.
plt.figure()
plt.plot(
    discarded_cells["FSC-A"],
    discarded_cells["SSC-A"],
    "k.",
    rasterized=True,
    label="discarded cells",
    markersize=3,
)
plt.plot(
    selected_cells["FSC-A"],
    selected_cells["SSC-A"],
    "m.",
    rasterized=True,
    label="selected cells",
    markersize=3,
)

# Fix formatting and restrict bounds.
plt.legend(loc="lower right", fontsize=6, markerscale=2)
plt.xlabel("forward scatter (a.u.)", fontsize=8)
plt.ylabel("side scatter (a.u.)", fontsize=8)
plt.xscale("log")
plt.yscale("log")
plt.tick_params(labelsize=8)
# Restrict bounds for aesthetic reasons.
plt.xlim([1e3, 1e5])
plt.ylim([1e3, 3e5])

# Save the figure.
plt.savefig("../figs/ch2_fig03_flow_cloud.pdf", bbox_inches="tight")
plt.close()

#%%
fig, ax = plt.subplots(2, 1, figsize=(6, 4), sharex=True)
# Now generate the example distributions.
colors_RBS1027 = sns.color_palette("Blues", n_colors=6)
colors_delta = sns.color_palette("Greens", n_colors=6)
IPTG_range = [0, 25, 50, 100, 500, 5000]
ax[0].plot(
    [],
    [],
    "v",
    markersize=10,
    markeredgecolor="k",
    markeredgewidth=1,
    markerfacecolor=colors_RBS1027[-1],
    label="mean",
)
ax[1].plot(
    [],
    [],
    "v",
    markersize=10,
    markeredgecolor="k",
    markeredgewidth=1,
    markerfacecolor=colors_delta[-1],
    label="mean",
)
ax[0].legend(title=r"repressors / cell = 260")
ax[1].legend(title=r"repressors / cell = 0")
for i, val in enumerate(IPTG_range):
    glob_RBS1027 = glob.glob(
        f"{datadir}/flow_csv/201608*O2*RBS1027*_%suM*IPTG.csv" % val
    )
    data_RBS1027 = pd.read_csv(glob_RBS1027[0])
    glob_delta = glob.glob(f"{datadir}/flow_csv/201608*O2*delta*_%suM*IPTG.csv" % val)
    data_delta = pd.read_csv(glob_delta[0])
    # Fit a two-dimensional gaussian to the data.
    gate_RBS1027 = phd.flow.gaussian_gate(data_RBS1027, alpha)
    gate_delta = phd.flow.gaussian_gate(data_delta, alpha)
    ax[0].hist(
        gate_RBS1027["FITC-A"],
        color=colors_RBS1027[i],
        alpha=0.5,
        bins=100,
        histtype="stepfilled",
        density=True,
    )
    mean_RBS1027 = np.mean(gate_RBS1027["FITC-A"])
    ax[0].plot(
        mean_RBS1027,
        1.3e-4,
        "v",
        markeredgecolor="k",
        markeredgewidth=1,
        markerfacecolor=colors_RBS1027[i],
        markersize=10,
    )
    ax[1].hist(
        gate_delta["FITC-A"],
        color=colors_delta[i],
        alpha=0.5,
        bins=100,
        histtype="stepfilled",
        density=True,
    )
    mean_delta = np.mean(gate_delta["FITC-A"])
    ax[1].plot(
        mean_delta,
        7.5e-5,
        "v",
        markeredgecolor="k",
        markeredgewidth=1,
        markerfacecolor=colors_delta[i],
        markersize=10,
    )

ax[0].yaxis.get_major_formatter().set_powerlimits((0, -1))
ax[1].yaxis.get_major_formatter().set_powerlimits((0, -1))
ax[0].xaxis.get_major_formatter().set_powerlimits((0, -1))
ax[1].xaxis.get_major_formatter().set_powerlimits((0, -1))

fig.text(0, 0.5, "frequency", fontsize=18, rotation="vertical")
ax[1].set_xlabel("total cell intensity (a.u.)", fontsize=18)

plt.savefig("../figs/ch2_fig03_flow_distributions.pdf")
plt.close()
#%%

# Now plot a full titration from this set.
# Read master CSV file
titration_data = pd.read_csv(f"{datadir}flow_master.csv")
# Keep only date data
titration_data = titration_data[
    (titration_data.date == 20160813) & (titration_data.rbs == "RBS1027")
].iloc[12:, :]

# Add "fake plot for legend"
plt.plot(
    [],
    [],
    "o",
    markersize=5,
    markeredgecolor="r",
    markerfacecolor="w",
    markeredgewidth=1,
    label="experimental data",
)
leg = plt.legend(
    loc="upper left",
    title=r"""repressors / cell = 260
        $\Delta\varepsilon_{RA} = -13.9\,k_BT$""",
    fontsize=7,
)
leg.get_title().set_fontsize("8")
plt.plot(
    titration_data.IPTG_uM / 1e6,
    titration_data.fold_change_A,
    "o",
    markeredgecolor="r",
    markerfacecolor="w",
    markeredgewidth=1,
    markersize=5,
)
plt.xscale("log")
plt.ylabel("fold-change", fontsize=8)
plt.xlabel("[IPTG] (M)", fontsize=8)
plt.xlim([1e-8, 1e-2])
plt.ylim([-0.01, 1.1])
plt.tick_params(labelsize=8)
plt.tight_layout()
plt.savefig("../figs/ch2_fig03_titration.pdf")
plt.close()

# %%
