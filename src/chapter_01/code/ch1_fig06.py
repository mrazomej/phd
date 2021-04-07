# %%
# For numerical analysis
import git
import enum
import numpy as np

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

repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir
# %%

# Define range of mRNA

# Define parameters
gamma = 1  # (a.u.)
r_gamma = 10  # mean mRNA count
r = r_gamma * gamma  # production rate

# Define time step
dt = 1e-2

# Define final time
t_final = 2.5

# Define mRNA bins
mRNA_range = np.arange(0, 25)

# Initialize matrix to save dynamics
pm = np.zeros([int(t_final / dt), len(mRNA_range)])

# Set initial condition
pm[0, 0] = 1

# Loop through time steps
for t in np.arange(1, int(t_final / dt)):
    # Loop through non-boundary conditions
    for m in mRNA_range[1:-1]:
        # Integrate master equation
        pm[t, m] = (
            pm[t - 1, m]
            + r * dt * pm[t - 1, m - 1]
            + gamma * dt * (m + 1) * pm[t - 1, m + 1]
            - r * dt * pm[t - 1, m]
            - gamma * dt * m * pm[t - 1, m]
        )
    # Deal with lower boundary
    pm[t, 0] = pm[t - 1, 0] + gamma * dt * pm[t - 1, 1] - r * dt * pm[t - 1, 0]

    # Deal with upper boundary
    pm[t, -1] = (
        pm[t - 1, m] + r * dt * pm[t - 1, -2] - gamma * dt * pm[t - 1, -1]
    )

# %%

# Initialize figure
fig, ax = plt.subplots(1, 1, figsize=[2.5, 3])

# Show probability distribution evolution heatmap
im = ax.matshow(pm, aspect="auto", cmap=cm.Blues_r)

# Remove ticks
ax.tick_params(size=4, color="white")

# Label axis
ax.set_ylabel(r"time (a.u.)")
ax.set_title(r"mRNA", size=8)

# Remove grid
ax.grid(False)

# Set ticks
xticks = [0, 4, 8, 12, 16, 20, 24]
ax.set_xticks(xticks)

# Set axis to left of plot for colorbar
cbar = ax.figure.colorbar(im, orientation="horizontal", pad=0.01)
cbar.outline.set_visible(False)
cbar.ax.set_xlabel("probability", rotation=0)#, va="bottom")
cbar.ax.tick_params(size=2, color="white")

# Define time points to plot
tpoints = [0, 25, 100, pm.shape[0]-5]

# Add reference points
ax.scatter([24] * len(tpoints), tpoints, color="white", s=1)

plt.savefig("../figs/ch1_fig06A.pdf", bbox_inches="tight")
# %%

# Initialize figure
fig, ax = plt.subplots(
    len(tpoints), 1, figsize=[2.5, 3], sharex=True, sharey=True
)

# initialize list to save axis
ax1 = list()
# Loop through time points
for i, t in enumerate(tpoints):
    ax[i].bar(
        mRNA_range,
        pm[t, :],
        1,
        tick_label=mRNA_range,
        edgecolor="#738FC1",
        color="#A9BFE3",
    )
    # Set title
    ax[i].text(10, 0.75, f"time: {t} (a.u.)", size=8)
    # Set axes limits
    ax[i].set_xlim(right=20)
    # ax[i].set_ylim(top=1)

    # Set x axis ticks
    xticks = [0, 4, 8, 12, 16, 20]
    ax[i].set_xticks(xticks)

    # Set y axis ticks on right side
    yticks = [0, 0.5, 1]

    ax[i].set_yticks([])
    ax[i].tick_params(size=4, color="white")

    ax1.append(ax[i].twinx())
    ax1[i].set_yticks(yticks)
    ax1[i].tick_params(size=4, color="white")

    

# Label plots
ax1[2].set_ylabel(r"probability $P(m, t)$")
ax[-1].set_xlabel(r"mRNA")

# Save figure
plt.savefig("../figs/ch1_fig06B.pdf", bbox_inches="tight")
# %%
