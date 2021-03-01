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

# Define time range
time = np.linspace(0, 5, 100)

# Define different intial conditions
mo = np.array([0, 5, 15, 20])

# Initialize figure
fig, ax = plt.subplots(1, 1, figsize=[2.5, 2.5])

# Define colors
colors = sns.color_palette("Blues_r", n_colors=len(mo) + 2)

# Loop through initial conditions and plot
for i, m_init in enumerate(mo):
    # Compute dynamics
    m = m_init * np.exp(-gamma * time) + r_gamma * (1 - np.exp(-gamma * time))
    # Plot trajectory
    ax.plot(time, m, color=colors[i], label=f"$m_o = {m_init}$")

# Add line indicating mean mRNA
ax.axhline(r_gamma, linestyle="--", color="black",
           label="$m_{ss} = r_m / \gamma_m$")

# Label plot
ax.set_xlabel(r"time $t$ (a.u.)")
ax.set_ylabel(r"mRNA count $m(t)$")

# Add legend
ax.legend(fontsize=7, frameon=False)

# Save figure
plt.savefig("../figs/ch1_fig01C.pdf", bbox_inches="tight")
# %%
