# %%
# For numerical analysis
import git
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
energies = [-10, -12.5, -15, -17.5]
repressors = np.logspace(0, np.log10(2000), 100)
Nns = 5E6

# Initialize figure
fig, ax = plt.subplots(1, 1, figsize=[2.5, 2.5])

# Define colors
colors = sns.color_palette("Greens_r", n_colors=len(energies) + 2)

# Loop through energies
for i, e in enumerate(energies):
    # Compute fold-change
    fc = (1 + repressors / Nns * np.exp(-e))**-1
    # Plot fold-change
    ax.loglog(repressors, fc, color=colors[i],
              label=f"{e}")

# Label plot
ax.set_xlabel(r"repressors/cell")
ax.set_ylabel(r"fold-change")

# Add legend
ax.legend(fontsize=7, frameon=False,
          title=r"$\Delta\epsilon_R$ ($k_BT$)", title_fontsize=7)

# Save figure
plt.savefig("../figs/ch1_fig04D.pdf", bbox_inches="tight")
# %%
