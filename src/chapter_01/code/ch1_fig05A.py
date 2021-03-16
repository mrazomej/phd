# %%
# For numerical analysis
import git
import numpy as np
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

repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir
# %%

# Define mean of distribution
lam = 2.5

# Define range of mRNA
mRNA = np.arange(0, 8)

# Compute PDF
pmf = scipy.stats.poisson.pmf(mRNA, lam)

# Initialize figure
fig, ax = plt.subplots(1, 1, figsize=[2.5, 2.5])

ax.bar(mRNA, pmf, 1, tick_label=mRNA, edgecolor="#738FC1", color="#A9BFE3")

# Label plot
ax.set_xlabel(r"mRNA")
ax.set_ylabel(r"probability $P(m, t)$")

# Set plotting limits
ax.set_ylim(top=0.3)

# Save figure
plt.savefig("../figs/ch1_fig05A.pdf", bbox_inches="tight")
# %%
