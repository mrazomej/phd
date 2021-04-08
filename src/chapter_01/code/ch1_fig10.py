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

# Initialize vector of probabilities
p = np.linspace(1E-5, 1, 1000)

# compute entropy
h_p = - p * np.log2(p) - (1 - p) * np.log2((1 - p))

# Initialize figure
fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))

# Plot entropy
ax.plot(p, h_p)

# Label axis
ax.set_xlabel(r"probability $p$")
ax.set_ylabel(r"entropy $H$ (bits)")

plt.savefig("../figs/ch1_fig10A.pdf", bbox_inches="tight")
# %%

# Define threshold
epsilon = 1E-4

# Define range of m values to evaluate
m = np.arange(0, 1000)

# Define range of <m>
mean_m = np.logspace(-1, 2, 100)

# Initialize array to save entropy
h_m = np.zeros(len(mean_m))
# Loop through mean values
for i, mm in enumerate(mean_m):
    # Evaluate Poisson distribution
    pm = scipy.stats.poisson.pmf(m, mm)
    # Select values with threshold higher than threshold
    pm = pm[pm > epsilon]

    # Compute entropy
    h_m[i] = -np.sum(pm * np.log2(pm))

# %%

# Initialize figure
fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5))

# Plot entropy
ax.plot(mean_m, h_m)

# Define point to highlight
mm = 2.5
# Compute distribution
pm = scipy.stats.poisson.pmf(m, mm)
# Select values that satisfy threshold
pm = pm[pm > epsilon]
# Compute entropy
hh = -np.sum(pm * np.log2(pm))

# Highlight point
ax.scatter(mm, hh, color="#D56C55", zorder=100)

# Label axis
ax.set_xlabel(r"mean mRNA $\langle m \rangle$")
ax.set_ylabel(r"entropy $H$ (bits)")

plt.savefig("../figs/ch1_fig10B.pdf", bbox_inches="tight")
# %%

# Define mean of distribution
lam = 2.5

# Define range of mRNA
mRNA = np.arange(0, 8)

# Compute PDF
pmf = scipy.stats.poisson.pmf(mRNA, lam)

# Initialize figure
fig, ax = plt.subplots(1, 1, figsize=[1.2, 1.2])

ax.bar(mRNA, pmf, 1, tick_label=mRNA, edgecolor="#738FC1", color="#A9BFE3")

# Label plot
ax.set_xlabel(r"mRNA")
ax.set_ylabel(r"probability $P(m)$")
ax.text(3, 0.25, r"$\langle m \rangle =$" + f"{lam}", fontsize=8)

# Set background color
ax.set_facecolor("#FFEDCE")

# Set plotting limits
ax.set_ylim(top=0.3)

plt.savefig("../figs/ch1_fig10B_inset.pdf", bbox_inches="tight")
# %%

# %%
