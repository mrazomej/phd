# %%
# For numerical analysis
import git
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

repo = git.Repo("./", search_parent_directories=True)
homedir = repo.working_dir
# %%


def mp_sample(n_samples, mean_m=2.5, mean_p=10):
    """
    Function to sample number of mRNA and number of proteins assuming a simple
    model where m ~ Poiss(mean_m) and p/m ~ Poiss(mean_p).
    Parameters
    ----------
    mean_m : float.
        Mean number of mRNA per cell
    mean_p : float.
        Mean number of protein per mRNA per cell
    """
    # Initialize array to save samples
    mp = np.zeros([2, n_samples])

    # Sample mRNA
    mp[0, :] = np.random.poisson(mean_m, n_samples)

    # Loop through each mRNA sample and sample protein
    for i in range(n_samples):
        mp[1, i] = np.random.poisson(mp[0, i] * mean_p)

    return pd.DataFrame(mp.T, columns=["mRNA/cell", "protein/cell"])


def shannon_entropy(p, thresh=1e-5):
    """
    Function to compute the Shannon entropy
    H(p) = - Σ_i pi * log(pi),
    or the joint entropy
    H(p) = - Σ_i Σ_j pij * log(pij)
    Parameters
    ----------
    p : array-like
        vector of probabilities
    base : float. Default = 2.
        base of the logarithm to be used when computing the entropy
    thresh : float. Default = 1E-5
        Threshold for values to consider zero. This is to avoid errors when
        computing 0 x log(0).
    Returns
    -------
    H : float.
        Shannon entropy in bits
    """
    # Flatten 2D array
    p = p.ravel()

    # Remove values below threshold
    p = p[p > thresh]

    # Compute Shannon entropy
    return -np.sum(p * np.log2(p))


# %%
np.random.seed(18)

# Define number of samples
n_samples = int(1e6)

# Generate  samples of mRNA and protein
mp1 = mp_sample(n_samples, mean_m=2.5, mean_p=10)

# %%

# Plot joint distribution for toy model
sns.jointplot(
    "mRNA/cell",
    "protein/cell",
    mp1[::75],
    height=2.5,
    s=2,
    linewidth=0.05,
    marginal_kws={'hist_kws': {'linewidth': 0}}
)
plt.savefig("../figs/ch1_fig11B.pdf", bbox_inches="tight")
# %%
