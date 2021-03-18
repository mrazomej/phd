# %%
# For numerical analysis
import enum
import git
import numpy as np

# For text
import re
import string
from numpy.core.function_base import linspace

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

# Define last paragraph of on the origin

text = "Thus, from the war of nature, from famine and death, the most exalted object of which we are capable of conceiving, namely, the production of the higher animals, directly follows. There is grandeur in this view of life, with its several powers, having been originally breathed into a few forms or into one; and that, whilst this planet has gone cycling on according to the fixed law of gravity, from so simple a beginning endless forms most beautiful and most wonderful have been, and are being, evolved."

# turn to lower case
text = text.lower()

# Remove punctuation
text = re.sub("\,", "", text)
text = re.sub("\.", "", text)
text = re.sub("\;", "", text)

# List alphabet
abc = list(string.ascii_lowercase + " ")

# Generate list for counts
abc_counts = np.zeros(len(abc))
# Loop through characters and count
for i, a in enumerate(abc):
    abc_counts[i] = text.count(a)

# Normalize count to probability
abc_counts = abc_counts / abc_counts.sum()

# Generate matrix for pairs
pair_counts = np.zeros([len(abc), len(abc)])
# Loop through first character
for i, a in enumerate(abc):
    # Loop through second character
    for j, b in enumerate(abc):
        pair_counts[i, j] = text.count(a + b)

# Normalize pair counts
pair_counts = pair_counts / pair_counts.sum()
# %%
# Define color
color = sns.color_palette("viridis", n_colors=2)[0]

# Initialize figure
fig = plt.figure(figsize=(0.5, 3.5))

# Plot frequency with size of squares representing frequency
plt.scatter(x=np.zeros_like(abc_counts), y=range(len(abc_counts)),
            s=abc_counts * 200, marker="s", color=color)

# Set y axis ticks as letters
plt.yticks(range(len(abc_counts)), abc)

# Turn off x axis
plt.xticks([])

# Turn off grid
plt.grid(b=None)

plt.savefig("../figs/ch1_fig08B.pdf", bbox_inches="tight")
# %%

# Initialize figure
fig = plt.figure(figsize=(3.5, 3.5))

# Loop through each row
for i, a in enumerate(abc):
    # Plot frequency with size of squares representing frequency
    plt.scatter(x=range(len(abc)), y=np.ones(len(abc)) * i,
                s=pair_counts[:, i] * 1000, marker="s", color=color)

# Set  axis ticks as letters
plt.xticks(range(len(abc_counts)), abc)
plt.yticks(range(len(abc_counts)), abc)

plt.xlabel("first letter")
plt.ylabel("second letter")

# Turn off grid
plt.grid(b=None)
plt.savefig("../figs/ch1_fig08C.pdf", bbox_inches="tight")
# %%
