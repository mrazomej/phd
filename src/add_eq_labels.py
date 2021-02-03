# %%
import numpy as np
import os
import glob
import re

# %%

# Define chapter
chapter = "chapter_05"
prefix = "ch5"

# List files
files = np.sort(glob.glob(f"{chapter}/*md"))

# %%
# Define regex pattern to find
pattern = re.compile("\$\$")

# Initialize equation counter
eq_num = 1

# Loop through files
for file in files:

    # Initialzie list to save lines
    eq_lines = list()
    
    # Open file and find lines with "$$"
    with open(file) as f:
        for i, line in enumerate(f):
            if pattern.search(line):
                eq_lines.append(i)

    # Select every other line
    eq_lines = eq_lines[1::2]

    # Open file and extract lines
    with open(file, "r") as f:
        lines = f.readlines()

    # Loop through equation lines
    for el in eq_lines:
        # Modify lines
        lines[el] = "$${#eq:" + f"{prefix}_eq" + str(eq_num).zfill(2) + "}\n"
        eq_num += 1

    # Write file
    with open(file, "w") as f:
        f.writelines(lines)
# %%
