# %%
import numpy as np
import os
import glob
import re

# %%

# Define chapter
chapter = "chapter_04"
prefix = "ch4"

# List files
files = np.sort(glob.glob(f"{chapter}/*md"))
# %%

# Define pattern to be found
# NOTE: (\S+) means match patterns with no spaces
ref_pattern = re.compile("\[\@Eq\:(\S+)\]") #\]")

# Initialize dictionary to save labels per file
f_labels = dict()

# Initialize global list of paterns
ref_list = []
# Loop through files
for f in files:
    print(f)
    # Open file and extract lines
    with open(f, "r") as file:
        lines = file.readlines()

    # Initialize dictionary to save lines with matching pattern
    eq_labels = dict()

    # Loop through lines and find the pattern
    for i, l in enumerate(lines):
        # Search pattern
        p = ref_pattern.findall(l)
        # If found pattern, store index and pattern
        if p:
            # Add labels to dictionary
            eq_labels[i] = p
            # Save labels to list
            ref_list = ref_list + p
    # Save file patterns
    f_labels[f] = eq_labels

# Keep unique list of labels
ref_list = list(set(ref_list))
# %%

# Loop through files to modify patterns
for j, f in enumerate(files):
    print(f)
    # Extract list of patterns in file
    f_pattern = f_labels[f]
    # check that the file doesn't have entries to skip it
    if len(f_pattern) == 0:
        continue

    # Open file and extract lines
    with open(f, "r") as file:
        lines = file.readlines()
    # Loop through dictionary with labels to modify pattern
    for key, value in f_pattern.items():
        # Loop through patterns
        for p in value:
            lines[key] = re.sub(
                f"\[\@Eq\:{p}\]", r"Eq. $\\ref{eq:" + p + "}$", lines[key]
            )

    # Write file
    with open(f, "w", encoding='utf-8') as file:
        file.writelines(lines)

# Establish pattern to modify labels
eq_pattern = re.compile("\$\$\{\#(\S+)\}")

# Loop through files
for f in files:
    print(f)
    # Open file and extract lines
    with open(f, "r") as file:
        lines = file.readlines()

    # Loop through lines and find the pattern
    for i, l in enumerate(lines):
        # Search pattern
        p = eq_pattern.findall(l)
        # If found pattern, store index and pattern
        if len(p) > 0: 
            lines[i] = r"\label{" + p[0] + "}\n$$\n"

    # Write file
    with open(f, "w", encoding='utf-8') as file:
        file.writelines(lines)


# %%