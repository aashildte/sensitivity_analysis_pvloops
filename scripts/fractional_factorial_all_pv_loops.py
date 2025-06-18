"""

Åshild Telle / University of Washington / 2025

This script extracts and plots all PV loops used for the FFD analysis.

"""

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from statannotations.Annotator import Annotator

from metrics import get_metrics, get_pv_loop, get_loops

import yaml

with open('mainfolder.yaml', 'r') as file:
    config = yaml.safe_load(file)
    mainfolder = config['input_data_org_fib']

runs = {
    "Run1": "BBBBBFFFF",
    "Run2": "FBBBBFBBB",
    "Run3": "BFBBBBFBB",
    "Run4": "FFBBBBBFF",
    "Run5": "BBFBBBBFB",
    "Run6": "FBFBBBFBF",
    "Run7": "BFFBBFBBF",
    "Run8": "FFFBBFFFB",
    "Run9": "BBBFBBBBF",
    "Run10": "FBBFBBFFB",
    "Run11": "BFBFBFBFB",
    "Run12": "FFBFBFFBF",
    "Run13": "BBFFBFFBB",
    "Run14": "FBFFBFBFF",
    "Run15": "BFFFBBFFF",
    "Run16": "FFFFBBBBB",
    "Run17": "BBBBFBBBB",
    "Run18": "FBBBFBFFF",
    "Run19": "BFBBFFBFF",
    "Run20": "FFBBFFFBB",
    "Run21": "BBFBFFFBF",
    "Run22": "FBFBFFBFB",
    "Run23": "BFFBFBFFB",
    "Run24": "FFFBFBBBF",
    "Run25": "BBBFFFFFB",
    "Run26": "FBBFFFBBF",
    "Run27": "BFBFFBFBF",
    "Run28": "FFBFFBBFB",
    "Run29": "BBFFFBBFF",
    "Run30": "FBFFFBFBB",
    "Run31": "BFFFFFBBB",
    "Run32": "FFFFFFFFF",
}

cmap = plt.cm.get_cmap("Reds", 42)

cases = ["P1", "P2", "P3"]
baseline = {"P1": {}, "P2": {}, "P3": {}}

fig, axes = plt.subplots(1, 3, figsize=(9.0, 2.5), sharey=True)

for cas, axis in zip(cases, axes):
    for key, value in runs.items():
        n = int(key[3:])
        color = cmap(np.random.randint(10, 42))
        fin = f"{mainfolder}/{cas}/{key}_{value}/cav.LA.csv"

        volume, pressure = get_pv_loop(fin)

        axis.plot(
            volume, pressure, "-", alpha=0.8, linewidth=0.8, color=color
        )  # , label=key)

for cas, axis in zip(cases, axes):
    fin = f"{mainfolder}/{cas}/baseline/cav.LA.csv"
    volume, pressure = get_pv_loop(fin)

    axis.plot(volume, pressure, "--", color="black", linewidth=2.0, label="Baseline")

# axes[-1].legend(ncols=3, bbox_to_anchor=(1.10, 1.0)) #, "Parameter changed (F)", "Parameter unchanged (B)"], bbox_to_anchor=(1.05, 1.0))

axes[0].set_xlim(85, 155)
axes[1].set_xlim(85, 155)
axes[2].set_xlim(25, 95)

for axis in axes:
    axis.set_xlabel("Volume (mL)")

axes[0].set_ylabel("Pressure (mmHg)")

axes[0].set_title("Patient 1")
axes[1].set_title("Patient 2")
axes[2].set_title("Patient 3")

plt.tight_layout()
plt.savefig("pv_loops_all_loops.png", dpi=300)
plt.show()
