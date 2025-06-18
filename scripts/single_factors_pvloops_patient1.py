"""

Åshild Telle / University of Washington / 2025

This script plots PV loops for Patient 1 simulations, for baseline versus
all factors investigated in the one factor at a time analysis.

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from collections import defaultdict
from metrics import get_metrics, get_pv_loop, get_loops

plt.rcParams.update({"mathtext.default": "regular"})

import yaml

with open('mainfolder.yaml', 'r') as file:
    config = yaml.safe_load(file)
    mainfolder = config['input_data_org_fib']

single_factors = [
    "baseline",
    "CV_L",
    "CV_T",
    "gK1",
    "gCaL",
    "gNa",
    "mux0.5",
    "Tax0.5",
    "stiffness_longitudinal",
    "stiffness_transverse",
]


def get_indices(fin):
    return np.load(fin, allow_pickle=True).item()


fin = "PA_indices_single_factors_original_fibrosis.npy"
PA_indices = get_indices(fin)

captions = [
    "baseline",
    r"Reduced $CV_L$",
    r"Reduced $CV_T$",
    r"$0.5 \times I_{K1}$",
    r"$0.5 \times I_{CaL}$",
    r"$0.6 \times I_{Na}$",
    r"$0.5 \times \mu$",
    r"$0.5 \times T_a$",
    r"2 $\times$ Load$_F$",
    r"2 $\times$ Load$_T$",
]


metrics = defaultdict(list)

fig, axes = plt.subplots(3, 3, figsize=(6.8, 5.0), sharex=True, sharey=True)

axes2 = list(axes[0]) + list(axes[1]) + list(axes[2])

cas = "P1"
subfolders = [f"{cas}/{factor}" for factor in single_factors]
fins = [f"{mainfolder}/{subfolder}/cav.LA.csv" for subfolder in subfolders]

# baseline case
volume, pressure = get_pv_loop(fins[0])

for axis in axes2:
    axis.plot(volume, pressure, "--", color="gray")
    # axis.plot(volume[60], pressure[60], "*", color="gray")

for fin, factor, caption, axis in zip(
    fins[1:], single_factors[1:], captions[1:], axes2
):
    volume, pressure = get_pv_loop(fin)
    axis.plot(volume, pressure, color="firebrick")
    # axis.plot(volume[60], pressure[60], "*", color="firebrick")
    axis.set_title(caption)

axes2[0].set_ylabel("Pressure (mmHg)")
axes2[3].set_ylabel("Pressure (mmHg)")
axes2[6].set_ylabel("Pressure (mmHg)")
axes2[2].legend(["Baseline", "Single\nfactor\nchanged"], bbox_to_anchor=(1.05, 1.0))
axes2[6].set_xlabel("Volume (ml)")
axes2[7].set_xlabel("Volume (ml)")
axes2[8].set_xlabel("Volume (ml)")
plt.tight_layout()

plt.savefig("OFAT_PV_loops.pdf")
plt.show()
