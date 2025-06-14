"""

Åshild Telle / University of Washington / 2025

This script plots volume, pressure, and pv loops for baseline,
impaired ICaL, imapired IK1, and fully fibrotic simulations

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# import seaborn as sns

from collections import defaultdict
from metrics import get_metrics, get_pv_loop, get_loops

plt.rcParams.update({"mathtext.default": "regular"})

mainfolder = "/data2/aashild/sensitivityanalysis/SA_gen2.2/original_fibrosis"

cases = ["AF2", "AF4", "AF5"]

single_factors = [
    "baseline",
    "gCaL",
    "gK1",
    "Run32_FFFFFFFFF",
]


def get_indices(fin):
    return np.load(fin, allow_pickle=True).item()


fin = "PA_indices_single_factors_original_fibrosis.npy"
PA_indices = get_indices(fin)

captions = [
    "baseline",
    r"$0.5 \times I_{CaL}$",
    r"$0.5 \times I_{K1}$",
    "Fully fibrotic",
]


metrics = defaultdict(list)

fig, axes = plt.subplots(3, 1, figsize=(3.0, 6.0))

cas = "AF2"
subfolders = [f"{cas}_{factor}_{factor}" for factor in single_factors]
fins = [f"{mainfolder}/{subfolder}/cav.LA.csv" for subfolder in subfolders]

# baseline case
volume, pressure = get_pv_loop(fins[0])

time = np.linspace(0, 1000, 2001)

axes[0].plot(time, volume, "--", color="black")
axes[1].plot(time, pressure, "--", color="black")
axes[2].plot(volume, pressure, "--", color="black")
print("Systole: ", np.argmin(volume)/2)

colors = ["cornflowerblue", "firebrick", "orange"]

for fin, factor, caption, color in zip(
    fins[1:], single_factors[1:], captions[1:], colors
):
    volume, pressure = get_pv_loop(fin)

    print("Systole: ", np.argmin(volume)/2)

    axes[0].plot(time, volume, color=color)
    axes[1].plot(time, pressure, color=color)
    axes[2].plot(volume, pressure, color=color)
    # axis.plot(volume[60], pressure[60], "*", color="firebrick")

axes[0].axvline(x=293 / 2, linestyle="-", color="lightgreen")

axes[0].set_xlabel("Time (ms)")
axes[1].set_xlabel("Time (ms)")
axes[2].set_xlabel("Volume (mL)")

axes[0].set_ylabel("Volume (mL)")
axes[1].set_ylabel("Pressure (mmHg)")
axes[2].set_ylabel("Pressure (mmHg)")

axes[0].legend(
    ["Baseline", r"$0.5 \times I_{CaL}$", r"$0.5 \times I_{K1}$", "Fully fibrotic"]
)

plt.tight_layout()

plt.savefig("significant_factors_volume_and_pressure_loops.pdf")
plt.show()
