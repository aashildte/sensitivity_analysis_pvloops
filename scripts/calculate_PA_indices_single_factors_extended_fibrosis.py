"""

Åshild Telle / University of Washington / 2025

This script was used to calculate the intersection indices between A and P loops
in the OFAT analysis, original fibrosis burden.

"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from metrics import get_loops, get_pv_loop, plot_loops, get_PA_indices


def get_indices(runs):
    indices = defaultdict(lambda: defaultdict(int))

    for run in runs:
        indices["AF5"][run] = 1

    return indices


mainfolder = "/data2/aashild/sensitivityanalysis/SA_gen2.2/extended_fibrosis"

cases = [
    "AF2",
    "AF4",
    "AF5",
]
runs = [
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
indices = get_indices(runs)

PA_indices = {}

for cas in cases:
    PA_indices[cas] = {}
    for run in runs:
        fin = f"{mainfolder}/{cas}_{run}_{run}/cav.LA.csv"
        print(cas, run)
        pressure, volume = get_pv_loop(fin)

        if cas == "AF2" and run == "baseline":
            A_index = 117
            P_index = 773
        elif cas == "AF2" and run == "CV_L":
            A_index = 117
            P_index = 772
        elif cas == "AF2" and run == "gCaL":
            A_index = 119
            P_index = 716
        elif cas == "AF2" and run == "Tax0.5":
            A_index = 116
            P_index = 762
        elif cas == "AF2" and run == "stiffness_transverse":
            A_index = 116
            P_index = 775
        elif cas == "AF4" and run == "gK1":
            A_index = 114
            P_index = 804
        elif False:
            # This code was used to find specific A_index and P_index if not found by loop intersection
            plt.plot(pressure, volume, "o-", markersize=2)
            plt.plot(pressure[A_index], volume[A_index], "o", markersize=6)
            plt.plot(pressure[P_index], volume[P_index], "o", markersize=6)
            plt.show()
            exit()
        else:
            index = indices[cas][run]
            A_index, P_index = get_PA_indices(fin, index)

        PA_indices[cas][run] = np.array([A_index, P_index])

fout = "PA_indices_single_factors_extended_fibrosis.npy"
np.save(fout, PA_indices)

# Plot for verification
fin = "PA_indices_single_factors_extended_fibrosis.npy"
PA_indices = np.load(fin, allow_pickle=True).item()

fig, axes = plt.subplots(
    len(cases), len(runs), figsize=(15, 3), sharex="row", sharey="row"
)

for i in range(len(cases)):
    axes[i][0].set_ylabel("Pressure (mmHg)")

for i in range(len(runs)):
    axes[-1][i].set_xlabel("Volume (mL)")


for cas, axs in zip(cases, axes):
    for run, axis in zip(runs, axs):
        fin = f"{mainfolder}/{cas}_{run}_{run}/cav.LA.csv"
        A_index, P_index = PA_indices[cas][run]

        P_loop_volume, P_loop_pressure, A_loop_volume, A_loop_pressure = get_loops(
            fin, A_index, P_index
        )
        plot_loops(P_loop_volume, P_loop_pressure, A_loop_volume, A_loop_pressure, axis)

for run, axis in zip(runs, axes[0]):
    axis.set_title(run)

plt.show()
