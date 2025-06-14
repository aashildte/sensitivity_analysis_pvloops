"""

Åshild Telle / University of Washington / 2025

This script was used to calculate the intersection indices between A and P loops
in the FFD analysis, original fibrosis burden.

"""

import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

from metrics import get_loops, get_pv_loop, plot_loops, get_PA_indices


def get_indices_fractional_factorial_org_fibrosis():
    indices = defaultdict(lambda: defaultdict(int))

    for i in range(1, 33):
        indices["AF5"][f"Run{i}"] = 1
    return indices


indices = get_indices_fractional_factorial_org_fibrosis()
mainfolder = "/data2/aashild/sensitivityanalysis/SA_gen2.2/original_fibrosis"

cases = [
    "AF2",
    "AF4",
    "AF5",
]
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

PA_indices = {}

import sys

r1 = sys.argv[1]  # "Run20"
c1 = "AF2"

for cas in cases:
    PA_indices[cas] = {}
    for key, value in runs.items():
        # if cas != c1 or key != r1:
        #    continue

        fin = f"{mainfolder}/{cas}_{key}_{value}_{key}_{value}/cav.LA.csv"
        pressure, volume = get_pv_loop(fin)

        if key == "Run1" and cas == "AF2":
            A_index = 116
            P_index = 765
        elif key == "Run3" and cas == "AF2":
            A_index = 116
            P_index = 764
        elif key == "Run4" and cas == "AF2":
            A_index = 116
            P_index = 772
        elif key == "Run5" and cas == "AF2":
            A_index = 116
            P_index = 783
        elif key == "Run6" and cas == "AF2":
            A_index = 116
            P_index = 775
        elif key == "Run8" and cas == "AF2":
            A_index = 116
            P_index = 771
        elif key == "Run9" and cas == "AF2":
            A_index = 119
            P_index = 730
        elif key == "Run10" and cas == "AF2":
            A_index = 119
            P_index = 729
        elif key == "Run11" and cas == "AF2":
            A_index = 120
            P_index = 730
        elif key == "Run12" and cas == "AF2":
            A_index = 119
            P_index = 729
        elif key == "Run13" and cas == "AF2":
            A_index = 118
            P_index = 741
        elif key == "Run14" and cas == "AF2":
            A_index = 117
            P_index = 746
        elif key == "Run15" and cas == "AF2":
            A_index = 117
            P_index = 745
        elif key == "Run16" and cas == "AF2":
            A_index = 118
            P_index = 743
        elif key == "Run17" and cas == "AF2":
            A_index = 117
            P_index = 770
        elif key == "Run18" and cas == "AF2":
            A_index = 116
            P_index = 764
        elif key == "Run21" and cas == "AF2":
            A_index = 116
            P_index = 773
        elif key == "Run23" and cas == "AF2":
            A_index = 115
            P_index = 771
        elif key == "Run25" and cas == "AF2":
            A_index = 119
            P_index = 728
        elif key == "Run26" and cas == "AF2":
            A_index = 119
            P_index = 728
        elif key == "Run27" and cas == "AF2":
            A_index = 119
            P_index = 728
        elif key == "Run28" and cas == "AF2":
            A_index = 119
            P_index = 728
        elif key == "Run29" and cas == "AF2":
            A_index = 117
            P_index = 744
        elif key == "Run30" and cas == "AF2":
            A_index = 118
            P_index = 739
        elif key == "Run32" and cas == "AF2":
            A_index = 117
            P_index = 743
        elif key == "Run4" and cas == "AF4":
            A_index = 116
            P_index = 790
        elif key == "Run7" and cas == "AF4":
            A_index = 117
            P_index = 802
        elif key == "Run19" and cas == "AF4":
            A_index = 117
            P_index = 787
        elif key == "Run22" and cas == "AF4":
            A_index = 117
            P_index = 796
        elif key == "Run25" and cas == "AF4":
            A_index = 119
            P_index = 737
        elif key == "Run26" and cas == "AF4":
            A_index = 120
            P_index = 738
        elif key == "Run28" and cas == "AF4":
            A_index = 120
            P_index = 738
        elif False:
            # This code was used to find specific A_index and P_index if not found by loop intersection
            plt.plot(pressure, volume, "o-", markersize=2)
            plt.plot(pressure[A_index], volume[A_index], "o", markersize=6)
            plt.plot(pressure[P_index], volume[P_index], "o", markersize=6)
            plt.show()
            exit()
        else:
            index = indices[cas][key]
            A_index, P_index = get_PA_indices(fin, index)

        PA_indices[cas][key] = np.array([A_index, P_index])

fout = "PA_indices_fractional_factorial_original_fibrosis.npy"
np.save(fout, PA_indices)

fin = "PA_indices_fractional_factorial_original_fibrosis.npy"
PA_indices = np.load(fin, allow_pickle=True).item()

fig, axes = plt.subplots(
    len(cases), len(runs), figsize=(20, 7), sharex="row", sharey="row"
)

for i in range(len(cases)):
    axes[i][0].set_ylabel("Pressure (mmHg)")

for i in range(len(runs)):
    axes[-1][i].set_xlabel("Volume (mL)")

captions = [run for run in runs.keys()]

for cas, axs in zip(cases, axes):
    baseline = f"{mainfolder}/{cas}_baseline_baseline/cav.LA.csv"

    volume, pressure = get_pv_loop(baseline)

    # for axis in axs:
    #    axis.plot(volume, pressure, color="gray")

    for (key, value), axis in zip(runs.items(), axs):
        # if cas != c1 or key != r1:
        #    continue

        fin = f"{mainfolder}/{cas}_{key}_{value}_{key}_{value}/cav.LA.csv"
        A_index, P_index = PA_indices[cas][key]

        P_loop_volume, P_loop_pressure, A_loop_volume, A_loop_pressure = get_loops(
            fin, A_index, P_index
        )
        plot_loops(P_loop_volume, P_loop_pressure, A_loop_volume, A_loop_pressure, axis)
        # except IndexError:
        #    pressure, volume = get_pv_loop(fin)
        #    axis.plot(volume, pressure)

for key, axis in zip(runs, axes[0]):
    axis.set_title(key)

plt.show()
