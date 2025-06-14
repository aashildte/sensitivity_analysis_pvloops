"""

Åshild Telle / University of Washington / 2025

This script plots extracts all metrics in the FFA analysis, and displayes average values in a colormap.

"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

from collections import defaultdict
from metrics import get_metrics, get_pv_loop, get_loops

plt.rcParams.update({"mathtext.default": "regular"})

cases = ["AF2", "AF4", "AF5"]

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


def get_indices(fin):
    return np.load(fin, allow_pickle=True).item()

mainfolder = f"/data2/aashild/sensitivityanalysis/SA_gen2.2"

cases = ["AF2", "AF4", "AF5"]
baseline = {"AF2": {}, "AF4": {}, "AF5": {}}

fin = "PA_indices_single_factors_original_fibrosis.npy"
PA_indices_single_factors = get_indices(fin)


for cas in cases:
    fin = f"{mainfolder}/original_fibrosis/{cas}_baseline_baseline/cav.LA.csv"
    volume, pressure = get_pv_loop(fin)
    p_start, a_start = PA_indices_single_factors[cas]["baseline"]
    _, _, A_loop_volume, A_loop_pressure = get_loops(fin, p_start, a_start)
    A_loop_area, booster, conduit, reservoir, pressure_difference = get_metrics(
        volume, pressure, A_loop_volume, A_loop_pressure
    )
    baseline[cas]["A-loop area"] = A_loop_area
    baseline[cas]["booster"] = booster
    baseline[cas]["reservoir"] = reservoir
    baseline[cas]["conduit"] = conduit
    baseline[cas]["pressure_difference"] = pressure_difference

captions = [
    "Baseline",
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


def cas2patient(cas):
    if cas == "AF2":
        return r"$P_1$"
    elif cas == "AF4":
        return r"$P_2$"
    elif cas == "AF5":
        return r"$P_3$"

for fibrosis_burden in ["original", "extended"]:

    fin = f"PA_indices_fractional_factorial_{fibrosis_burden}_fibrosis.npy"
    PA_indices_fractional_factorial_original = get_indices(fin)

    metrics = defaultdict(list)

    for cas in cases:
        subfolders = [f"{cas}_{k}_{v}_{k}_{v}" for (k, v) in runs.items()]

        fins = [f"{mainfolder}/{fibrosis_burden}_fibrosis/{subfolder}/cav.LA.csv" for subfolder in subfolders]

        for fin, run in zip(fins, runs):
            volume, pressure = get_pv_loop(fin)
            p_start, a_start = PA_indices_fractional_factorial_original[cas][run]
            _, _, A_loop_volume, A_loop_pressure = get_loops(fin, p_start, a_start)

            A_loop_area, booster, conduit, reservoir, pressure_difference = get_metrics(
                volume, pressure, A_loop_volume, A_loop_pressure
            )

            metrics["A_loop absolute"].append(A_loop_area)
            metrics["reservoir absolute"].append(reservoir)
            metrics["conduit absolute"].append(conduit)
            metrics["booster absolute"].append(booster)
            metrics["pressure_difference absolute"].append(pressure_difference)

            metrics["A_loop relative"].append(A_loop_area / baseline[cas]["A-loop area"])
            metrics["reservoir relative"].append(reservoir / baseline[cas]["reservoir"])
            metrics["conduit relative"].append(conduit / baseline[cas]["conduit"])
            metrics["booster relative"].append(booster / baseline[cas]["booster"])
            metrics["pressure_difference relative"].append(
                pressure_difference / baseline[cas]["pressure_difference"]
            )

            metrics["Run"].append(run)
            metrics["Patient"].append(cas2patient(cas))


    metrics = pd.DataFrame(data=metrics)


    metric_names = [
        "A_loop relative",
        "booster relative",
        "reservoir relative",
        "conduit relative",
        "pressure_difference relative",
    ]
    metrics_heatmap = np.zeros((32, 5))


    for metric in metric_names:
        for i, (k, v) in enumerate(runs.items()):
            print(k)
            data_caption = metrics[metrics["Run"] == k]

            # exit()
            for j, m in enumerate(metric_names):
                data_metrics = data_caption[m]
                print(m, data_metrics.mean())
                metrics_heatmap[i][j] = data_metrics.mean()


    metrics_percentage = []

    for i in range(len(metrics_heatmap)):
        metrics_percentage.append([])
        for j in range(len(metrics_heatmap[0])):
            val = metrics_heatmap[i][j]
            percent = 100 * (val - 1)
            if percent > 0:
                str_percent = f"{val:.3f}\n(+{percent:.0f}%)"
            else:
                str_percent = f"{val:.3f}\n({percent:.0f}%)"
            metrics_percentage[i].append(str_percent)

    plt.figure(figsize=(6, 20))
    colors = ["#C21E56", "#FFFFFF", "#3A5311"]
    cmap = LinearSegmentedColormap.from_list("name", colors)

    # cmap = "PiYG"
    metric_labels = [
        "A-loop\narea",
        "Booster\nfunction",
        "Reservoir\nfunction",
        "Conduit\nfunction",
        "Upstroke\npressure\ndifference",
    ]
    sns.heatmap(
        metrics_heatmap,
        cmap=cmap,
        annot=metrics_percentage,
        vmin=0.35,
        vmax=1.65,
        cbar_kws={"label": "Average change for all patients;\nrelative to baseline (-)"},
        annot_kws={"color": "black"},
        fmt="s",
    )
    plt.yticks([i + 0.5 for i in range(len(runs))], runs, rotation=0)
    plt.xticks(
        [i + 0.5 for i in range(len(metric_names))], metric_labels, rotation=0
    )


    plt.tight_layout()
    plt.savefig(
        f"heatmap_average_change_fractional_factorial_{fibrosis_burden}_fibrosis.pdf", dpi=300
    )
    plt.show()
