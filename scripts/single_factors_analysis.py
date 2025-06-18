"""

Åshild Telle / University of Washington / 2025

This script plots
- results from the one factor at a time analysis, for each patient individually.
- average values (across all patients) in a heatmap

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

import yaml

with open('mainfolder.yaml', 'r') as file:
    config = yaml.safe_load(file)
    mainfolder = config['input_data_org_fib']

cases = ["P1", "P2", "P3"]

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
    if cas == "P1":
        return r"$P_1$"
    elif cas == "P2":
        return r"$P_2$"
    elif cas == "P3":
        return r"$P_3$"


metrics = defaultdict(list)

for cas in cases:
    subfolders = [f"{cas}/{factor}" for factor in single_factors]

    fins = [f"{mainfolder}/{subfolder}/cav.LA.csv" for subfolder in subfolders]

    for fin, factor, caption in zip(fins, single_factors, captions):
        volume, pressure = get_pv_loop(fin)
        p_start, a_start = PA_indices[cas][factor]
        _, _, A_loop_volume, A_loop_pressure = get_loops(fin, p_start, a_start)

        A_loop_area, booster, conduit, reservoir, pressure_difference = get_metrics(
            volume, pressure, A_loop_volume, A_loop_pressure
        )
        if caption == "Baseline":
            baseline_A_loop_area = A_loop_area
            baseline_reservoir = reservoir
            baseline_conduit = conduit
            baseline_booster = booster
            baseline_pressure_difference = pressure_difference
            continue

        metrics["A_loop absolute"].append(A_loop_area)
        metrics["reservoir absolute"].append(reservoir)
        metrics["conduit absolute"].append(conduit)
        metrics["booster absolute"].append(booster)
        metrics["pressure_difference absolute"].append(pressure_difference)

        metrics["A_loop relative"].append(A_loop_area / baseline_A_loop_area)
        metrics["booster relative"].append(booster / baseline_booster)
        metrics["reservoir relative"].append(reservoir / baseline_reservoir)
        metrics["conduit relative"].append(conduit / baseline_conduit)
        metrics["pressure_difference relative"].append(
            pressure_difference / baseline_pressure_difference
        )

        metrics["Caption"].append(caption)
        metrics["Patient"].append(cas2patient(cas))

metrics = pd.DataFrame(data=metrics)

fig, axes = plt.subplots(5, 1, sharex=True, figsize=(4.4, 9.0))

colors_AJ = ["salmon", "tab:red", "maroon"]
colors_booster = ["navajowhite", "tab:orange", "xkcd:rust"]
colors_reservoir = ["palegreen", "tab:green", "darkolivegreen"]
colors_conduit = ["lightblue", "royalblue", "navy"]
colors_pressure_difference = ["violet", "darkorchid", "indigo"]

metric_names = [
    "A_loop relative",
    "booster relative",
    "reservoir relative",
    "conduit relative",
    "pressure_difference relative",
]
colors = ["tab:red", "tab:orange", "tab:green", "royalblue", "indigo"]
colors = [
    colors_AJ,
    colors_booster,
    colors_reservoir,
    colors_conduit,
    colors_pressure_difference,
]

metrics_heatmap = np.zeros((9, 5))

for metric, axis, color in zip(metric_names, axes, colors):
    axis.axhline(y=1, linestyle="-", color="darkgray", linewidth=1.0)

    sns.stripplot(
        data=metrics,
        x="Caption",
        y=metric,
        hue="Patient",
        palette=color,
        ax=axis,
        jitter=False,
    )
    axis.legend(bbox_to_anchor=(1.0, 1.05))

    axis.grid("on", linewidth=0.5, linestyle="--")

    axis.set_xticks(captions[1:], captions[1:], rotation=45, ha="right")
    axis.set_ylim(0.2, 1.3)
    axis.set_yticks([0.25, 0.5, 0.75, 1.0, 1.25], [0.25, 0.5, 0.75, 1.0, 1.25])

    for i, c in enumerate(captions[1:]):
        data_caption = metrics[metrics["Caption"] == c]
        for j, m in enumerate(metric_names):
            data_metrics = data_caption[m]
            if m == "conduit relative" and c == captions[-1]:
                print(c, np.min(data_metrics), np.max(data_metrics))

            metrics_heatmap[i][j] = data_metrics.mean()

axes[0].set_ylabel("A-loop area,\nrelative to baseline (-)")
axes[1].set_ylabel("Booster function,\nrelative to baseline (-)")
axes[2].set_ylabel("Reservoir function,\nrelative to baseline(-)")
axes[3].set_ylabel("Conduit function,\nrelative to baseline(-)")
axes[4].set_ylabel("Upstroke pressure diff.,\nrelative to baseline(-)")

# for axis in axes:
#    axis.legend(["Patient 1", "Patient 2", "Patient 3"])

plt.tight_layout()
plt.savefig("single_factor_relative_metrics.pdf")

metrics_percentage = []

for i in range(len(metrics_heatmap[0])):
    metrics_percentage.append([])
    for j in range(len(metrics_heatmap)):
        val = metrics_heatmap[j][i]
        percent = 100 * (val - 1)
        if percent > 0:
            str_percent = f"{val:.3f}\n(+{percent:.0f}%)"
        else:
            str_percent = f"{val:.3f}\n({percent:.0f}%)"
        metrics_percentage[i].append(str_percent)


plt.figure(figsize=(8.5, 3.0))
colors = ["#C21E56", "#FFFFFF", "#3A5311"]
cmap = LinearSegmentedColormap.from_list("name", colors)

metric_labels = [
    "A-loop area",
    "Booster function",
    "Reservoir function",
    "Conduit function",
    "Upstroke\npressure difference",
]
sns.heatmap(
    metrics_heatmap.transpose(),
    cmap=cmap,
    annot=metrics_percentage,
    vmin=0.25,
    vmax=1.75,
    cbar_kws={"label": "Average change for all patients;\nrelative to baseline (-)"},
    annot_kws={"color": "black"},
    fmt="s",
)
plt.xticks(
    [i + 0.5 for i in range(len(captions[1:]))], captions[1:], rotation=20, ha="right"
)
plt.yticks([i + 0.5 for i in range(len(metric_names))], metric_labels, rotation=0)

plt.tight_layout()
plt.savefig("heatmap_average_change.pdf", dpi=300)
plt.show()
