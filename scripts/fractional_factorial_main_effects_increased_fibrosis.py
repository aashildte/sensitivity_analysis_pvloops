"""

Åshild Telle / University of Washington / 2025

This script was used to perform FFD main effect analysis for extended fibrosis burden,
comparing it to original fibrosis burden.

"""

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from statannotations.Annotator import Annotator
from scipy.stats import ttest_ind
from scipy.stats import anderson
from matplotlib.colors import LinearSegmentedColormap

from metrics import get_metrics, get_pv_loop, get_loops


def annotate(p):
    if 5.00e-02 < p <= 1.00e00:
        return "ns"
    elif 1.00e-02 < p <= 5.00e-02:
        return "*"
    elif 1.00e-03 < p <= 1.00e-02:
        return "**"
    elif 1.00e-04 < p <= 1.00e-03:
        return "***"
    elif p <= 1.00e-04:
        return "****"


def analyze_metrics(folder, PA_indices, baseline):
    metrics = defaultdict(list)
    for cas in cases:
        for key, value in runs.items():
            fin = f"{folder}/{cas}_{key}_{value}_{key}_{value}/cav.LA.csv"

            volume, pressure = get_pv_loop(fin)
            p_start, a_start = PA_indices[cas][key]
            _, _, A_loop_volume, A_loop_pressure = get_loops(fin, p_start, a_start)

            A_loop_area, booster, conduit, reservoir, pressure_difference = get_metrics(
                volume, pressure, A_loop_volume, A_loop_pressure
            )

            for i, v in enumerate(value):
                metrics["Caption"] += [captions[i]]
                metrics["A-loop area"].append(
                    A_loop_area / baseline[cas]["A-loop area"]
                )
                metrics["booster"].append(booster / baseline[cas]["booster"])
                metrics["reservoir"].append(reservoir / baseline[cas]["reservoir"])
                metrics["conduit"].append(conduit / baseline[cas]["conduit"])
                metrics["pressure_difference"].append(
                    pressure_difference / baseline[cas]["pressure_difference"]
                )

                metrics["Category"] += [v]

                metrics["Run"] += [key]
                metrics["Case"] += [case_patient(cas)]

    df_metrics = pd.DataFrame(data=metrics)

    return df_metrics


mainfolder_o = "/data2/aashild/sensitivityanalysis/SA_gen2.2/original_fibrosis"
mainfolder_e = "/data2/aashild/sensitivityanalysis/SA_gen2.2/extended_fibrosis"
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


fin = "PA_indices_single_factors_original_fibrosis.npy"
PA_indices_single_factors = get_indices(fin)

fin = "PA_indices_fractional_factorial_original_fibrosis.npy"
PA_indices_fractional_factorial_original = get_indices(fin)

fin = "PA_indices_fractional_factorial_extended_fibrosis.npy"
PA_indices_fractional_factorial_extended = get_indices(fin)

num_factors = 9
num_metrics = 5

cases = ["AF2", "AF4", "AF5"]
baseline = {"AF2": {}, "AF4": {}, "AF5": {}}

for cas in cases:
    fin = f"{mainfolder_o}/{cas}_baseline_baseline/cav.LA.csv"
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
    r"CV$_L$",
    r"CV$_T$",
    r"$I_{K1}$",
    r"$I_{CaL}$",
    r"$I_{Na}$",
    r"$\mu$",
    r"$T_a$",
    r"Stiffness$_L$",
    r"Stiffness$_T$",
]


def case_patient(cas):
    if cas == "AF2":
        return "Patient 1"
    if cas == "AF5":
        return "Patient 2"
    if cas == "AF5":
        return "Patient 3"


palettes = [
    ["salmon", "maroon"],
    ["xkcd:faded orange", "xkcd:rust"],
    ["limegreen", "darkgreen"],
    ["lightblue", "navy"],
    ["violet", "purple"],
]

df_O = analyze_metrics(mainfolder_o, PA_indices_fractional_factorial_original, baseline)
df_E = analyze_metrics(mainfolder_e, PA_indices_fractional_factorial_extended, baseline)


fig, axes = plt.subplots(
    num_metrics, len(captions), sharex=True, sharey=True, figsize=(7, 8)
)
metrics = list(baseline["AF5"].keys())

significance = []

metrics_heatmap = np.zeros((9, 5))
for i in range(len(captions)):
    for m in range(num_metrics):
        print("Considering: ", captions[i], metrics[m])
        data_ca = df_E[df_E["Caption"] == captions[i]]
        data_ca_o = df_O[df_E["Caption"] == captions[i]]

        axes[m][i].axhline(y=1, linestyle="--", color="black", linewidth=0.5)
 
        if i==len(captions)-1:
            sns.pointplot(
                data=data_ca_o,
                x="Category",
                y=metrics[m],
                ax=axes[m][i],
                color="gray",
                zorder=1,
                hue_order=["B", "F"],
                errorbar="sd",
                label="Original\nfibrosis",
            )
            sns.pointplot(
                data=data_ca,
                x="Category",
                y=metrics[m],
                ax=axes[m][i],
                color=palettes[m][-1],
                zorder=2,
                hue_order=["B", "F"],
                errorbar="sd",
                label="Extended\nfibrosis",
            )
        else:
            sns.pointplot(
                data=data_ca_o,
                x="Category",
                y=metrics[m],
                ax=axes[m][i],
                color="gray",
                zorder=1,
                hue_order=["B", "F"],
                errorbar="sd",
            )
            sns.pointplot(
                data=data_ca,
                x="Category",
                y=metrics[m],
                ax=axes[m][i],
                color=palettes[m][-1],
                zorder=2,
                hue_order=["B", "F"],
                errorbar="sd",
            )

        axes[m][i].set(xlabel=captions[i])

        B_data = data_ca[data_ca["Category"] == "B"]
        F_data = data_ca[data_ca["Category"] == "F"]
        B_value = B_data[metrics[m]].mean()
        F_value = F_data[metrics[m]].mean()

        metrics_heatmap[i][m] = 100 * (F_value - B_value) / B_value

        pairs = [("B", "F")]
        annotator = Annotator(
            axes[m][i], pairs, data=data_ca, x="Category", y=metrics[m]
        )
        annotator.configure(test="t-test_ind", text_format="star", loc="outside")
        annotator.apply_and_annotate()

        group1 = data_ca[data_ca["Category"] == "F"][metrics[m]]
        group2 = data_ca_o[data_ca["Category"] == "F"][metrics[m]]

        t_stat, p = ttest_ind(group1, group2)

        if p < 0.01:
            significance.append(
                captions[i] + " " + metrics[m] + " " + str(p) + " " + annotate(p)
            )
        print("Considering: ", captions[i], metrics[m])
        print("\nFF comparison:")
        print("t, p", t_stat, p)

        group1 = data_ca[data_ca["Category"] == "F"][metrics[m]]
        group2 = data_ca[data_ca["Category"] == "F"][metrics[m]]

        print("statistical test: ")
        print(anderson(group1))
        print(anderson(group2))

print(significance)

axes[0][-1].legend(bbox_to_anchor=(1.1, 1.0))
axes[1][-1].legend(bbox_to_anchor=(1.1, 1.0))
axes[2][-1].legend(bbox_to_anchor=(1.1, 1.0))
axes[3][-1].legend(bbox_to_anchor=(1.1, 1.0))
axes[4][-1].legend(bbox_to_anchor=(1.1, 1.0))

axes[0][0].set_ylabel("A-loop area,\nrelative to baseline (-)")
axes[1][0].set_ylabel("Booster function,\nrelative to baseline (-)")
axes[2][0].set_ylabel("Reservoir function,\nrelative to baseline(-)")
axes[3][0].set_ylabel("Conduit function,\nrelative to baseline(-)")
axes[4][0].set_ylabel("Upstroke pressure diff.,\nrelative to baseline(-)")

plt.tight_layout()
plt.savefig("comparison_fibrosis_levels.pdf", dpi=300)
plt.show()

metrics_percentage = []

for j in range(len(metrics_heatmap[0])):
    metrics_percentage.append([])
    for i in range(len(metrics_heatmap)):
        val = metrics_heatmap[i][j]
        if abs(val) < 0.4:
            str_percent = f"{val:.0f}%"
        elif val > 0:
            str_percent = f"+{val:.0f}%"
        else:
            str_percent = f"{val:.0f}%"
        metrics_percentage[j].append(str_percent)

plt.figure(figsize=(7.0, 4.0))
colors = ["#C21E56", "#FFFFFF", "#3A5311"]
cmap = LinearSegmentedColormap.from_list("name", colors)
metric_labels = [
    "A-loop\narea",
    "Booster\nfunction",
    "Reservoir\nfunction",
    "Conduit\nfunction",
    "Upstroke\npressure\ndifference",
]
sns.heatmap(
    metrics_heatmap.transpose(),
    cmap=cmap,
    annot=metrics_percentage,
    vmin=-75,
    vmax=75,
    cbar_kws={"label": "Average change for all patients;\nrelative to baseline (-)"},
    annot_kws={"color": "black"},
    fmt="s",
)


captions = [
    r"Reduced CV$_L$",
    r"Reduced CV$_T$",
    r"$0.5 \times I_{K1}$",
    r"$0.5 \times I_{CaL}$",
    r"$0.6 \times I_{Na}$",
    r"$0.5 \times \mu$",
    r"$0.5 \times T_a$",
    r"$2 \times Load_L$",
    r"$2 \times Load_T$",
]

metric_labels = [
    "A-loop area",
    "Booster function",
    "Reservoir function",
    "Conduit function",
    "Upstroke pressure\ndifference",
]

plt.xticks([i + 0.5 for i in range(len(captions))], captions, rotation=20, ha="right")

plt.tight_layout()
plt.savefig("heatmap_main_effect_fractional_factorial_extended.pdf", dpi=300)
plt.show()
