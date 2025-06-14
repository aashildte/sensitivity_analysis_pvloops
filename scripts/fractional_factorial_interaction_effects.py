"""

Åshild Telle / University of Washington / 2025

This script was used to perform FFD interactive effect analysis.

"""

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns

import statsmodels.api as sm
from collections import defaultdict

from metrics import get_metrics, get_loops, plot_loops, get_pv_loop, analyze_metrics

mainfolder = "/data2/aashild/sensitivityanalysis/SA_gen2.2/original_fibrosis"

def get_PA_indices(fin):
    return np.load(fin, allow_pickle=True).item()


fin = "PA_indices_fractional_factorial_original_fibrosis.npy"
PA_indices_fractional_factorial = get_PA_indices(fin)

fin = "PA_indices_single_factors_original_fibrosis.npy"
PA_indices_single_factors = get_PA_indices(fin)


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

num_factors = 9

baseline = {"AF2": {}, "AF4": {}, "AF5": {}}

for cas in cases:
    fin = f"{mainfolder}/{cas}_baseline_baseline/cav.LA.csv"
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


df_O = analyze_metrics(
    mainfolder, cases, runs, captions, PA_indices_fractional_factorial, baseline
)

metrics = baseline["AF2"].keys()

interactions_FF = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
interactions_FB = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
interactions_BF = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))
interactions_BB = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))


def extract_mean(df, i, j, fib1, fib2, metric):
    df1 = df[df[captions[i]] == fib1]
    if i == j:
        return df1[metric].mean()

    df2 = df1[df1[captions[j]] == fib2]
    return df2[metric].mean()


for metric in metrics:
    for i in range(num_factors):
        for j in range(num_factors):
            FF_mean = extract_mean(df_O, i, j, "F", "F", metric)
            interactions_FF[metric][i][j].append(FF_mean)

            FB_mean = extract_mean(df_O, i, j, "F", "B", metric)
            interactions_FB[metric][i][j].append(FB_mean)

            BF_mean = extract_mean(df_O, i, j, "B", "F", metric)
            interactions_BF[metric][i][j].append(BF_mean)

            BB_mean = extract_mean(df_O, i, j, "B", "B", metric)
            interactions_BB[metric][i][j].append(BB_mean)


def plot_interactions(
    interactions_FF,
    interactions_FB,
    interactions_BF,
    interactions_BB,
    colors,
    axis_hm,
    axes,
):
    interaction_cross = np.full(
        shape=(num_factors, num_factors), fill_value=np.nan, dtype=float
    )

    max_synergy = (0, -1, -1)
    min_synergy = (0, -1, -1)

    for i in range(num_factors):
        for j in range(num_factors):
            avg_B = [interactions_BB[i][j][0], interactions_BF[i][j][0]]
            avg_F = [interactions_FB[i][j][0], interactions_FF[i][j][0]]

            if i == j:
                X = 0
            else:
                X =  np.arcsin(np.cross(np.array(avg_F), np.array(avg_B)) / (np.linalg.norm(avg_F)*np.linalg.norm(avg_B)) )
            interaction_cross[i][j] = X

            if X > max_synergy[0]:
                max_synergy = (X, i, j)
            if X < min_synergy[0]:
                min_synergy = (X, i, j)

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "mycmap", ["black", "white", colors[1]]
    )

    vminmax = 0.1  # max(max_synergy[0], -min_synergy[0])

    sns.heatmap(
        interaction_cross,
        cmap=cmap,
        ax=axis_hm,
        square=True,
        vmax=vminmax,
        vmin=-vminmax,
    )
    axis_hm.set_xticks([x + 0.5 for x in range(num_factors)], captions, rotation=45)
    axis_hm.set_yticks([x + 0.5 for x in range(num_factors)], captions, rotation=0)
    axis_hm.set_title(metric)

    X, i, j = max_synergy

    axes[0].set_title(f"{captions[i]}:{captions[j]}")

    sns.stripplot(
        data=df_O,
        x=captions[j],
        y=metric,
        hue=captions[i],
        palette=colors,
        alpha=0.3,
        order=["B", "F"],
        ax=axes[0],
        legend=None,
    )
    sns.pointplot(
        data=df_O,
        x=captions[j],
        y=metric,
        hue=captions[i],
        palette=colors,
        order=["B", "F"],
        hue_order=["B", "F"],
        ax=axes[0],
        errorbar=None,
        legend=None,
    )

    # axes[0].legend(bbox_to_anchor=(1.1, 1.0))

    X, i, j = min_synergy
    print(X, captions[i], captions[j])
    axes[1].set_title(f"{captions[i]}:{captions[j]}")

    sns.stripplot(
        data=df_O,
        x=captions[j],
        y=metric,
        hue=captions[i],
        palette=colors,
        alpha=0.3,
        order=["B", "F"],
        hue_order=["B", "F"],
        ax=axes[1],
        legend=None,
    )
    sns.pointplot(
        data=df_O,
        x=captions[j],
        y=metric,
        hue=captions[i],
        palette=colors,
        order=["B", "F"],
        hue_order=["B", "F"],
        ax=axes[1],
        errorbar=None,
        legend=None,
    )

color_pairs = [
    ["salmon", "maroon"],
    ["xkcd:faded orange", "xkcd:rust"],
    ["palegreen", "darkgreen"],
    ["lightblue", "navy"],
    ["violet", "purple"],
]


fig_hm, axes_hm = plt.subplots(len(metrics), figsize=(6, 12), sharex=True, sharey=True)

fig, axes = plt.subplots(len(metrics), 2, figsize=(3, 12), sharey=True)


for axis, axis_hm, metric, colors in zip(axes, axes_hm, metrics, color_pairs):
    plot_interactions(
        interactions_FF[metric],
        interactions_FB[metric],
        interactions_BF[metric],
        interactions_BB[metric],
        colors,
        axis_hm,
        axis,
    )

fig_hm.tight_layout()
fig_hm.savefig("heatmaps_interaction.pdf")
plt.show()
fig.tight_layout()
fig.savefig("interaction_curves2.pdf")
plt.show()
