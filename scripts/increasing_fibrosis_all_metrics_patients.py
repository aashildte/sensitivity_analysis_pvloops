"""

Åshild Telle / University of Washington / 2025

This script plots absolute metric values, for all patients, for increasing % of fibrosis.

"""

import numpy as np
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict

from metrics import get_metrics, get_pv_loop, get_loops


def analyze_metrics(folder, PA_indices_FF, PA_indices_single_factors, baseline, fib):
    metrics = defaultdict(list)
 
    runs = ["baseline", "Run32_FFFFFFFFF"]

    for cas, f_value in zip(cases, fib):
        for run in runs:
            fin = f"{folder}/{cas}_{run}_{run}/cav.LA.csv"
            key = run
            if "Run" in run:
                key = run.split("_")[0]
                p_start, a_start = PA_indices_FF[cas][key]
            else:
                p_start, a_start = PA_indices_single_factors[cas][run]
            volume, pressure = get_pv_loop(fin)
            _, _, A_loop_volume, A_loop_pressure = get_loops(fin, p_start, a_start)

            A_loop_area, booster, conduit, reservoir, pressure_difference = get_metrics(
                volume, pressure, A_loop_volume, A_loop_pressure
            )

            metrics["A-loop area"].append(A_loop_area)  # /baseline[cas]["A-loop area"])
            metrics["booster"].append(booster)  # /baseline[cas]["booster"])
            metrics["reservoir"].append(reservoir)  # /baseline[cas]["reservoir"])
            metrics["conduit"].append(conduit)  # /baseline[cas]["conduit"])
            metrics["pressure_difference"].append(
                pressure_difference
            )  # /baseline[cas]["pressure_difference"])

            metrics["Run"].append(key)
            metrics["Case"] += [cas]
            metrics["Fib"].append(f_value)

    df_metrics = pd.DataFrame(data=metrics)

    return df_metrics


mainfolder_o = "/data2/aashild/sensitivityanalysis/SA_gen2.2/original_fibrosis"
mainfolder_e = "/data2/aashild/sensitivityanalysis/SA_gen2.2/extended_fibrosis"

def get_indices(fin):
    return np.load(fin, allow_pickle=True).item()

fin = "PA_indices_single_factors_original_fibrosis.npy"
PA_indices_single_factors = get_indices(fin)

fin = "PA_indices_single_factors_extended_fibrosis.npy"
PA_indices_single_factors_extended = get_indices(fin)

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


colors = {}

colors["A-loop area"] = ["salmon", "tab:red", "maroon"]
colors["booster"] = ["navajowhite", "tab:orange", "xkcd:rust"]
colors["reservoir"] = ["palegreen", "tab:green", "darkolivegreen"]
colors["conduit"] = ["lightblue", "royalblue", "navy"]
colors["pressure_difference"] = ["violet", "darkorchid", "indigo"]

fib_O = [15.6, 23.9, 17.9]
fib_E = [15.6 * 1.5, 23.9 * 1.5, 17.9 * 1.5]

df_O = analyze_metrics(
    mainfolder_o,
    PA_indices_fractional_factorial_original,
    PA_indices_single_factors,
    baseline,
    fib_O,
)
df_E = analyze_metrics(
    mainfolder_e,
    PA_indices_fractional_factorial_extended,
    PA_indices_single_factors_extended,
    baseline,
    fib_E,
)

metric_names = ["A-loop area", "booster", "reservoir", "conduit", "pressure_difference"]

fig, axes = plt.subplots(1, 5, sharex=True, figsize=(13, 2.5))

for metric, axis in zip(metric_names, axes):
    for i, cas in enumerate(["AF2", "AF4", "AF5"]):
        data_org_case = df_O[df_O["Case"] == cas]
        data_org_case = data_org_case[data_org_case["Run"] == "Run32"]
        data_ext_case = df_E[df_E["Case"] == cas]
        data_ext_case = data_ext_case[data_ext_case["Run"] == "Run32"]

        fib_values = [0, data_org_case["Fib"].item(), data_ext_case["Fib"].item()]
        y_values = [
            baseline[cas][metric],
            data_org_case[metric].item(),
            data_ext_case[metric].item(),
        ]

        axis.plot(fib_values, y_values, "o-", color=colors[metric][i])
        axis.legend([r"$P_1$", r"$P_2$", r"$P_3$"], loc=1)

    axis.set_xlabel("Fibrosis burden (%)")

axes[0].set_ylabel("A-loop area (mJ)")
axes[1].set_ylabel("Booster function (-)")
axes[2].set_ylabel("Reservoir function (-)")
axes[3].set_ylabel("Conduit function (-)")
axes[4].set_ylabel("Upstroke pressure difference (-)")

plt.tight_layout()
plt.savefig("scatter_plots_fully_fibrotic.pdf", dpi=300)
plt.show()
