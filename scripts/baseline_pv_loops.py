"""

Åshild Telle / University of Washington / 2025

This script plots PV loops at baseline and calculates derived metrics for all 3 patients.
Note that metric annotations were added manually.

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from collections import defaultdict
from metrics import get_metrics, get_loops, plot_loops, get_pv_loop

plt.rcParams.update({"mathtext.default": "regular"})

import yaml

with open('mainfolder.yaml', 'r') as file:
    config = yaml.safe_load(file)
    mainfolder = config['input_data_org_fib']

cases = ["P1", "P2", "P3"]

def get_indices(fin):
    return np.load(fin, allow_pickle=True).item()

fin = "PA_indices_single_factors_original_fibrosis.npy"
PA_indices = get_indices(fin)

fig, axes = plt.subplots(1, 3, figsize=(10, 3.5), sharey=True)

metrics = defaultdict(list)

for axis, cas in zip(axes, cases):

    factor = "baseline"
    subfolder = f"{cas}/{factor}"
    fin = f"{mainfolder}/{subfolder}/cav.LA.csv"

    p_start, a_start = PA_indices[cas][factor]
    _, _, A_loop_volume, A_loop_pressure = get_loops(fin, p_start, a_start)

    axis.fill(A_loop_volume, A_loop_pressure, color="tab:red", alpha=0.4)

    volume, pressure = get_pv_loop(fin)
    A_loop_area, booster, conduit, reservoir, pressure_difference = get_metrics(
        volume, pressure, A_loop_volume, A_loop_pressure
    )

    print(cas)
    print(f"A-loop area: {A_loop_area:.2f}")
    print(f"Booster function: {booster:.2f}")
    print(f"Reservoir function: {reservoir:.2f}")
    print(f"Conduit function: {conduit:.2f}")
    print(f"Upstroke pressure difference: {pressure_difference:2f}")
    print("")

    booster_start = 60
    reservoir_start = np.argmin(volume)
    conduit_start = np.argmax(volume)

    booster_volume = volume[booster_start:reservoir_start]
    booster_pressure = pressure[booster_start:reservoir_start]

    reservoir_volume = volume[reservoir_start:conduit_start]
    reservoir_pressure = pressure[reservoir_start:conduit_start]

    conduit_volume = volume[conduit_start:] + volume[:booster_start]
    conduit_pressure = pressure[conduit_start:] + pressure[:booster_start]

    axis.plot(
        booster_volume, booster_pressure, "-", color="black", label="Booster phase"
    )
    axis.plot(
        reservoir_volume,
        reservoir_pressure,
        "--",
        color="black",
        label="Reservoir phase",
    )
    axis.plot(
        conduit_volume, conduit_pressure, ":", color="black", label="Conduit phase"
    )

    markersize = 10
    P = 7.0
    V = 0.95 * np.min(A_loop_volume)

    if cas == "AF5":
        V = 0.75 * np.min(A_loop_volume)

    vmax = booster_volume[0]
    vmin = booster_volume[-1]
    n = np.argmax(booster_pressure)

    dd = (5, 10)
    
    axis.plot(
        [vmax, vmax],
        [P, booster_pressure[0]],
        "--",
        dashes=dd,
        color="black",
        linewidth=1,
    )

    axis.plot(
        [vmin, vmin],
        [P - 1.5, reservoir_pressure[0]],
        "--",
        dashes=dd,
        color="black",
        linewidth=1,
    )
    axis.plot(
        [np.max(volume), np.max(volume)],
        [P - 1.5, reservoir_pressure[-1]],
        "--",
        dashes=dd,
        color="black",
        linewidth=1,
    )

    axis.plot(
        [V, booster_volume[n]],
        [booster_pressure[n], booster_pressure[n]],
        "--",
        dashes=dd,
        color="black",
        linewidth=1,
    )
    axis.plot(
        [V, vmax],
        [booster_pressure[0], booster_pressure[0]],
        "--",
        dashes=dd,
        color="black",
        linewidth=1,
    )

    # BOOSTER
    N = 100
    axis.plot(
        np.linspace(vmin, vmax, N), P * np.ones(N), "-", linewidth=2, color="tab:orange"
    )
    axis.plot([vmin, vmax], [P, P], "|", color="tab:orange", markersize=10)

    # CONDUIT
    axis.plot(
        np.linspace(vmax, np.max(volume), N),
        (P) * np.ones(N),
        ":",
        linewidth=2,
        color="royalblue",
    )
    axis.plot([vmax, np.max(volume)], [P, P], "|", color="royalblue", markersize=10)

    # RESERVOIR
    axis.plot(
        np.linspace(vmin, np.max(volume), N),
        (P - 1.5) * np.ones(N),
        "--",
        linewidth=2,
        color="tab:green",
    )
    axis.plot(
        [vmin, np.max(volume)],
        [P - 1.5, P - 1.5],
        "|",
        color="tab:green",
        markersize=10,
    )

    # PRESSURE
    axis.plot(
        np.ones(N) * V,
        np.linspace(booster_pressure[0], np.max(booster_pressure), N),
        "-",
        color="indigo",
        linewidth=2,
    )
    axis.plot(
        [V, V],
        [booster_pressure[0], np.max(booster_pressure)],
        marker=(2, 0, -90),
        color="indigo",
        markersize=markersize,
    )

    axis.set_xlabel("Volume (mL)")

axes[0].set_ylabel("Pressure (mmHg)")

axes[0].set_title("Patient 1")
axes[1].set_title("Patient 2")
axes[2].set_title("Patient 3")

axes[0].set_xlim(85, 155)
axes[1].set_xlim(85, 155)
axes[2].set_xlim(25, 95)


axis.legend(bbox_to_anchor=(1.05, 1.0))
plt.tight_layout()
plt.savefig("baseline_metrics.pdf")
plt.show()
