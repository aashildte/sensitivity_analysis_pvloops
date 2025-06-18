"""

Åshild Telle / University of Washington / 2025

This script plots volume, pressure, and pv loops for all 10 cycles (convergence plots).

"""

import numpy as np
import matplotlib.pyplot as plt

import yaml

with open('mainfolder.yaml', 'r') as file:
    config = yaml.safe_load(file)
    mainfolder = config['input_data_org_fib']

cases = ["P1", "P2", "P3"]
single_factors = ["baseline"]


def read_all_data(fin):
    pressure = []
    volume = []
    with open(fin, "r") as file:
        i = 0
        for line in file:
            if i < 2:
                i += 1  # skip header
                continue

            data = line.split(",")

            pressure.append(float(data[1]))
            volume.append(float(data[2]))

    return np.array(volume), np.array(pressure)


def cas2patient(cas):
    if cas == "P1":
        return "Patient 1"
    elif cas == "P2":
        return "Patient 2"
    return "Patient 3"


def get_indices(fin):
    return np.load(fin, allow_pickle=True).item()


fin = "PA_indices_single_factors_original_fibrosis.npy"
PA_indices = get_indices(fin)


fig, axes = plt.subplots(3, 3, figsize=(8.0, 6.0), sharey="row")

for j, cas in enumerate(cases):
    print(cas)
    factor = "baseline"
    subfolder = f"{cas}/{factor}"
    fin = f"{mainfolder}/{subfolder}/cav.LA.csv"

    volume, pressure = read_all_data(fin)

    # last_loop_volume = volume[-2001:]
    # last_loop_pressure = pressure[-2001:]

    volume_change = []
    pressure_change = []

    time = np.linspace(0, 1000, 2001)  # len(volume[50:]))

    for i in range(10):
        index1 = 50 + i * 2000
        index2 = index1 + 2000 + 1

        curr_volume = volume[index1:index2]
        curr_pressure = pressure[index1:index2]

        axes[0][j].plot(
            time, curr_pressure, color="red", alpha=i * 0.08 + 0.1, label=str(i + 1)
        )  # f"Cycle number {i+1}")
        axes[1][j].plot(
            time, curr_volume, color="black", alpha=i * 0.08 + 0.1, label=str(i + 1)
        )  # f"Cycle number {i+1}")

        if i > 0:
            volume_diff = (
                100 * np.max(curr_volume - last_loop_volume) / np.max(curr_volume)
            )
            pressure_diff = (
                100 * np.max(curr_pressure - last_loop_pressure) / np.max(curr_pressure)
            )

            if i == 9:
                print(i, "volume diff:", volume_diff)
                print(i, "pressure diff:", pressure_diff)

            volume_change.append(volume_diff)
            pressure_change.append(pressure_diff)

        last_loop_volume = curr_volume[:]
        last_loop_pressure = curr_pressure[:]

    axes[2][j].plot(
        [k + 2 for k in range(9)], pressure_change, "o-", markersize=4, color="red"
    )
    axes[2][j].plot(
        [k + 2 for k in range(9)], volume_change, "o-", markersize=4, color="black"
    )

    axes[0][j].set_xlabel("Time (ms)")
    axes[1][j].set_xlabel("Time (ms)")
    axes[2][j].set_xlabel("Cycle number")

    axes[0][j].set_title(cas2patient(cas))

    axes[2][j].set_xticks([2, 6, 10], [2, 6, 10])

plt.figtext(0.83, 0.95, "Cycle number")
plt.figtext(0.83, 0.64, "Cycle number")

axes[0][0].set_yticks([0, 10, 20], [0, 10, 20])
axes[2][0].set_yticks([0, 25, 50, 75, 100], [0, 25, 50, 75, 100])

axes[0][0].set_ylabel("Pressure (mmHg)")
axes[1][0].set_ylabel("Volume (mL)")
axes[2][0].set_ylabel("Maximum difference\nrelative to\nprevious cycle (%)")

axes[2][0].set_yscale("log")

axes[0][2].legend(bbox_to_anchor=(1.05, 1.05), ncol=2)
axes[1][2].legend(bbox_to_anchor=(1.05, 1.05), ncol=2)
axes[2][2].legend(["Pressure", "Volume"], bbox_to_anchor=(1.85, 1.05))

plt.tight_layout()

plt.savefig("convergence_plot.pdf")
plt.show()
