"""

Åshild Telle / University of Washington / 2025

This script contains keys functions for extracting PV loops, calciulating extracted metrics, etc.

"""


import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from collections import defaultdict


# source: https://stackoverflow.com/questions/65532031/how-to-find-number-of-self-intersection-points-on-2d-plot
def intersection(x1, x2, x3, x4, y1, y2, y3, y4):
    d = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    if d:
        xs = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / d
        ys = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / d
        if (
            xs >= min(x1, x2)
            and xs <= max(x1, x2)
            and xs >= min(x3, x4)
            and xs <= max(x3, x4)
        ):
            return xs, ys


def find_intersection_points(x, y):
    xs, ys = [], []
    pos_xs, pos_ys = [], []
    for i in range(len(x) - 1):
        for j in range(i - 1):
            if xs_ys := intersection(
                x[i], x[i + 1], x[j], x[j + 1], y[i], y[i + 1], y[j], y[j + 1]
            ):
                xs.append(xs_ys[0])
                ys.append(xs_ys[1])
                pos_xs.append(i)
                pos_ys.append(j)
    return xs, ys, pos_xs, pos_ys


# stolen from https://stackoverflow.com/questions/66784360/how-to-calculate-the-inside-area-of-a-closed-shape-given-a-set-of-coordinates-de
def PolyArea(x, y):
    return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


def toJoule(x, y):
    # 1 mm Hg mL = .1333 mJ
    return 0.1333 * PolyArea(x, y)


def get_pv_loop(fin):
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

    # plt.figure()
    # plt.plot(volume, pressure)

    # trunkate to the last loop
    N1 = -2001
    pressure = pressure[N1:]
    volume = volume[N1:]

    return volume, pressure


def get_loops(fin, P_loop_start, A_loop_start, plot=False):
    volume, pressure = get_pv_loop(fin)

    # + 2 for buffer so we get closed curves for both

    A_loop_pressure = pressure[P_loop_start : A_loop_start + 2]
    P_loop_pressure = pressure[A_loop_start:] + pressure[: P_loop_start + 2]

    A_loop_volume = volume[P_loop_start : A_loop_start + 2]
    P_loop_volume = volume[A_loop_start:] + volume[: P_loop_start + 2]

    if plot:
        plt.figure()
        plt.plot(P_loop_volume, P_loop_pressure)
        plt.plot(A_loop_volume, A_loop_pressure)
        plt.plot([P_loop_volume[:2]], [P_loop_pressure[:2]], "o")
        plt.plot([A_loop_volume[:2]], [A_loop_pressure[:2]], "*")
        plt.show()

    return P_loop_volume, P_loop_pressure, A_loop_volume, A_loop_pressure


def get_PA_indices(fin, intersection_index, plot=False):
    pressure, volume = get_pv_loop(fin)
    if plot:
        plt.figure()
        plt.plot(volume, pressure)
        # plt.plot(ys, xs, "o")
        # plt.plot([P_loop_volume[:2]], [P_loop_pressure[:2]], "o")
        # plt.plot([A_loop_volume[:2]], [A_loop_pressure[:2]], "*")
        plt.show()

    # sometimes this works better:
    # xs, ys, pos_xs, pos_ys = find_intersection_points(volume, pressure)
    ys, xs, pos_ys, pos_xs = find_intersection_points(pressure, volume)
    print("num intersections: ", len(xs))
    inters = [xs[intersection_index], ys[intersection_index]]

    P_loop_start = pos_xs[intersection_index]
    A_loop_start = pos_ys[intersection_index]

    return P_loop_start, A_loop_start


def plot_loops(P_loop_volume, P_loop_pressure, A_loop_volume, A_loop_pressure, axis):
    axis.plot(P_loop_volume, P_loop_pressure, color="black")
    axis.fill(P_loop_volume, P_loop_pressure, color="gray", alpha=0.4)
    
    axis.plot(A_loop_volume, A_loop_pressure, color="tab:red")
    axis.fill(A_loop_volume, A_loop_pressure, color="tab:red", alpha=0.4)


def get_metrics(total_volume, total_pressure, A_loop_volume, A_loop_pressure):
    AJ = toJoule(A_loop_volume, A_loop_pressure)

    vmax = np.max(total_volume)
    vmin = np.min(total_volume)
    pre_a = total_volume[60]

    booster = (pre_a - vmin) / pre_a
    conduit = (vmax - pre_a) / vmax
    reservoir = (vmax - vmin) / vmin

    P_active = np.max(A_loop_pressure)
    diff_P_active = (P_active - total_pressure[60]) / total_pressure[60]

    return AJ, booster, conduit, reservoir, diff_P_active


def case_patient(cas):
    if cas == "AF2":
        return "Patient 1"
    if cas == "AF4":
        return "Patient 2"
    if cas == "AF5":
        return "Patient 3"


def analyze_metrics(folder, cases, runs, captions, PA_indices, baseline):
    metrics = defaultdict(list)
    for cas in cases:
        for key, value in runs.items():
            fin = f"{folder}/{cas}/{key}_{value}/cav.LA.csv"

            volume, pressure = get_pv_loop(fin)
            p_start, a_start = PA_indices[cas][key]
            _, _, A_loop_volume, A_loop_pressure = get_loops(fin, p_start, a_start)

            A_loop_area, booster, conduit, reservoir, pressure_difference = get_metrics(
                volume, pressure, A_loop_volume, A_loop_pressure
            )

            metrics["A-loop area"].append(A_loop_area / baseline[cas]["A-loop area"])
            metrics["booster"].append(booster / baseline[cas]["booster"])
            metrics["reservoir"].append(reservoir / baseline[cas]["reservoir"])
            metrics["conduit"].append(conduit / baseline[cas]["conduit"])
            metrics["pressure_difference"].append(
                pressure_difference / baseline[cas]["pressure_difference"]
            )

            for i, v in enumerate(value):
                metrics[captions[i]] += [v]

            metrics["Run"] += [key]
            metrics["Case"] += [case_patient(cas)]
    df_metrics = pd.DataFrame(data=metrics)

    return df_metrics
