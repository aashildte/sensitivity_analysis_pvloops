"""

Åshild Telle / University of Washington / 2025

This script simply extracts volume and pressure for Patient 1; part of the pipeline presentation figure.

"""

import numpy as np
import matplotlib.pyplot as plt

from metrics import get_pv_loop

import yaml

with open('mainfolder.yaml', 'r') as file:
    config = yaml.safe_load(file)
    mainfolder = config['input_data_org_fib']

cas = "P1"
factor = "baseline"
subfolder = f"{cas}/{factor}"
fin = f"{mainfolder}/{subfolder}/cav.LA.csv"
volume, pressure = get_pv_loop(fin)
time = np.linspace(0, 1000, len(volume))

fig, axes = plt.subplots(1, 3, figsize=(9, 2.0))

axes[0].plot(time, volume, color="black")
axes[1].plot(time, pressure, color="black")
axes[2].plot(volume, pressure, color="black")

axes[0].set_xlabel("Time (ms)")
axes[1].set_xlabel("Time (ms)")
axes[2].set_xlabel("Volume (mL)")

axes[0].set_ylabel("Volume (mL)")
axes[1].set_ylabel("Pressure (mmHg)")
axes[2].set_ylabel("Pressure (mmHg)")

plt.tight_layout()
plt.savefig("volume_pressure_pv_loop.pdf", dpi=300)
plt.show()
