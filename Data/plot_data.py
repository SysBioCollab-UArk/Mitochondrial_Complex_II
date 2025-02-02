import numpy as np
import matplotlib.pyplot as plt
import os
import re

files = ["Flav_Time.csv", "Flav_Fumarate.csv", "Flav_FAD.csv"]
for file in files:
    print(file)

    data = np.genfromtxt(file, dtype=None, delimiter=",", names=True, encoding="utf_8_sig")
    print(data.dtype.names)

    expt_ids = np.unique([d["expt_id"] for d in data])

    plt.figure(constrained_layout=True)
    for expt_id in expt_ids:
        xvals = [d["xval"] for d in data if d["expt_id"] == expt_id]
        yvals = [d["yval"] for d in data if d["expt_id"] == expt_id]
        stderr = [d["stderr"] for d in data if d["expt_id"] == expt_id]
        plt.errorbar(xvals, yvals, stderr, label=expt_id)

    plt.legend(loc="best")

plt.show()
