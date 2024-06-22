import numpy as np
import matplotlib.pyplot as plt
import os
import re

paths = ['timecourse_no_AF2', 'timecourse_with_AF2']
labels = ['- AF2', '+ AF2']

for path, label in zip(paths, labels):
    files = sorted([file for file in os.listdir(path) if os.path.isfile(os.path.join(path, file))])
    file_avg = None
    file_se = None
    file_time = None
    for file in files:
        if re.search(r'_avg_.*\.csv$', file):
            file_avg = file
        elif re.search(r'_se_.*\.csv$', file):
            file_se = file
        elif re.search(r'_time_.*\.csv$', file):
            file_time = file

    # file_avg = 'complex_II_v3_exp_data_avg_0.csv'
    avg = np.genfromtxt(os.path.join(path, file_avg), dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
    #
    # file_se = 'complex_II_v3_exp_data_se_0.csv'
    se = np.genfromtxt(os.path.join(path, file_se), dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
    #
    # file_time = 'complex_II_v3_exp_data_time_0.csv'
    time = np.genfromtxt(os.path.join(path, file_time), dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
    #
    plt.errorbar(time['pct_flavinylation'], avg['pct_flavinylation'], se['pct_flavinylation'], ls='None', marker='o', ms=8,
                 capsize=10, label=label)
    plt.xlabel('Time (min)')
    plt.ylabel('% flavinylation')
    plt.legend(loc=0)
#
plt.show()

