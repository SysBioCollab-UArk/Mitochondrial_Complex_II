import matplotlib.pyplot as plt
from util import plot_drc_from_simdata, plot_tc_from_simdata
import numpy as np

basepath = 'SAVED'
directories = ['Flav_Fumarate_FAD_Time']  #['Flav_Fumarate', 'Flav_FAD', 'Flav_Fumarate_FAD']
run_pydream_filename = 'run_complex_II_pydream.py'

label_dict = {'pct_flavinylation': 'Flavinylation', 'fum_wAF2': 'With AF2', 'fum_noAF2': 'W/o AF2',
              'fad_wAF2': 'With AF2', 'fad_noAF2': 'W/o AF2', 'time_wAF2': 'With AF2', 'time_noAF2': 'W/o AF2',
              'percent': '%'}

# plot dose-response curves
expt_doses = [
    # # experiment 0
    # {'fum_wAF2': np.array([0, 0.00625, 0.0125, 0.025, 0.05, 0.1, 0.25, 1]) * 1000,
    #  'fum_noAF2': np.array([0, 0.05, 0.1, 0.25, 0.5, 1]) * 1000,
    #  'xlabel': r'Fumarate ($\mu$M)'},
    # # experiment 1
    # {'fad_wAF2': [0, 1, 2, 4, 6, 20, 40, 100],
    #  'fad_noAF2': [0, 5, 10, 20, 40, 60, 100],
    #  'xlabel': r'FAD ($\mu$M)'},
    # experiment 2
    [{'fum_wAF2': np.array([0, 0.00625, 0.0125, 0.025, 0.05, 0.1, 0.25, 1]) * 1000,
     'fum_noAF2': np.array([0, 0.05, 0.1, 0.25, 0.5, 1]) * 1000,
      'xlabel': r'Fumarate ($\mu$M)'},
     {'fad_wAF2': [0, 1, 2, 4, 6, 20, 40, 100],
     'fad_noAF2': [0, 5, 10, 20, 40, 60, 100],
     'xlabel': r'FAD ($\mu$M)'}]
]

plot_drc_from_simdata(basepath, directories, run_pydream_filename, expt_doses, label_dict)

# plot timecourses
tc_ids = ['time_wAF2', 'time_noAF2']
plot_tc_from_simdata(basepath, directories, run_pydream_filename, tc_ids, label_dict)
