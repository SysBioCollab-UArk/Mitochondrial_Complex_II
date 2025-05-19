from util import plot_drc
import numpy as np

basepath = 'SAVED'
directories = ['Flav_Fumarate', 'Flav_FAD']
run_pydream_filename = 'run_complex_II_pydream.py'
expts_doses = [
    {'fum_wAF2': np.array([0, 0.00625, 0.0125, 0.025, 0.05, 0.1, 0.25, 1]) * 1000,
     'fum_noAF2': np.array([0, 0.05, 0.1, 0.25, 0.5, 1]) * 1000,
     'xlabel': r'Fumarate ($\mu$M)'},
    {'fad_wAF2': [0, 1, 2, 4, 6, 20, 40, 100],
     'fad_noAF2': [0, 5, 10, 20, 40, 60, 100],
     'xlabel': r'FAD ($\mu$M)'},
]
label_dict = {'pct_flavinylation': 'Flavinylation', 'fum_wAF2': 'With AF2', 'fum_noAF2': 'W/o AF2'}

plot_drc(basepath, directories, run_pydream_filename, expts_doses, label_dict)
