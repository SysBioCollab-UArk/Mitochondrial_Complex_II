import numpy as np
import matplotlib.pyplot as plt
import os
import importlib

# path where are the files are (run file, data file, and PyDREAM output files)
path = os.getcwd()

# import everything from run_complex_II_pydream.py in the path
run_pydream_file = os.path.join(path,'run_complex_II_pydream.py')
import_string = run_pydream_file.replace('/', '.').replace('\\', '.').rstrip('.py')
module = importlib.import_module(import_string)  # Import the module
# Update the current namespace with everything from the module
globals().update({k: v for k, v in vars(module).items() if not k.startswith('_')})

exp_data_file = os.path.join(path, os.path.basename(globals()['exp_data_file']))
expt_data = np.genfromtxt(exp_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
expt_ids = list(dict.fromkeys(expt_data['expt_id']))  # get unique expt_ids and retain order of appearance
print(expt_data.dtype.names)
print(expt_ids)

sim_data_file = os.path.join(path, 'SIM_DATA.csv')
sim_data = np.genfromtxt(sim_data_file, dtype=None, delimiter=',', names=True, encoding="utf_8_sig")
print(sim_data.dtype.names)

#####
conc_perturb = []  # TODO: Make a list of dictionaries like this in the run_*_pydream file for dose-response expts
for conc_dict, perturb in zip(['conc_dicarb', 'conc_fad'], ['fum', 'fad']):
    if conc_dict in globals():
        conc_perturb.append({})
        conc_perturb[-1]['%s_wAF2' % perturb] = list(globals()[conc_dict]['wAF2'])
        conc_perturb[-1]['%s_noAF2' % perturb] = list(globals()[conc_dict]['noAF2'])
#####

xlabels = dict(zip(['fum', 'fad'], ['Fumarate', 'FAD']))
start = 0
for cp in conc_perturb:
    plt.figure(constrained_layout=True)
    xlabel = None
    for expt_id in expt_ids:
        print(expt_id)
        if expt_id in cp.keys():
            conc = cp[expt_id]
            # sim data
            yval_min = sim_data['yval_min'][start:start + len(conc)]
            yval_max = sim_data['yval_max'][start:start + len(conc)]
            p = plt.plot(conc, (yval_min + yval_max)/2, ls='--', label='%s (sim)' % expt_id)
            plt.fill_between(conc, yval_min, yval_max, alpha=0.25, color=p[0].get_color(), label='x')
            # expt data
            avg = expt_data['average'][start:start + len(conc)]
            stderr = expt_data['stderr'][start:start + len(conc)]
            plt.errorbar(conc, avg, yerr=stderr, fmt='o', ms=8, capsize=6, color=p[0].get_color(),
                         label='%s (expt)' % expt_id)
            #
            start += len(conc)

            # get the xlabel (Fumarate or FAD)
            if xlabel is None:
                for key in xlabels.keys():
                    if key in expt_id:
                        xlabel = xlabels[key]
                        break

    plt.xlabel('%s (μM)' % xlabel)
    plt.ylabel('% flavinylation')
    plt.legend(loc='best')

    # merge line and fill_between legend handles
    handles, labels = plt.gca().get_legend_handles_labels()
    n_sims = len([label for label in labels if 'sim' in label])
    new_handles = [(handles[n], handles[n + 1]) for n in range(0, n_sims * 2, n_sims)] + list(handles[n_sims * 2:])
    new_labels = [labels[n] for n in range(0, n_sims*2, n_sims)] + list(labels[n_sims*2:])
    plt.legend(new_handles, new_labels, loc='best')

plt.savefig('fig_PyDREAM_dose_response.png')

plt.show()
