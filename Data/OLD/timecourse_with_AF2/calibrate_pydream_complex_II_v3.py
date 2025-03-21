"""
Generated by pydream_it
PyDREAM run script for complex_II_v3.py 
"""
from pydream.core import run_dream
from pysb.simulator import ScipyOdeSimulator
import numpy as np
from pydream.parameters import SampledParam
from pydream.convergence import Gelman_Rubin
from scipy.stats import norm, uniform
from complex_II_v3 import model
import os
import re

# DREAM Settings
# Number of chains - should be at least 3.
nchains = 5
# Number of iterations
niterations = 10000

# Initialize PySB solver object for running simulations. Simulation timespan should match experimental data.
files = sorted(os.listdir(''))
exp_time_files = [f for f in files if re.search(r'exp_data_time_\d+', f)]
experiments_time = [np.genfromtxt(file, delimiter=',', names=True) for file in exp_time_files]
n_experiments = len(experiments_time)
tspan = []
tspan_mask = []
for exp_time in experiments_time:
    tspan.append([])
    for name in exp_time.dtype.names:
        tspan[-1] += [t for t in exp_time[name] if not np.isnan(t)]
    tspan[-1] = sorted(list(set(tspan[-1])))  # get a common set of time points for simulations
    tspan_mask.append({})  # for each species, need to mark which time points we have data for
    for name in exp_time.dtype.names:
        tspan_mask[-1][name] = [False] * len(tspan[-1])
        for i in range(len(tspan[-1])):
            if tspan[-1][i] in exp_time[name]:
                tspan_mask[-1][name][i] = True
solver = ScipyOdeSimulator(model)
parameters_idxs = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71]
rates_mask = [False, False, False, False, False, False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False, False]
param_values = np.array([p.value for p in model.parameters])

# USER must add commands to import/load any experimental data for use in the likelihood function!
exp_avg_files = [f for f in files if re.search(r'exp_data_avg_\d+', f)]
experiments_avg = [np.genfromtxt(file, delimiter=',', names=True) for file in exp_avg_files]
exp_se_files = [f for f in files if re.search(r'exp_data_se_\d+', f)]
experiments_se = [np.genfromtxt(file, delimiter=',', names=True) for file in exp_se_files]
like_data = []
for exp_avg, exp_se in zip(experiments_avg, experiments_se):
    like_data.append({})
    for name in exp_avg.dtype.names:
        # remove any nans, which will happen if the time points are different for different species
        avg = [e for e in exp_avg[name] if not np.isnan(e)]
        se = [e for e in exp_se[name] if not np.isnan(e)]
        like_data[-1][name] = norm(loc=avg, scale=se)


# USER must define a likelihood function!
def likelihood(position):
    y = np.copy(position)
    logp_data = [0] * n_experiments
    for n in range(n_experiments):
        param_values[rates_mask] = 10 ** y
        sim = solver.run(tspan=tspan[n], param_values=param_values).all
        for sp in like_data[n].keys():
            logp_data[n] += np.sum(like_data[n][sp].logpdf(sim[sp][tspan_mask[n][sp]]))
        if np.isnan(logp_data[n]):
            logp_data[n] = -np.inf
    return sum(logp_data)


sampled_params_list = list()
sampled_params_list.append(SampledParam(norm, loc=np.log10(1000.0), scale=2.0))  # kf_a_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1000.0), scale=2.0))  # kf_a_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(10.0), scale=2.0))  # kf_a_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1000.0), scale=2.0))  # kf_a_fadnc_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_fadnc_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_fadnc_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(10.0), scale=2.0))  # kr_a_fadnc_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(10.0), scale=2.0))  # kf_a_fadnc_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_fadnc_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1000.0), scale=2.0))  # kf_a_dicarb_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_dicarb_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_dicarb_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_dicarb_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_dicarb_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_dicarb_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_af2_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1000.0), scale=2.0))  # kr_a_af2_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_af2_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_af2_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_af2_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_af2_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_af4_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1000.0), scale=2.0))  # kr_a_af4_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_af4_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_af4_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_af4_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_af4_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_fadnc_dicarb_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_fadnc_dicarb_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_fadnc_dicarb_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_fadnc_dicarb_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_fadnc_af2_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_fadnc_af2_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_fadnc_af2_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_fadnc_af2_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_fadnc_af4_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_fadnc_af4_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_fadnc_af4_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_fadnc_af4_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_dicarb_af2_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_dicarb_af2_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_dicarb_af2_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_dicarb_af2_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_dicarb_af4_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_dicarb_af4_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_dicarb_af4_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_dicarb_af4_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_af2_af4_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_af2_af4_binds_fadnc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_af2_af4_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kr_a_af2_af4_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1000000.0), scale=2.0))  # kcat_a_fadnc_dicarb_to_fadc
sampled_params_list.append(SampledParam(norm, loc=np.log10(60.0), scale=2.0))  # Km_a_fadnc_dicarb_to_fadc
sampled_params_list.append(SampledParam(uniform, loc=np.log10(4.0)-0.25, scale=0.5))  # n_Hill
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # k_a_fadnc_dicarb_af2_to_fadc
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # kf_a_fadc_dicarb_af2_binds_af4
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # k_a_fadnc_dicarb_af4_binds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # k_a_fadnc_af2_af4_binds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # k_a_fadc_dicarb_unbinds_dicarb
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # k_a_fadc_af2_unbinds_af2
sampled_params_list.append(SampledParam(norm, loc=np.log10(1.0), scale=2.0))  # k_a_fadc_af4_unbinds_af4

if __name__ == '__main__':

    sampled_params, log_ps = run_dream(parameters=sampled_params_list,
                                       likelihood=likelihood,
                                       niterations=niterations,
                                       nchains=nchains,
                                       multitry=False,
                                       gamma_levels=4,
                                       adapt_gamma=True,
                                       history_thin=1,
                                       model_name='dreamzs_%dchain' % nchains,
                                       verbose=True)
    total_iterations = niterations
    burnin = int(total_iterations / 2)
    # Save sampling output (sampled parameter values and their corresponding logps).
    for chain in range(len(sampled_params)):
        np.save('dreamzs_%dchain_sampled_params_chain_%d_%d' %
                (nchains, chain, total_iterations), sampled_params[chain])
        np.save('dreamzs_%dchain_logps_chain_%d_%d' % (nchains, chain, total_iterations), log_ps[chain])
    old_samples = sampled_params

    # Check convergence and continue sampling if not converged
    GR = Gelman_Rubin(sampled_params)
    print('At iteration: ', total_iterations, ' GR = ', GR)
    np.savetxt('dreamzs_%dchain_GelmanRubin_iteration_%d.txt' % (nchains, total_iterations), GR)
    if np.any(GR > 1.2):
        starts = [sampled_params[chain][-1, :] for chain in range(nchains)]
        converged = False
        while not converged:
            total_iterations += niterations
            burnin += niterations
            sampled_params, log_ps = run_dream(parameters=sampled_params_list,
                                               likelihood=likelihood,
                                               niterations=niterations,
                                               nchains=nchains,
                                               start=starts,
                                               multitry=False,
                                               gamma_levels=4,
                                               adapt_gamma=True,
                                               history_thin=1,
                                               model_name='dreamzs_%dchain' % nchains,
                                               verbose=True,
                                               restart=True)
            for chain in range(len(sampled_params)):
                np.save('dreamzs_%dchain_sampled_params_chain_%d_%d' %
                        (nchains, chain, total_iterations), sampled_params[chain])
                np.save('dreamzs_%dchain_logps_chain_%d_%d' % (nchains, chain, total_iterations), log_ps[chain])
            old_samples = [np.concatenate((old_samples[chain], sampled_params[chain])) for chain in range(nchains)]
            GR = Gelman_Rubin(old_samples)
            print('At iteration: ', total_iterations, ' GR = ', GR)
            np.savetxt('dreamzs_%dchain_GelmanRubin_iteration_%d.txt' % (nchains, total_iterations), GR)
            if np.all(GR < 1.2):
                converged = True

    # Plot output
    try:
        from pydream_it import plot_param_dist, plot_log_likelihood, plot_time_courses, \
            get_unique_samples_for_simulation

        total_iterations = len(old_samples[0])
        # parameter distributions
        print('Plotting parameter distributions')
        samples = np.concatenate(tuple([old_samples[chain][burnin:, :] for chain in range(nchains)]))
        for n in range(n_experiments):
            samples_n = samples  # TODO: Add code to get samples for the nth experiment
            plot_param_dist(samples_n, [model.parameters[i].name for i in parameters_idxs],
                            suffix='_exp_%d' % n)
        # log likelihood
        print('Plotting log-likelihoods')
        log_ps = []
        n_files = int(total_iterations / niterations)
        for chain in range(nchains):
            log_ps.append(np.concatenate(
                tuple(np.load('dreamzs_%dchain_logps_chain_%d_%d.npy' % (nchains, chain, niterations * (i+1))).flatten()
                      for i in range(n_files))))
        plot_log_likelihood(log_ps, cutoff=2)
        # time courses
        print('Plotting time courses')
        log_ps = np.concatenate(tuple(log_ps[i][burnin:] for i in range(nchains)))
        for n in range(n_experiments):
            print('Experiment %d' % n)
            tspan = np.linspace(tspan[n][0], tspan[n][-1], int((tspan[n][-1] - tspan[n][0]) * 10 + 1))
            samples_n = samples  # TODO: Add code to get samples for the nth experiment
            samples_n, counts = get_unique_samples_for_simulation(samples_n, log_ps, cutoff=2)
            param_values = np.array([param_values] * len(samples_n))
            for i in range(len(param_values)):
                param_values[i][parameters_idxs] = 10 ** samples_n[i]
            print('Running %d simulations' % len(param_values))
            output_all = solver.run(tspan=tspan, param_values=param_values).all
            plot_time_courses(experiments_avg[n].dtype.names, tspan, output_all, counts=counts,
                              exp_data=(experiments_time[n], experiments_avg[n], experiments_se[n]),
                              suffix='_exp_%d' % n)
        print('DONE')

    except ImportError:
        pass

else:
    run_kwargs = {'parameters': sampled_params_list, 'likelihood': likelihood, 'niterations': niterations,
                  'nchains': nchains, 'multitry': False, 'gamma_levels': 4, 'adapt_gamma': True, 'history_thin': 1,
                  'model_name': 'dreamzs_%dchain' % nchains, 'verbose': True}
