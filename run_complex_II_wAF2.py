from complex_II_v3 import model
from pysb.simulator import ScipyOdeSimulator
from param_calibration import *
from SIM_PROTOCOLS.sim_protocols import SequentialInjections

model.parameters['BCD_init'].value = 0

solver = ScipyOdeSimulator(model)
sim_protocol = SequentialInjections(solver, t_equil=100, perturb_day_amount={"Dicarb(a=None)": (0, 10000)})

custom_priors = {'n_Hill': ('uniform', 0.5)}
no_sample = ['A_init', 'FAD_init', 'Dicarb_init', 'AF2_init', 'AF4_init', 'BCD_init', 'kf_a_binds_af4',
             'kr_a_binds_af4', 'kf_a_fadnc_binds_af4', 'kr_a_fadnc_binds_af4', 'kf_a_dicarb_binds_af4',
             'kr_a_dicarb_binds_af4', 'kf_a_af2_binds_af4', 'kr_a_af2_binds_af4', 'kf_a_af4_binds_fadnc',
             'kr_a_af4_binds_fadnc', 'kf_a_af4_binds_dicarb', 'kr_a_af4_binds_dicarb', 'kf_a_af4_binds_af2',
             'kr_a_af4_binds_af2', 'kf_a_fadnc_dicarb_binds_af4', 'kr_a_fadnc_dicarb_binds_af4',
             'kf_a_fadnc_af2_binds_af4', 'kr_a_fadnc_af2_binds_af4', 'kf_a_fadnc_af4_binds_dicarb',
             'kr_a_fadnc_af4_binds_dicarb', 'kf_a_fadnc_af4_binds_af2', 'kr_a_fadnc_af4_binds_af2',
             'kf_a_dicarb_af2_binds_af4', 'kr_a_dicarb_af2_binds_af4', 'kf_a_dicarb_af4_binds_fadnc',
             'kr_a_dicarb_af4_binds_fadnc', 'kf_a_dicarb_af4_binds_af2', 'kr_a_dicarb_af4_binds_af2',
             'kf_a_af2_af4_binds_fadnc', 'kr_a_af2_af4_binds_fadnc', 'kf_a_af2_af4_binds_dicarb',
             'kr_a_af2_af4_binds_dicarb', 'kf_a_fadc_dicarb_af2_binds_af4', 'k_a_fadnc_dicarb_af4_binds_af2',
             'k_a_fadnc_af2_af4_binds_dicarb', 'k_a_fadc_af4_unbinds_af4', 'k_a_fadnc_binds_bcd', 'k_a_fadc_binds_bcd']

# for x in [p.name for p in model.parameters if p.name not in no_sample]:
#     print(x)
# quit()

exp_data_file = 'Data/Complex_II_Experiment_Data_WithAF2.csv'

if __name__ == '__main__':
    # from pysb import ANY
    # from pysb.util import alias_model_components
    #
    # # run a pre-simulation
    # alias_model_components()
    # Obs_to_plot = {
    #     'Free_SDHA': A(af2=None, af4=None, fad=None, dicarb=None),
    #     'SDHA_SDHAF2': A(af2=ANY, af4=None, fad=None, dicarb=None),
    #     'SDHA_FAD': A(af2=None, af4=None, fad=ANY, dicarb=None),
    #     'SDHA_SDHAF2_FAD': A(af2=ANY, af4=None, fad=ANY, dicarb=None),
    #     'Free_SDHAF2': AF2(a=None),
    #     'Free_FAD': FAD(a=None)
    # }
    #
    # tspan = np.linspace(0, 1, 101)
    # output = solver.run(tspan)
    # for obs in Obs_to_plot.keys():
    #     plt.plot(tspan, output.observable(Obs_to_plot[obs]), lw=2, label=obs)
    # plt.xlabel("time (min)")
    # plt.ylabel("concentration (uM)")
    # plt.yscale('log')
    # plt.legend(loc='lower left', ncol=2, bbox_to_anchor=(0.1, 0.1))
    # plt.tight_layout()
    #
    # # run a simulation using the SimulationProtocol object
    # plt.figure()
    # data = np.genfromtxt(exp_data_file, dtype=None, delimiter=',',  names=True, encoding="utf_8_sig")
    # tspan = np.linspace(0, 23, 231)
    # output = sim_protocol.run(tspan, [p.value for p in model.parameters])
    # plt.plot(tspan, output['pct_flavinylation'], lw=2, label='simulation')
    # plt.errorbar([d['time'] for d in data], [d['average'] for d in data], yerr=[d['stderr'] for d in data], fmt='ro',
    #              capsize=6, label='experiment')
    # plt.xlabel('time (min)')
    # plt.ylabel('pct_flavinylation')
    # plt.legend(loc=0)
    # plt.tight_layout()
    #
    # plt.show()

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      sim_protocol,
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5, plot_results=True)
