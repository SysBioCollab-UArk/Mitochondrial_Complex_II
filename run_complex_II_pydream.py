from complex_II_v3 import model
from pysb.simulator import ScipyOdeSimulator
from param_calibration import *
from SIM_PROTOCOLS.sim_protocols import *

model.parameters['AF2_init'].value = 0
model.parameters['BCD_init'].value = 0

solver = ScipyOdeSimulator(model)

conc_dicarb = {
    'wAF2' : [0, 0.00625, 0.0125, 0.025, 0.05, 0.1, 0.25, 1],  # with AF2
    'noAF2': [0, 0.05, 0.1, 0.25, 0.5, 1]  # without AF2
}

time_perturb_value_wAF2=[{-100: ("AF2(a=None)", 10.4), 0: ("Dicarb(a=None)", conc)} for conc in conc_dicarb['wAF2']]
sim_protocol_wAF2 = ParallelExperiments(solver, t_equil=None, time_perturb_value=time_perturb_value_wAF2)

time_perturb_value_noAF2=[{0: ("Dicarb(a=None)", conc)} for conc in conc_dicarb['noAF2']]
sim_protocol_noAF2 = ParallelExperiments(solver, t_equil=100, time_perturb_value=time_perturb_value_noAF2)

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

exp_data_file = os.path.join('Data', 'Flav_Fumarate.csv')

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      [sim_protocol_wAF2, sim_protocol_noAF2],
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5, plot_results=True,
                   plot_tc_args={'separate_plots': False, 'save_sim_data': True})
