from complex_II_v3 import model
from pysb.simulator import ScipyOdeSimulator
from param_calibration import *
from SIM_PROTOCOLS.sim_protocols import *

model.parameters['AF2_init'].value = 0
model.parameters['BCD_init'].value = 0

solver = ScipyOdeSimulator(model)

# Experimental protocol (Maklashina et al., 2022)
# -----------------------------------------------
# 1: Mix SDHA (6uM) + FAD with or without SDHAF2 (10.4 uM), let equilibrate
# 2: Add Dicarb at t=0 to initiate reaction
# NOTE: FAD = 100 uM in Fig. 3B and varies between 0-100 uM in Fig. 3C of Maklashina et al. (2022)
# NOTE: Dicarb = 10 mM (10,000 uM) in Fig. 3C and varies between 0-1 mM in Fig. 3B of Maklashina et al. (2022)

###### Experiments adding fumarate ######
conc_dicarb = {
    'wAF2' : np.array([0, 0.00625, 0.0125, 0.025, 0.05, 0.1, 0.25, 1]) * 1000,  # uM, with AF2
    'noAF2': np.array([0, 0.05, 0.1, 0.25, 0.5, 1]) * 1000  # uM, without AF2
}

# With AF2
tpv_dicarb_wAF2=[{
    -100: ("AF2(a=None)", 10.4),
    0: ("Dicarb(a=None)", conc)
} for conc in conc_dicarb['wAF2']]
sim_protocol_dicarb_wAF2 = ParallelExperiments(solver, t_equil=None, time_perturb_value=tpv_dicarb_wAF2)

# Without AF2
tpv_dicarb_noAF2=[{
    0: ("Dicarb(a=None)", conc)
} for conc in conc_dicarb['noAF2']]
sim_protocol_dicarb_noAF2 = ParallelExperiments(solver, t_equil=100, time_perturb_value=tpv_dicarb_noAF2)

###### Experiments adding FAD ######
conc_fad = {
    'wAF2' : [0, 1, 2, 4, 6, 20, 40, 100],  # uM, with AF2
    'noAF2': [0, 5, 10, 20, 40, 60, 100]  # uM, without AF2
}

# With AF2
tpv_fad_wAF2=[{
    -100: [("AF2(a=None)", 10.4), ("FAD(a=None, state='nc')", conc)],
    0: ("Dicarb(a=None)", 10000)
} for conc in conc_fad['wAF2']]
sim_protocol_fad_wAF2 = ParallelExperiments(solver, t_equil=None, time_perturb_value=tpv_fad_wAF2)

# Without AF2
tpv_fad_noAF2=[{
    -100: ("FAD(a=None, state='nc')", conc),
    0: ("Dicarb(a=None)", 10000)
} for conc in conc_fad['noAF2']]
sim_protocol_fad_noAF2 = ParallelExperiments(solver, t_equil=100, time_perturb_value=tpv_fad_noAF2)

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
                                      [sim_protocol_dicarb_wAF2, sim_protocol_dicarb_noAF2],
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5, plot_results=True,
                   plot_tc_args={'separate_plots': False, 'save_sim_data': True})
