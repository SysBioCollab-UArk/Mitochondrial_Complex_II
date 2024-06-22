from complex_II_v3 import model
from param_calibration import *


def sim_protocol(solver,tspan,param_values):
    output = solver.run(tspan=tspan, param_values=param_values).all
    return output


custom_priors = {'n_Hill': ('uniform', 0.5)}
no_sample = ['A_init', 'FAD_init', 'Dicarb_init', 'AF2_init', 'AF4_init', 'BCD_init', 'kf_a_binds_af2',
             'kr_a_binds_af2', 'kf_a_fadnc_binds_af2', 'kr_a_fadnc_binds_af2', 'kf_a_dicarb_binds_af2',
             'kr_a_dicarb_binds_af2',  'kf_a_af2_binds_fadnc', 'kr_a_af2_binds_fadnc', 'kf_a_af2_binds_dicarb',
             'kr_a_af2_binds_dicarb', 'kf_a_af2_binds_af4', 'kr_a_af2_binds_af4', 'kf_a_af4_binds_af2',
             'kr_a_af4_binds_af2', 'kf_a_fadnc_dicarb_binds_af2', 'kr_a_fadnc_dicarb_binds_af2',
             'kf_a_fadnc_af2_binds_dicarb', 'kr_a_fadnc_af2_binds_dicarb', 'kf_a_fadnc_af2_binds_af4',
             'kr_a_fadnc_af2_binds_af4', 'kf_a_fadnc_af4_binds_af2', 'kr_a_fadnc_af4_binds_af2',
             'kf_a_dicarb_af2_binds_fadnc', 'kr_a_dicarb_af2_binds_fadnc', 'kf_a_dicarb_af2_binds_af4',
             'kr_a_dicarb_af2_binds_af4', 'kf_a_dicarb_af4_binds_af2', 'kr_a_dicarb_af4_binds_af2',
             'kf_a_af2_af4_binds_fadnc', 'kr_a_af2_af4_binds_fadnc', 'kf_a_af2_af4_binds_dicarb',
             'kr_a_af2_af4_binds_dicarb', 'k_a_fadnc_dicarb_af2_to_fadc', 'kf_a_fadc_dicarb_af2_binds_af4',
             'k_a_fadnc_dicarb_af4_binds_af2', 'k_a_fadnc_af2_af4_binds_dicarb', 'k_a_fadc_af2_unbinds_af2',
             'k_a_fadnc_binds_bcd', 'k_a_fadc_binds_bcd']


exp_data_file = 'Data/Complex_II_Experiment_Data_NoAF2.csv'

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      sim_protocol,
                                      priors=custom_priors,
                                      no_sample=no_sample)
    calibrator.run(niterations=50000, nchains=5, plot_results=True)
