from complex_II_v3 import model
from pysb.simulator import ScipyOdeSimulator
from param_calibration import *

solver = ScipyOdeSimulator(model)
sim_protocol = SimulationProtocol(solver)

custom_priors = {'n_Hill': ('uniform', 0.5)}
no_sample = ['A_init', 'FAD_init', 'Dicarb_init', 'AF2_init', 'AF4_init', 'BCD_init', 'k_a_fadnc_binds_bcd',
             'k_a_fadc_binds_bcd']


exp_data_file = 'Data/Complex_II_Experiment_Data_WithAF2.csv'

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      sim_protocol,
                                      priors=custom_priors,
                                      no_sample=no_sample)
    calibrator.run(niterations=50000, nchains=5, plot_results=True)
