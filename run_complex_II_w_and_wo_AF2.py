from complex_II_v3 import model
from pysb.simulator import ScipyOdeSimulator
from param_calibration import *
from SIM_PROTOCOLS.sim_protocols import SequentialInjections

model.parameters['AF2_init'].value = 0
model.parameters['BCD_init'].value = 0

solver = ScipyOdeSimulator(model)
sim_protocol_noAF2 = SimulationProtocol(solver)
sim_protocol_withAF2 = SequentialInjections(solver, perturb_day_amount={"AF2(a=None)": (0,100)})

custom_priors = {'n_Hill': ('uniform', 0.5)}
no_sample = ['A_init', 'FAD_init', 'Dicarb_init', 'AF2_init', 'AF4_init', 'BCD_init', 'k_a_fadnc_binds_bcd',
             'k_a_fadc_binds_bcd']

exp_data_file = 'Data/Complex_II_Experiment_Data_With_And_WO_AF2.csv'

if __name__ == '__main__':

    calibrator = ParameterCalibration(model,
                                      exp_data_file,
                                      [sim_protocol_withAF2, sim_protocol_noAF2],
                                      priors=custom_priors,
                                      no_sample=no_sample)

    calibrator.run(niterations=50000, nchains=5, plot_results=True)
