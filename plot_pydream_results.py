from run_complex_II_w_and_wo_AF2 import *
import glob
import os

# Create plots using output files already generated by PyDREAM using the calibrator.create_figures() function

# get the existing PyDREAM output files
path = os.getcwd()  # path to where PyDREAM generated files are
logps_files = glob.glob(os.path.join(path, 'dreamzs*logps*'))  # need to pass ALL 'logps' files
samples_files = glob.glob(os.path.join(path, 'dreamzs*params*'))  # need to pass ALL 'params' files

# create the ParameterCalibration object (copy from the original file where PyDREAM was run from)
calibrator = ParameterCalibration(model,
                                  exp_data_file,
                                  [sim_protocol_wAF2, sim_protocol_noAF2],
                                  priors=custom_priors,
                                  no_sample=no_sample)

# call the 'create_figures' function
calibrator.create_figures(logps_files, samples_files, show_plots=True, plot_tc_args={'separate_plots': False})
