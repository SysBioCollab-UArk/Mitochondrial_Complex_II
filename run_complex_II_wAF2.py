from complex_II_v3 import model
from pysb.simulator import ScipyOdeSimulator
from param_calibration import *
from SIM_PROTOCOLS.sim_protocols import SequentialInjections

model.parameters['BCD_init'].value = 0

solver = ScipyOdeSimulator(model)
sim_protocol = SequentialInjections(solver, t_equil=100, perturb_day_amount={"Dicarb(a=None)": (0, 10000)})

custom_priors = {'n_Hill': ('uniform', 0.5)}
no_sample = ['A_init', 'FAD_init', 'Dicarb_init', 'AF2_init', 'AF4_init', 'BCD_init', 'k_a_fadnc_binds_bcd',
             'k_a_fadc_binds_bcd']

exp_data_file = 'Data/Complex_II_Experiment_Data_WithAF2.csv'

if __name__ == '__main__':
    from pysb import Observable, ANY
    from pysb.util import alias_model_components
    import numpy as np
    import matplotlib.pyplot as plt

    alias_model_components()
    Observable('Free_SDHA', A(af2=None, af4=None, fad=None, dicarb=None))
    Observable('SDHA_SDHAF2', A(af2=ANY, af4=None, fad=None, dicarb=None))
    Observable('SDHA_FAD', A(af2=None, af4=None, fad=ANY, dicarb=None))
    Observable('SDHA_DICARB', A(af2=None, af4=None, fad=None, dicarb=ANY))
    Observable('SDHA_SDHAF2_FAD', A(af2=ANY, af4=None, fad=ANY, dicarb=None))
    Observable('SDHA_SDHA2_FAD', A(af2=ANY, af4=None, fad=None, dicarb=ANY))
    Observable('SDHA_FAD_DICARB', A(af2=None, af4=None, fad=ANY, dicarb=ANY))
    Observable('SDHA_SDHAF2_FAD_DICARB', A(af2=ANY, af4=None, fad=ANY, dicarb=ANY))
    Observable('Free_SDHAF2', AF2(a=None))
    Observable('Free_FAD', FAD(a=None))
    Observable('Free_Dicarb', Dicarb(a=None))

    Obs_to_plot = ['Free_SDHA', 'SDHA_SDHAF2', 'SDHA_FAD', 'SDHA_DICARB', 'SDHA_SDHAF2_FAD','SDHA_SDHA2_FAD',
                   'SDHA_FAD_DICARB','SDHA_SDHAF2_FAD_DICARB','Free_SDHAF2','Free_FAD','Free_Dicarb']

    tspan = np.linspace(0, 1, 101) #23, 231
    print(model.parameters)
    for ic in model.initial_conditions :
        print(ic)
    for sp in model.species:
        print(sp)
    output = solver.run(tspan)
    #output = sim_protocol.run(tspan, [p.value for p in model.parameters])
    print(output.observables["Free_SDHA"])
    print(output.observables["Free_FAD"])
    for obs in Obs_to_plot :
        plt.plot(tspan, output.observables[obs], lw=2, label = obs)
    plt.xlabel("time (min)")
    plt.ylabel("concentration (uM)")
    plt.legend(loc=0)

    plt.show()

    # calibrator = ParameterCalibration(model,
    #                                   exp_data_file,
    #                                   sim_protocol,
    #                                   priors=custom_priors,
    #                                   no_sample=no_sample)
    #
    # calibrator.run(niterations=50000, nchains=5, plot_results=True)
