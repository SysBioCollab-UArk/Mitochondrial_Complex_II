from complex_II_v3 import model
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from pysb.util import alias_model_components

alias_model_components()

data = np.genfromtxt("DATA/Complex_II_Experiment_Data_NoAF2.csv", dtype=None, delimiter=",", names=True)
print(data.dtype.names)

model.parameters['AF2_init'].value = 0

sim = ScipyOdeSimulator(model, verbose=True)

#pre-equilibration simulation
tequil = np.linspace(-100, 0, 101)
Result = sim.run(tequil)

plt.figure("species")
for i, sp in enumerate(model.species):
    print(i, sp)
    plt.plot(tequil, Result.all["__s%d" % i], label = "s%d" % i)
plt.yscale("log")
plt.ylim(bottom = 1e-5, top = 1e3)
plt.xlabel("time, min")
plt.ylabel("concentration, uM")
plt.legend(loc = 0, ncol = 2, bbox_to_anchor = (1,1))
plt.tight_layout()

plt.figure("flav")
plt.plot(tequil, Result.expressions["pct_flavinylation"], lw=2)
plt.xlabel("time, min")
plt.ylabel("flavinylation, %")
plt.tight_layout()

#add dicarboxylate

tspan = np.linspace(tequil[-1], tequil[-1] + 80, 801)
initials = Result.species[-1]
initials[2] = 10000
Result = sim.run(tspan, initials=initials)

plt.figure("species")
for i, sp in enumerate(model.species):
    plt.plot(tspan, Result.all["__s%d" % i])
plt.tight_layout()

plt.figure("flav")
plt.plot(tspan, Result.expressions["pct_flavinylation"], lw=2)
plt.errorbar(data['time'], data['average'], yerr=data['stderr'], fmt='o', capsize=6)
plt.tight_layout()

plt.figure('obs')
plt.plot(tspan, Result.observable(A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc')), lw=2,
         label='Obs_A_FADnc')
plt.plot(tspan, Result.observable(A(af2=None, af4=None, fad=None, dicarb=1) % Dicarb(a=1)), lw=2,
         label='Obs_A_Dicarb')
for obs_name in ['Obs_A_FADnc_Dicarb', 'FADc', 'SDHA_tot']:
    plt.plot(tspan, Result.observables[obs_name], lw=2, label=obs_name)
plt.xlabel("time, min")
plt.ylabel("concentration, uM")
plt.legend(loc=0)
plt.tight_layout()

plt.show()
