from complex_II_v3 import model
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

model.parameters['AF2_init'].value = 0

sim = ScipyOdeSimulator(model, verbose=True)

#pre-equilibration simulation
tequil = np.linspace(-10, 0, 101)
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
plt.tight_layout()

plt.show()
