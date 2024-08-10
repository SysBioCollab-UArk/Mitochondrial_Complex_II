from complex_II_v3 import model
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

tspan = np.linspace(0, 1, 101)

sim = ScipyOdeSimulator(model, tspan, verbose=True)

Result = sim.run()
print(Result.all.dtype.names)
# print(model.expressions)
for i, sp in enumerate(model.species):
    print(i, sp)
    plt.plot(tspan, Result.all["__s%d" % i], label = "s%d" % i)
plt.yscale("log")
plt.xlabel("time, min")
plt.ylabel("concentration, uM")
plt.legend(loc = 0, ncol = 2, bbox_to_anchor = (1,1))
plt.tight_layout()
plt.show()