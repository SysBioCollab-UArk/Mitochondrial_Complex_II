from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from pysb.export.sbml import SbmlExporter

Model()

Monomer('SDHA', ['sdhaf2', 'sdhaf4'])
Monomer('SDHAF2', ['sdha'])
Monomer('SDHAF4', ['sdha'])

Parameter('SDHA_init', 100)
Parameter('SDHAF2_init', 100)
Parameter('SDHAF4_init', 100)

Initial(SDHA(sdhaf2=None, sdhaf4=None), SDHA_init)
Initial(SDHAF2(sdha=None), SDHAF2_init)
Initial(SDHAF4(sdha=None), SDHAF4_init)

Parameter('kf_sdha_free_binds_sdhaf2', 1)
Parameter('kr_sdha_free_binds_sdhaf2', 1)

Rule('SDHA_free_binds_SDHAF2',
     SDHA(sdhaf2=None, sdhaf4=None) + SDHAF2(sdha=None) | SDHA(sdhaf2=1, sdhaf4=None) % SDHAF2(sdha=1),
     kf_sdha_free_binds_sdhaf2, kr_sdha_free_binds_sdhaf2)

Parameter('kf_sdha_sdhaf4_binds_sdhaf2', 0.01)
Parameter('kr_sdha_sdhaf4_binds_sdhaf2', 100)

Rule('SDHA_SDHAF4_binds_SDHAF2',
     SDHA(sdhaf2=None, sdhaf4=ANY) + SDHAF2(sdha=None) | SDHA(sdhaf2=1, sdhaf4=ANY) % SDHAF2(sdha=1),
     kf_sdha_sdhaf4_binds_sdhaf2, kr_sdha_sdhaf4_binds_sdhaf2)

Parameter('kf_sdha_free_binds_sdhaf4', 0.1)
Parameter('kr_sdha_free_binds_sdhaf4', 0.01)

Rule('SDHA_free_binds_SDHAF4',
     SDHA(sdhaf2=None, sdhaf4=None) + SDHAF4(sdha=None) | SDHA(sdhaf2=None, sdhaf4=1) % SDHAF4(sdha=1),
     kf_sdha_free_binds_sdhaf4, kr_sdha_free_binds_sdhaf4)

Parameter('kf_sdha_sdhaf2_binds_sdhaf4', 10)
Parameter('kr_sdha_sdhaf2_binds_sdhaf4', 0.1)

Rule('SDHA_SDHAF2_binds_SDHAF4',
     SDHA(sdhaf2=ANY, sdhaf4=None) + SDHAF4(sdha=None) | SDHA(sdhaf2=ANY, sdhaf4=1) % SDHAF4(sdha=1),
     kf_sdha_sdhaf2_binds_sdhaf4, kr_sdha_sdhaf2_binds_sdhaf4)

Observable('SDHA_free', SDHA(sdhaf2=None, sdhaf4=None))
Observable('SDHAF2_free', SDHAF2(sdha=None))
Observable('SDHAF4_free', SDHAF4(sdha=None))
Observable('SDHA_SDHAF2', SDHA(sdhaf2=ANY, sdhaf4=None))
Observable('SDHA_SDHAF4', SDHA(sdhaf2=None, sdhaf4=ANY))
Observable('SDHA_SDHAF2_SDHAF4', SDHA(sdhaf2=ANY, sdhaf4=ANY))

# export model in SBML format

# sbml = SbmlExporter(model)
# with open('complexII.xml', 'w') as f:
#     print(sbml.export(level=(2,4)), file=f)

# simulation commands

tspan = np.linspace(0, 0.1, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
out = sim.run()

# for i, sp in enumerate(model.species):
#     print(i, sp)
# print()
# for i, rxn in enumerate(model.reactions):
#     print(i, rxn)

for obs in model.observables:
    plt.plot(tspan, out.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('time', fontsize=12)
plt.ylabel('amount', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='best', fontsize=12)

plt.tight_layout()
plt.show()
