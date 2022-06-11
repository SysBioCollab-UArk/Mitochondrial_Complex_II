from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt
from pysb.export.sbml import SbmlExporter

# TODO: Initial amounts and second-order rate constants

Model()

Monomer('A', ['af2', 'af4', 'fad', 'dicarb'])
Monomer('AF2', ['a'])
Monomer('AF4', ['a'])
Monomer('FAD', ['a', 'state'], {'state': ['nc', 'c']})
Monomer('Dicarb', ['a'])
Monomer('A_FADc_BCD')

Parameter('A_init', 100)
Parameter('AF2_init', 100)
Parameter('AF4_init', 100)
Parameter('FAD_init', 100)
Parameter('Dicarb_init', 100)

Initial(A(af2=None, af4=None, fad=None, dicarb=None), A_init)
Initial(AF2(a=None), AF2_init)
Initial(AF4(a=None), AF4_init)
Initial(FAD(a=None, state='nc'), FAD_init)
Initial(Dicarb(a=None), Dicarb_init)

Parameter('kf_a_binds_fad', 10)  # *
Parameter('kr_a_binds_fad', 1)
Parameter('kf_a_binds_af2', 10)  # *
Parameter('kr_a_binds_af2', 1)
Parameter('kf_a_binds_af4', 1)
Parameter('kr_a_binds_af4', 1)
Parameter('kf_a_binds_dicarb', 1)
Parameter('kr_a_binds_dicarb', 1)

# Free A binding rules

Rule('A_binds_FADnc',
     A(af2=None, af4=None, fad=None, dicarb=None) + FAD(a=None, state='nc') |
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc'),
     kf_a_binds_fad, kr_a_binds_fad)

Rule('A_binds_AF2',
     A(af2=None, af4=None, fad=None, dicarb=None) + AF2(a=None) |
     A(af2=1, af4=None, fad=None, dicarb=None) % AF2(a=1),
     kf_a_binds_af2, kr_a_binds_af2)

Rule('A_binds_AF4',
     A(af2=None, af4=None, fad=None, dicarb=None) + AF4(a=None) |
     A(af2=None, af4=1, fad=None, dicarb=None) % AF4(a=1),
     kf_a_binds_af4, kr_a_binds_af4)

Rule('A_binds_Dicarb',
     A(af2=None, af4=None, fad=None, dicarb=None) + Dicarb(a=None) |
     A(af2=None, af4=None, fad=None, dicarb=1) % Dicarb(a=1),
     kf_a_binds_dicarb, kr_a_binds_dicarb)

# A-FAD(nc) binding rules

Parameter('kf_a_fadnc_binds_dicarb', 1)
Parameter('kr_a_fadnc_binds_dicarb', 1)
Parameter('kf_a_fadnc_binds_af4', 1)
Parameter('kr_a_fadnc_binds_af4', 10)  # *
Parameter('kf_a_fadnc_binds_af2', 10)  # *
Parameter('kr_a_fadnc_binds_af2', 1)

Rule('A_FADnc_binds_Dicarb',
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') + Dicarb(a=None) |
     A(af2=None, af4=None, fad=1, dicarb=2) % FAD(a=1, state='nc') % Dicarb(a=2),
     kf_a_fadnc_binds_dicarb, kr_a_fadnc_binds_dicarb)

Rule('A_FADnc_binds_AF4',
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') + AF4(a=None) |
     A(af2=None, af4=2, fad=1, dicarb=None) % FAD(a=1, state='nc') % AF4(a=2),
     kf_a_fadnc_binds_af4, kr_a_fadnc_binds_af4)

Rule('A_FADnc_binds_AF2',
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') + AF2(a=None) |
     A(af2=2, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') % AF2(a=2),
     kf_a_fadnc_binds_af2, kr_a_fadnc_binds_af2)

# A-AF2 binding rules

Parameter('kf_a_af2_binds_fadnc', 10)  # *
Parameter('kr_a_af2_binds_fadnc', 1)

Rule('A_AF2_binds_FADnc',
     A(af2=1, af4=None, fad=None, dicarb=None) % AF2(a=1) + FAD(a=None, state='nc') |
     A(af2=1, af4=None, fad=2, dicarb=None) % AF2(a=1) % FAD(a=2, state='nc'),
     kf_a_af2_binds_fadnc, kr_a_af2_binds_fadnc)

# A-AF4 binding rules

Parameter('kf_a_af4_binds_dicarb', 1)
Parameter('kr_a_af4_binds_dicarb', 1)
Parameter('kf_a_af4_binds_af2', 1)
Parameter('kr_a_af4_binds_af2', 1)
Parameter('kf_a_af4_binds_fadnc', 1)
Parameter('kr_a_af4_binds_fadnc', 1)

Rule('A_AF4_binds_Dicarb',
     A(af2=None, af4=1, fad=None, dicarb=None) % AF4(a=1) + Dicarb(a=None) |
     A(af2=None, af4=1, fad=None, dicarb=2) % AF4(a=1) % Dicarb(a=2),
     kf_a_af4_binds_dicarb, kr_a_af4_binds_dicarb)

Rule('A_AF4_binds_AF2',
     A(af2=None, af4=1, fad=None, dicarb=None) % AF4(a=1) + AF2(a=None) |
     A(af2=2, af4=1, fad=None, dicarb=None) % AF4(a=1) % AF2(a=2),
     kf_a_af4_binds_af2, kr_a_af4_binds_af2)

Rule('A_AF4_binds_FAD',
     A(af2=None, af4=1, fad=None, dicarb=None) % AF4(a=1) + FAD(a=None, state='nc') |
     A(af2=None, af4=1, fad=2, dicarb=None) % AF4(a=1) % FAD(a=2, state='nc'),
     kf_a_af4_binds_fadnc, kr_a_af4_binds_fadnc)

# A-Dicarb rules

Parameter('k_a_dicarb_degrade', 1)

# TODO: This irreversible degradation rule is standing in for the 'A-Dicarb -> etc' reaction in the sketch
Rule('A_Dicarb_degrade',
     A(af2=None, af4=None, fad=None, dicarb=1) % Dicarb(a=1) >> None,
     k_a_dicarb_degrade)

# Additional rules

Parameter('k_a_fadnc_dicarb_to_a_fadc', 10)  # *

# TODO: Assuming the dicarboxylate is released in the reaction
Rule('A_FADnc_Dicarb_to_A_FADc',
     A(af2=None, af4=None, fad=1, dicarb=2) % FAD(a=1, state='nc') % Dicarb(a=2) >>
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c') + Dicarb(a=None),
     k_a_fadnc_dicarb_to_a_fadc)

Parameter('k_a_af2_fadnc_to_a_af2_fadnc_dicarb', 10)  # *
Parameter('k_a_af2_fadnc_to_a_af2_fadc', 10)  # *

Rule('A_AF2_FADnc_to_A_AF2_FADnc_Dicarb',
     A(af2=2, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') % AF2(a=2) + Dicarb(a=None) >>
     A(af2=2, af4=None, fad=1, dicarb=3) % FAD(a=1, state='nc') % AF2(a=2) % Dicarb(a=3),
     k_a_af2_fadnc_to_a_af2_fadnc_dicarb)

Rule('A_AF2_FADnc_to_A_AF2_FADc',
     A(af2=2, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') % AF2(a=2) >>
     A(af2=2, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c') % AF2(a=2),
     k_a_af2_fadnc_to_a_af2_fadc)

Parameter('k_a_af2_fadnc_dicarb_to_a_af2_fadc_dicarb', 10)  # *
Parameter('kf_a_af2_fadc_dicarb_to_a_af2_af4_fadc', 10)  # *
Parameter('kr_a_af2_fadc_dicarb_to_a_af2_af4_fadc', 1)
Parameter('kf_a_af2_fadc_to_a_af2_af4_fadc', 1)
Parameter('kr_a_af2_fadc_to_a_af2_af4_fadc', 1)
Parameter('kf_a_af2_af4_fadc_to_a_fadc', 100)  # *
Parameter('kr_a_af2_af4_fadc_to_a_fadc', 1)
Parameter('k_a_fadc_to_a_fadc_bcd', 1)

Rule('A_AF2_FADnc_Dicarb_to_A_AF2_FADc_Dicarb',
     A(af2=2, af4=None, fad=1, dicarb=3) % FAD(a=1, state='nc') % AF2(a=2) % Dicarb(a=3) >>
     A(af2=2, af4=None, fad=1, dicarb=3) % FAD(a=1, state='c') % AF2(a=2) % Dicarb(a=3),
     k_a_af2_fadnc_dicarb_to_a_af2_fadc_dicarb)

Rule('A_AF2_FADc_Dicarb_to_A_AF2_AF4_FADc',
     A(af2=2, af4=None, fad=1, dicarb=3) % FAD(a=1, state='c') % AF2(a=2) % Dicarb(a=3) + AF4(a=None) |
     A(af2=2, af4=3, fad=1, dicarb=None) % FAD(a=1, state='c') % AF2(a=2) % AF4(a=3) + Dicarb(a=None),
     kf_a_af2_fadc_dicarb_to_a_af2_af4_fadc, kr_a_af2_fadc_dicarb_to_a_af2_af4_fadc)

Rule('A_AF2_FADc_to_A_AF2_AF4_FADc',
     A(af2=2, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c') % AF2(a=2) + AF4(a=None) |
     A(af2=2, af4=3, fad=1, dicarb=None) % FAD(a=1, state='c') % AF2(a=2) % AF4(a=3),
     kf_a_af2_fadc_to_a_af2_af4_fadc, kr_a_af2_fadc_to_a_af2_af4_fadc)

Rule('A_AF2_AF4_FADc_to_A_FADc',
     A(af2=2, af4=3, fad=1, dicarb=None) % FAD(a=1, state='c') % AF2(a=2) % AF4(a=3) |
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c') + AF2(a=None) + AF4(a=None),
     kf_a_af2_af4_fadc_to_a_fadc, kr_a_af2_af4_fadc_to_a_fadc)

# TODO: This irreversible rule is standing in for the 'A-FAD(c) -> A-FAD(c)_BCD' reaction in the schematic
Rule('A_FADc_to_A_FADc_BCD',
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c') >> A_FADc_BCD(),
     k_a_fadc_to_a_fadc_bcd)

Observable('A_free', A(af2=None, af4=None, fad=None, dicarb=None))
# Observable('FADnc_free', FAD(a=None, state='nc'))
# Observable('AF2_free', AF2(a=None))
# Observable('AF4_free', AF4(a=None))
# Observable('Dicarb_free', Dicarb(a=None))
Observable('A_FADnc', A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc'))
Observable('A_AF2', A(af2=ANY, af4=None, fad=None, dicarb=None))
Observable('A_AF4', A(af2=None, af4=ANY, fad=None, dicarb=None))
Observable('A_Dicarb', A(af2=None, af4=None, fad=None, dicarb=ANY))
Observable('A_FADnc_Dicarb', A(af2=None, af4=None, fad=1, dicarb=ANY) % FAD(a=1, state='nc'))
Observable('A_AF2_FADnc', A(af2=ANY, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc'))
Observable('A_AF2_FADnc_Dicarb', A(af2=ANY, af4=None, fad=1, dicarb=ANY) % FAD(a=1, state='nc'))
Observable('A_AF2_FADc_Dicarb', A(af2=ANY, af4=None, fad=1, dicarb=ANY) % FAD(a=1, state='c'))
Observable('A_AF2_FADc', A(af2=ANY, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c'))
Observable('A_AF2_AF4_FADc', A(af2=ANY, af4=ANY, fad=1, dicarb=None) % FAD(a=1, state='c'))
Observable('A_FADc', A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c'))
Observable('A_FADc_BCD_total', A_FADc_BCD())

# export model in SBML format

# sbml = SbmlExporter(model)
# with open('complexII.xml', 'w') as f:
#     print(sbml.export(level=(2,4)), file=f)

# simulation commands

tspan = np.linspace(0, 1, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
out = sim.run()

# for i, sp in enumerate(model.species):
#     print(i, sp)
# print()
# for i, rxn in enumerate(model.reactions):
#     print(i, rxn)

plt.figure(figsize=(9.6, 4.8))
for i, obs in enumerate(model.observables):
     ls = '-' if i < 10 else '--'
     plt.plot(tspan, out.observables[obs.name], lw=2, ls=ls, label=obs.name)
plt.xlabel('time', fontsize=12)
plt.ylabel('amount', fontsize=12)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(loc='best', fontsize=12, bbox_to_anchor=(1.05, 1))

plt.tight_layout()
plt.show()
