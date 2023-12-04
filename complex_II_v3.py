from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

# TODO: Volume factor for 2nd-order rate constants
# TODO: Made a plot similar to Fig. 3A from Maklashina et al. (2022). Next time, do the same for Figs. 3B and 3C

Model()

Monomer('A', ['af2', 'af4', 'fad', 'dicarb'])
Monomer('FAD', ['a', 'state'], {'state': ['nc', 'c']})
Monomer('Dicarb', ['a'])
Monomer('AF2', ['a'])
Monomer('AF4', ['a'])
Monomer('BCD', ['a'])  # BCD and AF2 bind to the same site on A

Parameter('A_init', 100)
Parameter('FAD_init', 100)
Parameter('Dicarb_init', 100)
Parameter('AF2_init', 100)
Parameter('AF4_init', 100)
Parameter('BCD_init', 100)

Initial(A(af2=None, af4=None, fad=None, dicarb=None), A_init)
Initial(FAD(a=None, state='nc'), FAD_init)
Initial(Dicarb(a=None), Dicarb_init)
Initial(AF2(a=None), AF2_init)
Initial(AF4(a=None), AF4_init)
Initial(BCD(a=None), BCD_init)

# Free A binding rules

Parameter('kf_a_binds_fadnc', 1000)  #1 *
Parameter('kr_a_binds_fadnc', 1) #1000
Parameter('kf_a_binds_dicarb', 1000) #1
Parameter('kr_a_binds_dicarb', 1)
Parameter('kf_a_binds_af2', 10)  # *
Parameter('kr_a_binds_af2', 1)
Parameter('kf_a_binds_af4', 1)
Parameter('kr_a_binds_af4', 1)

# TODO: increase the binding rate for this rule to get more A_FADnc_Dicarb
Rule('A_binds_FADnc',
     A(af2=None, af4=None, fad=None, dicarb=None) + FAD(a=None, state='nc') |
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc'),
     kf_a_binds_fadnc, kr_a_binds_fadnc)

# TODO: increase the binding rate for this rule to get more A_FADnc_Dicarb
Rule('A_binds_Dicarb',
     A(af2=None, af4=None, fad=None, dicarb=None) + Dicarb(a=None) |
     A(af2=None, af4=None, fad=None, dicarb=1) % Dicarb(a=1),
     kf_a_binds_dicarb, kr_a_binds_dicarb)

Rule('A_binds_AF2',
     A(af2=None, af4=None, fad=None, dicarb=None) + AF2(a=None) |
     A(af2=1, af4=None, fad=None, dicarb=None) % AF2(a=1),
     kf_a_binds_af2, kr_a_binds_af2)

Rule('A_binds_AF4',
     A(af2=None, af4=None, fad=None, dicarb=None) + AF4(a=None) |
     A(af2=None, af4=1, fad=None, dicarb=None) % AF4(a=1),
     kf_a_binds_af4, kr_a_binds_af4)

# A-FAD(nc) binding rules

Parameter('kf_a_fadnc_binds_dicarb', 1000) #1
Parameter('kr_a_fadnc_binds_dicarb', 1)
Parameter('kf_a_fadnc_binds_af4', 1)
Parameter('kr_a_fadnc_binds_af4', 10)  # *
Parameter('kf_a_fadnc_binds_af2', 10)  # *
Parameter('kr_a_fadnc_binds_af2', 1)

# TODO: increase the binding rate for this rule to get more A_FADnc_Dicarb
Rule('A_FADnc_binds_Dicarb',
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') + Dicarb(a=None) |
     A(af2=None, af4=None, fad=1, dicarb=2) % FAD(a=1, state='nc') % Dicarb(a=2),
     kf_a_fadnc_binds_dicarb, kr_a_fadnc_binds_dicarb)

Rule('A_FADnc_binds_AF2',
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') + AF2(a=None) |
     A(af2=2, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') % AF2(a=2),
     kf_a_fadnc_binds_af2, kr_a_fadnc_binds_af2)

Rule('A_FADnc_binds_AF4',
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') + AF4(a=None) |
     A(af2=None, af4=2, fad=1, dicarb=None) % FAD(a=1, state='nc') % AF4(a=2),
     kf_a_fadnc_binds_af4, kr_a_fadnc_binds_af4)

# A-Dicarb rules

Parameter('kf_a_dicarb_binds_fadnc', 1000) #1
Parameter('kr_a_dicarb_binds_fadnc', 1) #1000
Parameter('kf_a_dicarb_binds_af2', 1)
Parameter('kr_a_dicarb_binds_af2', 1)
Parameter('kf_a_dicarb_binds_af4', 1)
Parameter('kr_a_dicarb_binds_af4', 1)

# TODO: increase the binding rate for this rule to get more A_FADnc_Dicarb
Rule('A_Dicarb_binds_FADnc',
     A(af2=None, af4=None, fad=None, dicarb=1) % Dicarb(a=1) + FAD(a=None, state='nc') |
     A(af2=None, af4=None, fad=2, dicarb=1) % Dicarb(a=1) % FAD(a=2, state='nc'),
     kf_a_dicarb_binds_fadnc, kr_a_dicarb_binds_fadnc)

Rule('A_Dicarb_binds_AF2',
     A(af2=None, af4=None, fad=None, dicarb=1) % Dicarb(a=1) + AF2(a=None) |
     A(af2=2, af4=None, fad=None, dicarb=1) % Dicarb(a=1) % AF2(a=2),
     kf_a_dicarb_binds_af2, kr_a_dicarb_binds_af2)

Rule('A_Dicarb_binds_AF4',
     A(af2=None, af4=None, fad=None, dicarb=1) % Dicarb(a=1) + AF4(a=None) |
     A(af2=None, af4=2, fad=None, dicarb=1) % Dicarb(a=1) % AF4(a=2),
     kf_a_dicarb_binds_af4, kr_a_dicarb_binds_af4)

# A-AF2 binding rules

Parameter('kf_a_af2_binds_fadnc', 1)  # *
Parameter('kr_a_af2_binds_fadnc', 1000)
Parameter('kf_a_af2_binds_dicarb', 1)
Parameter('kr_a_af2_binds_dicarb', 1)
Parameter('kf_a_af2_binds_af4', 1)
Parameter('kr_a_af2_binds_af4', 1)

Rule('A_AF2_binds_FADnc',
     A(af2=1, af4=None, fad=None, dicarb=None) % AF2(a=1) + FAD(a=None, state='nc') |
     A(af2=1, af4=None, fad=2, dicarb=None) % AF2(a=1) % FAD(a=2, state='nc'),
     kf_a_af2_binds_fadnc, kr_a_af2_binds_fadnc)

Rule('A_AF2_binds_Dicarb',
     A(af2=1, af4=None, fad=None, dicarb=None) % AF2(a=1) + Dicarb(a=None) |
     A(af2=1, af4=None, fad=None, dicarb=2) % AF2(a=1) % Dicarb(a=2),
     kf_a_af2_binds_dicarb, kr_a_af2_binds_dicarb)

Rule('A_AF2_binds_AF4',
     A(af2=1, af4=None, fad=None, dicarb=None) % AF2(a=1) + AF4(a=None) |
     A(af2=1, af4=2, fad=None, dicarb=None) % AF2(a=1) % AF4(a=2),
     kf_a_af2_binds_af4, kr_a_af2_binds_af4)

# A-AF4 binding rules

Parameter('kf_a_af4_binds_fadnc', 1)
Parameter('kr_a_af4_binds_fadnc', 1000)
Parameter('kf_a_af4_binds_dicarb', 1)
Parameter('kr_a_af4_binds_dicarb', 1)
Parameter('kf_a_af4_binds_af2', 1)
Parameter('kr_a_af4_binds_af2', 1)

Rule('A_AF4_binds_FADnc',
     A(af2=None, af4=1, fad=None, dicarb=None) % AF4(a=1) + FAD(a=None, state='nc') |
     A(af2=None, af4=1, fad=2, dicarb=None) % AF4(a=1) % FAD(a=2, state='nc'),
     kf_a_af4_binds_fadnc, kr_a_af4_binds_fadnc)

Rule('A_AF4_binds_Dicarb',
     A(af2=None, af4=1, fad=None, dicarb=None) % AF4(a=1) + Dicarb(a=None) |
     A(af2=None, af4=1, fad=None, dicarb=2) % AF4(a=1) % Dicarb(a=2),
     kf_a_af4_binds_dicarb, kr_a_af4_binds_dicarb)

Rule('A_AF4_binds_AF2',
     A(af2=None, af4=1, fad=None, dicarb=None) % AF4(a=1) + AF2(a=None) |
     A(af2=2, af4=1, fad=None, dicarb=None) % AF4(a=1) % AF2(a=2),
     kf_a_af4_binds_af2, kr_a_af4_binds_af2)

# A-FADnc-Dicarb binding rules

Parameter('kf_a_fadnc_dicarb_binds_af2', 1)
Parameter('kr_a_fadnc_dicarb_binds_af2', 1)
Parameter('kf_a_fadnc_dicarb_binds_af4', 1)
Parameter('kr_a_fadnc_dicarb_binds_af4', 1)

Rule('A_FADnc_Dicarb_binds_AF2',
     A(af2=None, af4=None, fad=1, dicarb=ANY) % FAD(a=1, state='nc') + AF2(a=None) |
     A(af2=2, af4=None, fad=1, dicarb=ANY) % FAD(a=1, state='nc') % AF2(a=2),
     kf_a_fadnc_dicarb_binds_af2, kr_a_fadnc_dicarb_binds_af2)

Rule('A_FADnc_Dicarb_binds_AF4',
     A(af2=None, af4=None, fad=1, dicarb=ANY) % FAD(a=1, state='nc') + AF4(a=None) |
     A(af2=None, af4=2, fad=1, dicarb=ANY) % FAD(a=1, state='nc') % AF4(a=2),
     kf_a_fadnc_dicarb_binds_af4, kr_a_fadnc_dicarb_binds_af4)

# A-FADnc-AF2 binding rules

Parameter('kf_a_fadnc_af2_binds_dicarb', 1)
Parameter('kr_a_fadnc_af2_binds_dicarb', 1)
Parameter('kf_a_fadnc_af2_binds_af4', 1)
Parameter('kr_a_fadnc_af2_binds_af4', 1)

Rule('A_FADnc_AF2_binds_Dicarb',
     A(af2=3, af4=None, fad=1, dicarb=None) % AF2(a=3) % FAD(a=1, state='nc') + Dicarb(a=None) |
     A(af2=3, af4=None, fad=1, dicarb=2) % AF2(a=3) % FAD(a=1, state='nc') % Dicarb(a=2),
     kf_a_fadnc_af2_binds_dicarb, kr_a_fadnc_af2_binds_dicarb)

Rule('A_FADnc_AF2_binds_AF4',
     A(af2=3, af4=None, fad=1, dicarb=None) % AF2(a=3) % FAD(a=1, state='nc') + AF4(a=None) |
     A(af2=3, af4=2, fad=1, dicarb=None) % AF2(a=3) % FAD(a=1, state='nc') % AF4(a=2),
     kf_a_fadnc_af2_binds_af4, kr_a_fadnc_af2_binds_af4)

# A-FADnc-AF4 binding rules

Parameter('kf_a_fadnc_af4_binds_dicarb', 1)
Parameter('kr_a_fadnc_af4_binds_dicarb', 1)
Parameter('kf_a_fadnc_af4_binds_af2', 1)
Parameter('kr_a_fadnc_af4_binds_af2', 1)

Rule('A_FADnc_AF4_binds_Dicarb',
     A(af2=None, af4=ANY, fad=1, dicarb=None) % FAD(a=1, state='nc') + Dicarb(a=None) |
     A(af2=None, af4=ANY, fad=1, dicarb=2) % FAD(a=1, state='nc') % Dicarb(a=2),
     kf_a_fadnc_af4_binds_dicarb, kr_a_fadnc_af4_binds_dicarb)

Rule('A_FADnc_AF4_binds_AF2',
     A(af2=None, af4=ANY, fad=1, dicarb=None) % FAD(a=1, state='nc') + AF2(a=None) |
     A(af2=2, af4=ANY, fad=1, dicarb=None) % FAD(a=1, state='nc') % AF2(a=2),
     kf_a_fadnc_af4_binds_af2, kr_a_fadnc_af4_binds_af2)

# A-Dicarb-AF2 binding rules

Parameter('kf_a_dicarb_af2_binds_fadnc', 1)
Parameter('kr_a_dicarb_af2_binds_fadnc', 1)
Parameter('kf_a_dicarb_af2_binds_af4', 1)
Parameter('kr_a_dicarb_af2_binds_af4', 1)

Rule('A_Dicarb_AF2_binds_FADnc',
     A(af2=2, af4=None, fad=None, dicarb=ANY) % AF2(a=2) + FAD(a=None, state='nc') |
     A(af2=2, af4=None, fad=1, dicarb=ANY) % AF2(a=2) % FAD(a=1, state='nc'),
     kf_a_dicarb_af2_binds_fadnc, kr_a_dicarb_af2_binds_fadnc)

Rule('A_Dicarb_AF2_binds_AF4',
     A(af2=2, af4=None, fad=None, dicarb=ANY) % AF2(a=2) + AF4(a=None) |
     A(af2=2, af4=1, fad=None, dicarb=ANY) % AF2(a=2) % AF4(a=1),
     kf_a_dicarb_af2_binds_af4, kr_a_dicarb_af2_binds_af4)

# A-Dicarb-AF4 binding rules

Parameter('kf_a_dicarb_af4_binds_fadnc', 1)
Parameter('kr_a_dicarb_af4_binds_fadnc', 1)
Parameter('kf_a_dicarb_af4_binds_af2', 1)
Parameter('kr_a_dicarb_af4_binds_af2', 1)

Rule('A_Dicarb_AF4_binds_FADnc',
     A(af2=None, af4=ANY, fad=None, dicarb=ANY) + FAD(a=None, state='nc') |
     A(af2=None, af4=ANY, fad=1, dicarb=ANY) % FAD(a=1, state='nc'),
     kf_a_dicarb_af4_binds_fadnc, kr_a_dicarb_af4_binds_fadnc)

Rule('A_Dicarb_AF4_binds_AF2',
     A(af2=None, af4=ANY, fad=None, dicarb=ANY) + AF2(a=None) |
     A(af2=1, af4=ANY, fad=None, dicarb=ANY) % AF2(a=1),
     kf_a_dicarb_af4_binds_af2, kr_a_dicarb_af4_binds_af2)

# A-AF2-AF4 binding rules

Parameter('kf_a_af2_af4_binds_fadnc', 1)
Parameter('kr_a_af2_af4_binds_fadnc', 1)
Parameter('kf_a_af2_af4_binds_dicarb', 1)
Parameter('kr_a_af2_af4_binds_dicarb', 1)

Rule('A_AF2_AF4_binds_FADnc',
     A(af2=2, af4=ANY, fad=None, dicarb=None) % AF2(a=2) + FAD(a=None, state='nc') |
     A(af2=2, af4=ANY, fad=1, dicarb=None) % AF2(a=2) % FAD(a=1, state='nc'),
     kf_a_af2_af4_binds_fadnc, kr_a_af2_af4_binds_fadnc)

Rule('A_AF2_AF4_binds_Dicarb',
     A(af2=2, af4=ANY, fad=None, dicarb=None) % AF2(a=2) + Dicarb(a=None) |
     A(af2=2, af4=ANY, fad=None, dicarb=1) % AF2(a=2) % Dicarb(a=1),
     kf_a_af2_af4_binds_dicarb, kr_a_af2_af4_binds_dicarb)

##########

# A-FADnc-Dicarb state change

# TODO: This is the rate constant to change to model flavinylation in the absence of AF2
#  (SdhE in Maklashina et al. 2022)
#Parameter('k_a_fadnc_dicarb_to_fadc', 5)

# NOTE: Dicarb dissociates from the complex upon FAD state change (thin black arrow on schematic)
# Rule('A_FADnc_Dicarb_to_FADc',
#      A(af2=None, af4=None, fad=1, dicarb=2) % FAD(a=1, state='nc') % Dicarb(a=2) >>
#      A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c') + Dicarb(a=None),
#      k_a_fadnc_dicarb_to_fadc)

Parameter('kcat_a_fadnc_dicarb_to_fadc', 1E6)
Parameter('Km_a_fadnc_dicarb_to_fadc', 60)
Parameter('n_Hill', 4)
Observable("Obs_A_FADnc_Dicarb", A(af2=None, af4=None, fad=1, dicarb=2) % FAD(a=1, state='nc') % Dicarb(a=2))
Expression("rate_A_FADnc_Dicarb_to_FADc",
           kcat_a_fadnc_dicarb_to_fadc * Obs_A_FADnc_Dicarb ** (n_Hill-1) /
           (Km_a_fadnc_dicarb_to_fadc ** n_Hill + Obs_A_FADnc_Dicarb ** n_Hill))
# NOTE: Dicarb dissociates from the complex upon FAD state change (thin black arrow on schematic)
Rule('A_FADnc_Dicarb_to_FADc',
     A(af2=None, af4=None, fad=1, dicarb=2) % FAD(a=1, state='nc') % Dicarb(a=2) >>
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c') + Dicarb(a=None),
     rate_A_FADnc_Dicarb_to_FADc)


# A-FADnc-Dicarb-AF2 state change

Parameter('k_a_fadnc_dicarb_af2_to_fadc', 1)

Rule('A_FADnc_Dicarb_AF2_to_FADc',
     A(af2=2, af4=None, fad=1, dicarb=ANY) % AF2(a=2) % FAD(a=1, state='nc') >>
     A(af2=2, af4=None, fad=1, dicarb=ANY) % AF2(a=2) % FAD(a=1, state='c'),
     k_a_fadnc_dicarb_af2_to_fadc)

# A-FADc-Dicarb-AF2 binds AF4
# NOTE: unbinding rule is defined later

Parameter('kf_a_fadc_dicarb_af2_binds_af4', 1)

Rule('A_FADc_Dicarb_AF2_binds_AF4',
     A(af2=3, af4=None, fad=1, dicarb=ANY) % AF2(a=3) % FAD(a=1, state='c') + AF4(a=None) >>
     A(af2=3, af4=2, fad=1, dicarb=ANY) % AF2(a=3) % FAD(a=1, state='c') % AF4(a=2),
     kf_a_fadc_dicarb_af2_binds_af4)

##########

# A-FADnc-Dicarb-AF4 binding + state change rule

Parameter('k_a_fadnc_dicarb_af4_binds_af2', 1)

Rule('A_FADnc_Dicarb_AF4_binds_AF2',
     A(af2=None, af4=ANY, fad=1, dicarb=ANY) % FAD(a=1, state='nc') + AF2(a=None) >>
     A(af2=2, af4=ANY, fad=1, dicarb=ANY) % FAD(a=1, state='c') % AF2(a=2),
     k_a_fadnc_dicarb_af4_binds_af2)

# A-FADnc-AF2-AF4 binding + state change rule

Parameter('k_a_fadnc_af2_af4_binds_dicarb', 1)

Rule('A_FADnc_AF2_AF4_binds_Dicarb',
     A(af2=1, af4=ANY, fad=2, dicarb=None) % AF2(a=1) % FAD(a=2, state='nc') + Dicarb(a=None) >>
     A(af2=1, af4=ANY, fad=2, dicarb=3) % AF2(a=1) % FAD(a=2, state='c') % Dicarb(a=3),
     k_a_fadnc_af2_af4_binds_dicarb)

##########

# A-FADc unbinding rules

Parameter('k_a_fadc_dicarb_unbinds_dicarb', 1)
Parameter('k_a_fadc_af2_unbinds_af2', 1)
Parameter('k_a_fadc_af4_unbinds_af4', 1)

# NOTE: Sequence of unbinding events is currently unknown
# Currently modeling dissociation events as independent - May change this in the future

Rule('A_FADc_Dicarb_unbinds_Dicarb',
     A(fad=1, dicarb=2) % FAD(a=1, state='c') % Dicarb(a=2) >>
     A(fad=1, dicarb=None) % FAD(a=1, state='c') + Dicarb(a=None),
     k_a_fadc_dicarb_unbinds_dicarb)

Rule('A_FADc_AF2_unbinds_AF2',
     A(af2=2, fad=1) % FAD(a=1, state='c') % AF2(a=2) >>
     A(af2=None, fad=1) % FAD(a=1, state='c') + AF2(a=None),
     k_a_fadc_af2_unbinds_af2)

Rule('A_FADc_AF4_unbinds_AF4',
     A(af4=2, fad=1) % FAD(a=1, state='c') % AF4(a=2) >>
     A(af4=None, fad=1) % FAD(a=1, state='c') + AF4(a=None),
     k_a_fadc_af4_unbinds_af4)

##########

# Formation of partially active CII (w/ non-covalent FAD)

Parameter('k_a_fadnc_af2_binds_bcd', 1)
Parameter('k_a_fadnc_af4_binds_bcd', 1)

Rule('A_FADnc_binds_BCD',
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') + BCD(a=None) >>
     A(af2=2, af4=None, fad=1, dicarb=None) % FAD(a=1, state='nc') % BCD(a=2),
     k_a_fadnc_af2_binds_bcd)

# Rule('A_FADnc_AF2_binds_BCD',
#      A(af2=ANY, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + BCD(a=None) >>
#      A(af2=ANY, af4=None, fad=1, dicarb=None, bcd=2) % FAD(a=1, state='nc') % BCD(a=2),
#      k_a_fadnc_af2_binds_bcd)

# Rule('A_FADnc_AF4_binds_BCD',
#      A(af2=None, af4=ANY, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + BCD(a=None) >>
#      A(af2=None, af4=ANY, fad=1, dicarb=None, bcd=2) % FAD(a=1, state='nc') % BCD(a=2),
#      k_a_fadnc_af4_binds_bcd)

# Formation of active CII

Parameter('k_a_fadc_dicarb_binds_bcd', 1)
Parameter('k_a_fadc_binds_bcd', 1)

# Rule('A_FADc_Dicarb_binds_BCD',
#      A(af2=None, af4=None, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='c') + BCD(a=None) >>
#      A(af2=None, af4=None, fad=1, dicarb=ANY, bcd=2) % FAD(a=1, state='c') % BCD(a=2),
#      k_a_fadc_dicarb_binds_bcd)

Rule('A_FADc_binds_BCD',
     A(af2=None, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c') + BCD(a=None) >>
     A(af2=2, af4=None, fad=1, dicarb=None) % FAD(a=1, state='c') % BCD(a=2),
     k_a_fadc_binds_bcd)

if __name__ == "__main__" :

     # Observables for various complexes within the pathway
     obs_to_plot = [
          Observable('Free_SDHA', A(af2=None, af4=None, fad=None, dicarb=None)),
          Observable('SDHA_without_FAD', A(fad=None)),
          Observable('SDHA_FADnc', A(fad=1) % FAD(a=1, state='nc')),
          Observable('SDHA_FADc', A(fad=1) % FAD(a=1, state='c')),
          Observable('CII_with_FADnc', A(fad=1, af2=2) % FAD(a=1, state='nc') % BCD(a=2)),
          Observable('Active_CII', A(fad=1, af2=2) % FAD(a=1, state='c') % BCD(a=2))
     ]

     # Observables for creating plots like in Maklashina et al. (2022)
     Observable('FADc', FAD(state='c'))
     Observable('FAD_tot', FAD())

     # simulation commands

     tspan = np.linspace(0, 30, 301)
     sim = ScipyOdeSimulator(model, tspan, verbose=True)
     out = sim.run(initials={BCD(a=None): 0})

     # Plot various complexes within the pathway
     for obs in obs_to_plot:
         plt.plot(tspan, out.observables[obs.name], lw=2, label=obs.name)
     plt.xlabel('time')
     plt.ylabel('amount')
     plt.legend(loc='best')
     plt.tight_layout()

     # Plot % flavinylation
     plt.figure('flavinylation')
     plt.plot(tspan, out.observables['FADc'] / out.observables['FAD_tot'] * 100, color='b', lw=2, label='+ AF2')
     plt.ylim(top=100)
     plt.xlabel('time')
     plt.ylabel('flavinylation, %')
     plt.tight_layout()

     # Change AF2 concentration to zero and run again
     out = sim.run(initials={BCD(a=None): 0, AF2(a=None): 0})
     plt.figure('flavinylation')
     plt.plot(tspan, out.observables['FADc'] / out.observables['FAD_tot'] * 100, color='k', lw=2, label='- AF2')
     plt.legend(loc=0)

     # Run Simulations with different initial amounts of Dicarb.
     plt.figure('dicarb')
     # Todo: loop over different initial amounts of dicarb and run simulations with and without AF2
     tspan = np.linspace(0, 30, 301)
     flav1 = []
     flav2 = []
     dicarb = np.append(np.arange(0,10,1),np.arange(10,100,10))
     for d in dicarb :
          out1 = sim.run( tspan = tspan, initials={BCD(a=None): 0, Dicarb(a=None): d})
          #plt.plot(tspan, out1.observables['FADc'] / out1.observables['FAD_tot'] * 100, color='b', lw=2, label='+ AF2')
          flav1.append(out1.observables['FADc'][-1] / out1.observables['FAD_tot'][-1] * 100)
          out2 = sim.run( tspan = tspan, initials={BCD(a=None): 0, Dicarb(a=None): d, AF2(a=None): 0})
          flav2.append(out2.observables['FADc'][-1] / out2.observables['FAD_tot'][-1] * 100)
     plt.plot(dicarb, flav1, "-o", color='b', lw=2, label='+ AF2')
     plt.plot(dicarb, flav2, "-o", color='k', lw=2, label='- AF2')
     plt.xlabel('dicarboxylate, conc')
     plt.ylabel('flavinylation, %')
     plt.legend(loc=0)
     plt.tight_layout()


     # Run Simulations with different initial amounts of FAD.
     plt.figure('FAD')
     # Todo: loop over different initial amounts of dicarb and run simulations with and without AF2
     tspan = np.linspace(0, 30, 301)
     flav1 = []
     flav2 = []
     fadc = []
     fadtot = []
     fad = np.append(np.arange(0, 10, 1), np.arange(10, 100, 10))
     for f in fad:
          out1 = sim.run(tspan=tspan, initials={BCD(a=None): 0, FAD(a=None, state="nc"): f})
          # plt.plot(tspan, out1.observables['FADc'] / out1.observables['FAD_tot'] * 100, color='b', lw=2, label='+ AF2')
          flav1.append(out1.observables['FADc'][-1] / out1.observables['FAD_tot'][-1] * 100)
          out2 = sim.run(tspan=tspan, initials={BCD(a=None): 0, FAD(a=None, state="nc"): f, AF2(a=None): 0})
          ####
          plt.figure("rate_exp")
          #plt.plot(tspan, out2.expressions["rate_A_FADnc_Dicarb_to_FADc"] * out2.observables["Obs_A_FADnc_Dicarb"], lw = 2, label = "fad = %g" % f)
          plt.plot(tspan,out2.observables["Obs_A_FADnc_Dicarb"], lw = 2, label = "fad = %g" % f)
          plt.legend(loc = 0)
          plt.figure('FAD')
          ####
          flav2.append(out2.observables['FADc'][-1] / out2.observables['FAD_tot'][-1] * 100)
          fadc.append(out2.observables['FADc'][-1])
          fadtot.append(out2.observables['FAD_tot'][-1])
     plt.plot(fad, flav1, "-o", color='b', lw=2, label='+ AF2')
     plt.plot(fad, flav2, "-o", color='k', lw=2, label='- AF2')
     plt.xlabel('FAD, conc')
     plt.ylabel('flavinylation, %')
     plt.ylim(bottom = -5, top = 105)
     plt.legend(loc=0)
     plt.tight_layout()

     plt.figure()
     plt.plot(fad, fadc, label = "fadc")
     plt.plot(fad, fadtot, label = "fadtot")
     plt.legend(loc=0)



     plt.show()
