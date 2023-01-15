from pysb import *
from pysb.simulator import ScipyOdeSimulator
import numpy as np
import matplotlib.pyplot as plt

# TODO: Volume factor for 2nd-order rate constants

Model()

Monomer('A', ['af2', 'af4', 'fad', 'dicarb', 'bcd'])
Monomer('FAD', ['a', 'state'], {'state': ['nc', 'c']})
Monomer('Dicarb', ['a'])
Monomer('AF2', ['a'])
Monomer('AF4', ['a'])
# TODO: BCD and AF2 share a binding site to A
Monomer('BCD', ['a'])

# Monomer('A_FADc_BCD')

Parameter('A_init', 100)
Parameter('FAD_init', 100)
Parameter('Dicarb_init', 100)
Parameter('AF2_init', 100)
Parameter('AF4_init', 100)
Parameter('BCD_init', 100)

Initial(A(af2=None, af4=None, fad=None, dicarb=None, bcd=None), A_init)
Initial(FAD(a=None, state='nc'), FAD_init)
Initial(Dicarb(a=None), Dicarb_init)
Initial(AF2(a=None), AF2_init)
Initial(AF4(a=None), AF4_init)
Initial(BCD(a=None), BCD_init)

# Free A binding rules

Parameter('kf_a_binds_fadnc', 10)  # *
Parameter('kr_a_binds_fadnc', 1)
Parameter('kf_a_binds_dicarb', 1)
Parameter('kr_a_binds_dicarb', 1)
Parameter('kf_a_binds_af2', 10)  # *
Parameter('kr_a_binds_af2', 1)
Parameter('kf_a_binds_af4', 1)
Parameter('kr_a_binds_af4', 1)

Rule('A_binds_FADnc',
     A(af2=None, af4=None, fad=None, dicarb=None, bcd=None) + FAD(a=None, state='nc') |
     A(af2=None, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc'),
     kf_a_binds_fadnc, kr_a_binds_fadnc)

Rule('A_binds_Dicarb',
     A(af2=None, af4=None, fad=None, dicarb=None, bcd=None) + Dicarb(a=None) |
     A(af2=None, af4=None, fad=None, dicarb=1, bcd=None) % Dicarb(a=1),
     kf_a_binds_dicarb, kr_a_binds_dicarb)

Rule('A_binds_AF2',
     A(af2=None, af4=None, fad=None, dicarb=None, bcd=None) + AF2(a=None) |
     A(af2=1, af4=None, fad=None, dicarb=None, bcd=None) % AF2(a=1),
     kf_a_binds_af2, kr_a_binds_af2)

Rule('A_binds_AF4',
     A(af2=None, af4=None, fad=None, dicarb=None, bcd=None) + AF4(a=None) |
     A(af2=None, af4=1, fad=None, dicarb=None, bcd=None) % AF4(a=1),
     kf_a_binds_af4, kr_a_binds_af4)

# A-FAD(nc) binding rules

Parameter('kf_a_fadnc_binds_dicarb', 1)
Parameter('kr_a_fadnc_binds_dicarb', 1)
Parameter('kf_a_fadnc_binds_af4', 1)
Parameter('kr_a_fadnc_binds_af4', 10)  # *
Parameter('kf_a_fadnc_binds_af2', 10)  # *
Parameter('kr_a_fadnc_binds_af2', 1)

Rule('A_FADnc_binds_Dicarb',
     A(af2=None, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + Dicarb(a=None) |
     A(af2=None, af4=None, fad=1, dicarb=2, bcd=None) % FAD(a=1, state='nc') % Dicarb(a=2),
     kf_a_fadnc_binds_dicarb, kr_a_fadnc_binds_dicarb)

Rule('A_FADnc_binds_AF2',
     A(af2=None, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + AF2(a=None) |
     A(af2=2, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') % AF2(a=2),
     kf_a_fadnc_binds_af2, kr_a_fadnc_binds_af2)

Rule('A_FADnc_binds_AF4',
     A(af2=None, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + AF4(a=None) |
     A(af2=None, af4=2, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') % AF4(a=2),
     kf_a_fadnc_binds_af4, kr_a_fadnc_binds_af4)

# A-Dicarb rules

Parameter('kf_a_dicarb_binds_fadnc', 1)
Parameter('kr_a_dicarb_binds_fadnc', 1)
Parameter('kf_a_dicarb_binds_af2', 1)
Parameter('kr_a_dicarb_binds_af2', 1)
Parameter('kf_a_dicarb_binds_af4', 1)
Parameter('kr_a_dicarb_binds_af4', 1)

Rule('A_Dicarb_binds_FADnc',
     A(af2=None, af4=None, fad=None, dicarb=1, bcd=None) % Dicarb(a=1) + FAD(a=None, state='nc') |
     A(af2=None, af4=None, fad=2, dicarb=1, bcd=None) % Dicarb(a=1) % FAD(a=2, state='nc'),
     kf_a_dicarb_binds_fadnc, kr_a_dicarb_binds_fadnc)

Rule('A_Dicarb_binds_AF2',
     A(af2=None, af4=None, fad=None, dicarb=1, bcd=None) % Dicarb(a=1) + AF2(a=None) |
     A(af2=2, af4=None, fad=None, dicarb=1, bcd=None) % Dicarb(a=1) % AF2(a=2),
     kf_a_dicarb_binds_af2, kr_a_dicarb_binds_af2)

Rule('A_Dicarb_binds_AF4',
     A(af2=None, af4=None, fad=None, dicarb=1, bcd=None) % Dicarb(a=1) + AF4(a=None) |
     A(af2=None, af4=2, fad=None, dicarb=1, bcd=None) % Dicarb(a=1) % AF4(a=2),
     kf_a_dicarb_binds_af4, kr_a_dicarb_binds_af4)

# A-AF2 binding rules

Parameter('kf_a_af2_binds_fadnc', 10)  # *
Parameter('kr_a_af2_binds_fadnc', 1)
Parameter('kf_a_af2_binds_dicarb', 1)
Parameter('kr_a_af2_binds_dicarb', 1)
Parameter('kf_a_af2_binds_af4', 1)
Parameter('kr_a_af2_binds_af4', 1)

Rule('A_AF2_binds_FADnc',
     A(af2=1, af4=None, fad=None, dicarb=None, bcd=None) % AF2(a=1) + FAD(a=None, state='nc') |
     A(af2=1, af4=None, fad=2, dicarb=None, bcd=None) % AF2(a=1) % FAD(a=2, state='nc'),
     kf_a_af2_binds_fadnc, kr_a_af2_binds_fadnc)

Rule('A_AF2_binds_Dicarb',
     A(af2=1, af4=None, fad=None, dicarb=None, bcd=None) % AF2(a=1) + Dicarb(a=None) |
     A(af2=1, af4=None, fad=None, dicarb=2, bcd=None) % AF2(a=1) % Dicarb(a=2),
     kf_a_af2_binds_dicarb, kr_a_af2_binds_dicarb)

Rule('A_AF2_binds_AF4',
     A(af2=1, af4=None, fad=None, dicarb=None, bcd=None) % AF2(a=1) + AF4(a=None) |
     A(af2=1, af4=2, fad=None, dicarb=None, bcd=None) % AF2(a=1) % AF4(a=2),
     kf_a_af2_binds_af4, kr_a_af2_binds_af4)

# A-AF4 binding rules

Parameter('kf_a_af4_binds_fadnc', 1)
Parameter('kr_a_af4_binds_fadnc', 1)
Parameter('kf_a_af4_binds_dicarb', 1)
Parameter('kr_a_af4_binds_dicarb', 1)
Parameter('kf_a_af4_binds_af2', 1)
Parameter('kr_a_af4_binds_af2', 1)

Rule('A_AF4_binds_FADnc',
     A(af2=None, af4=1, fad=None, dicarb=None, bcd=None) % AF4(a=1) + FAD(a=None, state='nc') |
     A(af2=None, af4=1, fad=2, dicarb=None, bcd=None) % AF4(a=1) % FAD(a=2, state='nc'),
     kf_a_af4_binds_fadnc, kr_a_af4_binds_fadnc)

Rule('A_AF4_binds_Dicarb',
     A(af2=None, af4=1, fad=None, dicarb=None, bcd=None) % AF4(a=1) + Dicarb(a=None) |
     A(af2=None, af4=1, fad=None, dicarb=2, bcd=None) % AF4(a=1) % Dicarb(a=2),
     kf_a_af4_binds_dicarb, kr_a_af4_binds_dicarb)

Rule('A_AF4_binds_AF2',
     A(af2=None, af4=1, fad=None, dicarb=None, bcd=None) % AF4(a=1) + AF2(a=None) |
     A(af2=2, af4=1, fad=None, dicarb=None, bcd=None) % AF4(a=1) % AF2(a=2),
     kf_a_af4_binds_af2, kr_a_af4_binds_af2)

# A-FADnc-Dicarb binding rules

Parameter('kf_a_fadnc_dicarb_binds_af2', 1)
Parameter('kr_a_fadnc_dicarb_binds_af2', 1)
Parameter('kf_a_fadnc_dicarb_binds_af4', 1)
Parameter('kr_a_fadnc_dicarb_binds_af4', 1)

Rule('A_FADnc_Dicarb_binds_AF2',
     A(af2=None, af4=None, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='nc') + AF2(a=None) |
     A(af2=2, af4=None, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='nc') % AF2(a=2),
     kf_a_fadnc_dicarb_binds_af2, kr_a_fadnc_dicarb_binds_af2)

Rule('A_FADnc_Dicarb_binds_AF4',
     A(af2=None, af4=None, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='nc') + AF4(a=None) |
     A(af2=None, af4=2, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='nc') % AF4(a=2),
     kf_a_fadnc_dicarb_binds_af4, kr_a_fadnc_dicarb_binds_af4)

# A-FADnc-AF2 binding rules

Parameter('kf_a_fadnc_af2_binds_dicarb', 1)
Parameter('kr_a_fadnc_af2_binds_dicarb', 1)
Parameter('kf_a_fadnc_af2_binds_af4', 1)
Parameter('kr_a_fadnc_af2_binds_af4', 1)

Rule('A_FADnc_AF2_binds_Dicarb',
     A(af2=ANY, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + Dicarb(a=None) |
     A(af2=ANY, af4=None, fad=1, dicarb=2, bcd=None) % FAD(a=1, state='nc') % Dicarb(a=2),
     kf_a_fadnc_af2_binds_dicarb, kr_a_fadnc_af2_binds_dicarb)

Rule('A_FADnc_AF2_binds_AF4',
     A(af2=ANY, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + AF4(a=None) |
     A(af2=ANY, af4=2, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') % AF4(a=2),
     kf_a_fadnc_af2_binds_af4, kr_a_fadnc_af2_binds_af4)

# A-FADnc-AF4 binding rules

Parameter('kf_a_fadnc_af4_binds_dicarb', 1)
Parameter('kr_a_fadnc_af4_binds_dicarb', 1)
Parameter('kf_a_fadnc_af4_binds_af2', 1)
Parameter('kr_a_fadnc_af4_binds_af2', 1)

Rule('A_FADnc_AF4_binds_Dicarb',
     A(af2=None, af4=ANY, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + Dicarb(a=None) |
     A(af2=None, af4=ANY, fad=1, dicarb=2, bcd=None) % FAD(a=1, state='nc') % Dicarb(a=2),
     kf_a_fadnc_af4_binds_dicarb, kr_a_fadnc_af4_binds_dicarb)

Rule('A_FADnc_AF4_binds_AF2',
     A(af2=None, af4=ANY, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + AF2(a=None) |
     A(af2=2, af4=ANY, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') % AF2(a=2),
     kf_a_fadnc_af4_binds_af2, kr_a_fadnc_af4_binds_af2)

# A-Dicarb-AF2 binding rules

Parameter('kf_a_dicarb_af2_binds_fadnc', 1)
Parameter('kr_a_dicarb_af2_binds_fadnc', 1)
Parameter('kf_a_dicarb_af2_binds_af4', 1)
Parameter('kr_a_dicarb_af2_binds_af4', 1)

Rule('A_Dicarb_AF2_binds_FADnc',
     A(af2=ANY, af4=None, fad=None, dicarb=ANY, bcd=None) + FAD(a=None, state='nc') |
     A(af2=ANY, af4=None, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='nc'),
     kf_a_dicarb_af2_binds_fadnc, kr_a_dicarb_af2_binds_fadnc)

Rule('A_Dicarb_AF2_binds_AF4',
     A(af2=ANY, af4=None, fad=None, dicarb=ANY, bcd=None) + AF4(a=None) |
     A(af2=ANY, af4=1, fad=None, dicarb=ANY, bcd=None) % AF4(a=1),
     kf_a_dicarb_af2_binds_af4, kr_a_dicarb_af2_binds_af4)

# A-Dicarb-AF4 binding rules

Parameter('kf_a_dicarb_af4_binds_fadnc', 1)
Parameter('kr_a_dicarb_af4_binds_fadnc', 1)
Parameter('kf_a_dicarb_af4_binds_af2', 1)
Parameter('kr_a_dicarb_af4_binds_af2', 1)

Rule('A_Dicarb_AF4_binds_FADnc',
     A(af2=None, af4=ANY, fad=None, dicarb=ANY, bcd=None) + FAD(a=None, state='nc') |
     A(af2=None, af4=ANY, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='nc'),
     kf_a_dicarb_af4_binds_fadnc, kr_a_dicarb_af4_binds_fadnc)

Rule('A_Dicarb_AF4_binds_AF2',
     A(af2=None, af4=ANY, fad=None, dicarb=ANY, bcd=None) + AF2(a=None) |
     A(af2=1, af4=ANY, fad=None, dicarb=ANY, bcd=None) % AF2(a=1),
     kf_a_dicarb_af4_binds_af2, kr_a_dicarb_af4_binds_af2)

# A-AF2-AF4 binding rules

Parameter('kf_a_af2_af4_binds_fadnc', 1)
Parameter('kr_a_af2_af4_binds_fadnc', 1)
Parameter('kf_a_af2_af4_binds_dicarb', 1)
Parameter('kr_a_af2_af4_binds_dicarb', 1)

Rule('A_AF2_AF4_binds_FADnc',
     A(af2=ANY, af4=ANY, fad=None, dicarb=None, bcd=None) + FAD(a=None, state='nc') |
     A(af2=ANY, af4=ANY, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc'),
     kf_a_af2_af4_binds_fadnc, kr_a_af2_af4_binds_fadnc)

Rule('A_AF2_AF4_binds_Dicarb',
     A(af2=ANY, af4=ANY, fad=None, dicarb=None, bcd=None) + Dicarb(a=None) |
     A(af2=ANY, af4=ANY, fad=None, dicarb=1, bcd=None) % Dicarb(a=1),
     kf_a_af2_af4_binds_dicarb, kr_a_af2_af4_binds_dicarb)

##########

# A-FADnc-Dicarb state change

Parameter('k_a_fadnc_dicarb_to_fadc', 1)

# NOTE: Dicarb dissociates from the complex upon FAD state change (thin black arrow on schematic)
Rule('A_FADnc_Dicarb_to_FADc',
     A(af2=None, af4=None, fad=1, dicarb=2, bcd=None) % FAD(a=1, state='nc') % Dicarb(a=2) >>
     A(af2=None, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='c') + Dicarb(a=None),
     k_a_fadnc_dicarb_to_fadc)

# A-FADnc-Dicarb-AF2 state change

Parameter('k_a_fadnc_dicarb_af2_to_fadc', 1)

Rule('A_FADnc_Dicarb_AF2_to_FADc',
     A(af2=ANY, af4=None, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='nc') >>
     A(af2=ANY, af4=None, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='c'),
     k_a_fadnc_dicarb_af2_to_fadc)

# A-FADc-Dicarb-AF2 binds AF4
# NOTE: unbinding rule is defined later

Parameter('kf_a_fadc_dicarb_af2_binds_af4', 1)

Rule('A_FADc_Dicarb_AF2_binds_AF4',
     A(af2=ANY, af4=None, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='c') + AF4(a=None) >>
     A(af2=ANY, af4=2, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='c') % AF4(a=2),
     kf_a_fadc_dicarb_af2_binds_af4)

##########

# A-FADnc-Dicarb-AF4 binding + state change rule

Parameter('k_a_fadnc_dicarb_af4_binds_af2', 1)

Rule('A_FADnc_Dicarb_AF4_binds_AF2',
     A(af2=None, af4=ANY, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='nc') + AF2(a=None) >>
     A(af2=2, af4=ANY, fad=1, dicarb=ANY, bcd=None) % FAD(a=1, state='c') % AF2(a=2),
     k_a_fadnc_dicarb_af4_binds_af2)

# A-FADnc-AF2-AF4 binding + state change rule

Parameter('k_a_fadnc_af2_af4_binds_dicarb', 1)

Rule('A_FADnc_AF2_AF4_binds_Dicarb',
     A(af2=ANY, af4=ANY, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + Dicarb(a=None) >>
     A(af2=ANY, af4=ANY, fad=1, dicarb=2, bcd=None) % FAD(a=1, state='c') % Dicarb(a=2),
     k_a_fadnc_af2_af4_binds_dicarb)

##########

# A-FADc unbinding rules

Parameter('k_a_fadc_dicarb_unbinds_dicarb', 1)
Parameter('k_a_fadc_af2_unbinds_af2', 1)
Parameter('k_a_fadc_af4_unbinds_af4', 1)

# NOTE: Sequence of unbinding events is currently unknown
# Currently modeling dissociation events as independent - May change this in the future

Rule('A_FADc_Dicarb_unbinds_Dicarb',
     A(fad=1, dicarb=2, bcd=None) % FAD(a=1, state='c') % Dicarb(a=2) >>
     A(fad=1, dicarb=None, bcd=None) % FAD(a=1, state='c') + Dicarb(a=None),
     k_a_fadc_dicarb_unbinds_dicarb)

Rule('A_FADc_AF2_unbinds_AF2',
     A(af2=2, fad=1, bcd=None) % FAD(a=1, state='c') % AF2(a=2) >>
     A(af2=None, fad=1, bcd=None) % FAD(a=1, state='c') + AF2(a=None),
     k_a_fadc_af2_unbinds_af2)

Rule('A_FADc_AF4_unbinds_AF4',
     A(af4=2, fad=1, bcd=None) % FAD(a=1, state='c') % AF4(a=2) >>
     A(af4=None, fad=1, bcd=None) % FAD(a=1, state='c') + AF4(a=None),
     k_a_fadc_af4_unbinds_af4)

##########

# Formation of partially active CII (w/ non-covalent FAD)

Parameter('k_a_fadnc_af2_binds_bcd', 1)
Parameter('k_a_fadnc_af4_binds_bcd', 1)

Rule('A_FADnc_binds_BCD',
     A(af2=None, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='nc') + BCD(a=None) >>
     A(af2=None, af4=None, fad=1, dicarb=None, bcd=2) % FAD(a=1, state='nc') % BCD(a=2),
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
     A(af2=None, af4=None, fad=1, dicarb=None, bcd=None) % FAD(a=1, state='c') + BCD(a=None) >>
     A(af2=None, af4=None, fad=1, dicarb=None, bcd=2) % FAD(a=1, state='c') % BCD(a=2),
     k_a_fadc_binds_bcd)

Observable('Free_SDHA', A(af2=None, af4=None, fad=None, dicarb=None, bcd=None))
Observable('SDHA_without_FAD', A(fad=None, bcd=None))
Observable('SDHA_FADnc', A(fad=1, bcd=None) % FAD(a=1, state='nc'))
Observable('SDHA_FADc', A(fad=1, bcd=None) % FAD(a=1, state='c'))
Observable('CII_with_FADnc', A(fad=1, bcd=ANY) % FAD(a=1, state='nc'))
Observable('Active_CII', A(fad=1, bcd=ANY) % FAD(a=1, state='c'))

print(len(model.parameters))
print(model.parameters)

# simulation commands

tspan = np.linspace(0, 10, 101)
sim = ScipyOdeSimulator(model, tspan, verbose=True)
out = sim.run()

# print()
# print(len(model.species))
# quit()

# for i, sp in enumerate(model.species):
#     print(i, sp)
# print()
# for i, rxn in enumerate(model.reactions):
#     print(i, rxn)

for obs in model.observables:
    plt.plot(tspan, out.observables[obs.name], lw=2, label=obs.name)
plt.xlabel('time')
plt.ylabel('amount')
plt.legend(loc='best')

plt.tight_layout()
plt.show()
