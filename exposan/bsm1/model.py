# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Joy Zhang <joycheung1994@gmail.com>
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter, AttrGetter
from exposan import bsm1 as bm


__all__ = ('model_bsm1',)


bsm1 = bm.bsm1
s = bm.system
model_bsm1 = qs.Model(system=bsm1)

########## Add Uncertainty Parameters ##########
param = model_bsm1.parameter
# Give ±10% of the baseline
get_uniform_w_frac = lambda b, frac: shape.Uniform(lower=b*(1-frac), upper=b*(1+frac))

# Flow rates
# b = s.Q
# D = get_uniform_w_frac(b, 0.1)
# PE = s.PE
# inf_kwargs = s.inf_kwargs
# @param(name='Influent flowrate', element=PE, kind='coupled', units='m3/d',
#        baseline=b, distribution=D)
# def set_Q_inf(i):
#     PE.set_flow_by_concentration(i, **inf_kwargs)

b = s.Q_was
D = get_uniform_w_frac(b, 0.1)
C1 = s.C1
@param(name='Waste sludge flowrate', element=C1, kind='coupled', units='m3/d',
       baseline=b, distribution=D)
def set_Q_was(i):
    C1.wastage = i

# b = s.Q_ras
# D = get_uniform_w_frac(b, 0.1)
# @param(name='Return sludge flowrate', element=C1, kind='coupled', units='m3/d',
#        baseline=b, distribution=D)
# def set_Q_ras(i):
#     C1.underflow = i

Q, V_an, V_ae = s.Q, s.V_an, s.V_ae
b = 1
D = shape.Uniform(lower=0.75, upper=1)
@param(name='Sludge recycling as a fraction of influent', element=C1, 
       kind='coupled', units='', baseline=b, distribution=D)
def set_Q_ras(i):
    C1.underflow = Q*i

b = 3
D = shape.Uniform(lower=2.25, upper=3.75)
O3 = s.O3
@param(name='Internal recirculation rate as a fraction of influent', element=O3, 
       kind='coupled', units='', baseline=b, distribution=D)
def set_Q_intr(i):
    intr = i/(1+i+C1._Qras/Q)
    O3.split = [intr, 1-intr]
    
b = V_an * 2 / Q * 24
D = shape.Uniform(lower=2.34, upper=b)
A1, A2 = s.A1, s.A2
@param(name='Anoxic zone hydraulic retention time', element=A1, 
       kind='coupled', units='hr', baseline=b, distribution=D)
def set_A1_A2_HRT(i):
    A1._V_max = A2._V_max = i / 24 * Q / 2

b = V_ae * 3 / Q * 24
D = shape.Uniform(lower=4.68, upper=b)
O1, O2 = s.O1, s.O2
@param(name='Aerobic zone hydraulic retention time', element=O1, 
       kind='coupled', units='hr', baseline=b, distribution=D)
def set_O1_O2_O3_HRT(i):
    O1._V_max = O2._V_max = O3._V_max = i / 24 * Q / 3
# Kinetic parameters based on
# Sin, G.; Gernaey, K. V.; Neumann, M. B.; van Loosdrecht, M. C. M.; Gujer, W.
# Uncertainty Analysis in WWTP Model Applications:
# A Critical Discussion Using an Example from Design.
# Water Research 2009, 43 (11), 2894–2906.
# https://doi.org/10.1016/j.watres.2009.03.048.
asm1 = s.asm1

param_ranges = {
    'mu_H': (4, 3, 5, '/d'), # default, min, max, unit
    'K_S': (10, 5, 15, 'g COD/m3'),
    'K_OH': (0.2, 0.1, 0.3, 'g O2/m3'),
    'K_NO': (0.5, 0.25, 0.75, 'g N/m3'),
    'b_H': (0.3, 0.285, 0.315, '/d'),
    'mu_A': (0.5, 0.475, 0.525, '/d'),
    'K_NH': (1, 0.5, 1.5, 'g N/m3'),
    'K_O_A': (0.4, 0.3, 0.5, 'g COD/m3'),
    'b_A': (0.05, 0.04, 0.06, '/d'),
    'eta_g': (0.8, 0.6, 1, ''),
    'k_a': (0.05, 0.03, 0.08, 'm3/g COD/d'),
    'k_h': (3, 2.25, 3.75, 'g X_S/g X_BH COD/d'),
    'K_X': (0.1, 0.075, 0.125, 'g X_S/g X_BH COD'),
    'eta_h': (0.8, 0.6, 1, ''),
    'Y_H': (0.67, 0.64, 0.7, 'g COD/g COD'),
    'Y_A': (0.24, 0.23, 0.25, 'g COD/g N'),
    'i_XB': (0.08, 0.04, 0.12, 'g N/g COD'),
    'i_XP': (0.06, 0.057, 0.063, 'g N/g COD'),
    # 'f_P': (0.75, 0.7, 0.95, 'g TSS/g COD'),
    'fr_SS_COD': (0.75, 0.7, 0.95, 'g TSS/g COD'),
    }

# Use a loop to add the parameters in batch
asm1_baseline_param = asm1.parameters
for name, vals in param_ranges.items():
    b, lb, ub, unit = vals
    set_val = asm1_baseline_param.get(name)
    if set_val != b:
        warn(f'The value for parameter `{name}` set in the system script is {set_val}, '
             f'different from the provided baseline of {b}, is this intentional?')
    D = shape.Triangle(lower=lb, midpoint=b, upper=ub)
    setter = DictAttrSetter(asm1, '_parameters', keys=(name,))
    param(setter=setter, name=name, element=A1, kind='coupled', units=unit,
          baseline=b, distribution=D)

# Set f_P by setting the since f_P is calculataed from f_Pobs and Y_H,
# f_P = f_Pobs*(1-Y_H)/(1-f_Pobs*Y_H)
# the calcualted baseline is 0.2*(1-0.67)/(1-0.2*0.67)=0.076,
# close enough to the set 0.08
b = 0.2
D = shape.Triangle(lower=0.15, midpoint=b, upper=0.25)
@param(name='f_Pobs', element=A1, kind='coupled', units='',
       baseline=b, distribution=D)
def set_f_Pobs(i):
    Y_H = asm1.Y_H
    asm1.f_P = i*(1-Y_H) / (1-i*Y_H)


# Aeration
aer1 = s.aer1
b = aer1.KLa
D = shape.Uniform(lower=180, upper=360)
@param(name='O1 and O2 KLa', element=O1, kind='coupled', units='',
       baseline=b, distribution=D)
def set_O1_O2_KLa(i):
    aer1.KLa = i

aer2 = s.aer2
b = aer2.KLa
D = get_uniform_w_frac(b, 0.1)
@param(name='O3 KLa', element=O3, kind='coupled', units='',
       baseline=b, distribution=D)
def set_O3_KLa(i):
    aer2.KLa = i
    
b = 8.0
D = get_uniform_w_frac(b, 0.1)
@param(name='Saturation DO', element=O1, kind='coupled', units='mg/L',
       baseline=b, distribution=D)
def set_DOsat(i):
    aer1.DOsat = i
    aer2.DOsat = i

########## Add Evaluation Metrics ##########
metric = model_bsm1.metric

# Use a function to batch-add metrics
# concentrations of all of the components for the clarifier effluent (seconday effluent),
# cmps = s.cmps
# def add_conc_as_metrics(ws):
#     for cmp in cmps:
#         getter = AttrGetter(ws, attr='iconc',
#                             hook=lambda iconc, ID: iconc[ID],
#                             hook_param=(cmp.ID,))
#         metric(getter=getter, name=f'{ws.ID} {cmp.ID}', 
#                units='mg/L', element='WasteStreams')

SE, RAS, WAS, RE = s.SE, s.RAS, s.WAS, s.RE
# for ws in (SE, RAS, WAS, RE):
#     add_conc_as_metrics(ws)

# Effluent composite variables and daily sludge production
for i in ('COD', 'TN', 'TKN'):
    metric(getter=AttrGetter(SE, attr=i), name=i, 
           units='mg/L', element='Effluent')
metric(getter=SE.get_TSS, name='TSS', units='mg/L', element='Effluent')

def get_daily_sludge_production():
    return WAS.get_TSS() * 1e-3 * WAS.F_vol * 24
metric(getter=get_daily_sludge_production, 
       name='Daily sludge production', units='kg TSS/d', element='WAS')