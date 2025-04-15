# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
        
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

'''
import numpy as np
from qsdsan.sanunits import CSTR, PFR, CompletelyMixedMBR, AerobicDigester
from qsdsan.processes import DiffusedAeration
from warnings import warn

__all__ = ('get_aeration_kLa', 
           'get_crossflow_kLa',
           'get_aeration_demand',
           'aer_kwargs',
           'plantwide_aeration_demand',
           'plantwide_aeration_energy',
           )

isa = isinstance
def get_aeration_kLa(unit, DO_sat):
    if isa(unit, CSTR):
        DO = unit._aeration
        if DO is None: return 0
        if 0 < DO < DO_sat : return unit._cache_OTR/(DO_sat - DO)
        elif DO >= DO_sat: return np.inf
        else: return 0
    elif isa(unit, PFR):
        kLas = []
        for DO, otr in zip(unit.DO_setpoints, unit._cache_OTR):
            if 0 < DO < DO_sat : kLas.append(otr/(DO_sat - DO))
            elif DO >= DO_sat: kLas.append(np.inf)
            else: kLas.append(0)
        return np.array(kLas)
    else:
        raise TypeError(f'unrecognized reactor type {type(unit)}')

aer_kwargs = {
    'fine_bubble': dict(
        beta=0.95, 
        theta=1.024,
        F=0.6,
        SOTE=0.28, 
        fine_pore=True
        ),
    'cross_flow': dict(
        alpha=0.65, 
        beta=0.95, 
        theta=1.024,
        SOTE=0.1,
        fine_pore=False
        ),
    'digester': dict(
        beta=0.95, 
        theta=1.024,
        F=1.0,
        SOTE=0.1425,
        fine_pore=False,
        )
    }

ft2m = 0.3048
cfm2cmd = 40.77625909247999

def get_crossflow_kLa(unit, Q_air, DOsat_s20=8.0, T_air=13,):
    '''
    Estimate kLa [d^(-1)] for cross-flow air in membrane bioreactor.

    Parameters
    ----------
    unit : :class:`CompletelyMixedMBR`
        Membrane bioreactor unit object.
    Q_air : float
        Cross flow air flowrate, in ft3/min.
    T_air : float, optional
        Inlet air temperature, in degree Celcius.

    See Also
    --------
    :class:`qsdsan.processes.DiffusedAeration`
    
    '''
    if not isa(unit, CompletelyMixedMBR):
        warn(f'`get_crossflow_kLa` may not be applicable to {type(unit)}')
    Q_air *= cfm2cmd
    cfa = DiffusedAeration('cfa', DO_ID=unit.DO_ID, V=unit.V_max, Q_air=Q_air, 
                           DOsat_s20=DOsat_s20, d_submergence=4.3, 
                           T_water=unit.outs[0].T, T_air=T_air+273.15, 
                           **aer_kwargs['cross_flow'])
    return cfa.kLa
    

def _setup_diffused_aeration(unit, kLa, DOsat_s20, T_air,
                             tank_depth, diffuser_height, 
                             aeration_type, **kwargs):
    # breakpoint()
    if kLa is None: 
        kLa = get_aeration_kLa(unit, DOsat_s20)
        if isa(unit, CompletelyMixedMBR):
            Q_cfa = kwargs.pop('Q_cfa', 20000)
            k_cfa = get_crossflow_kLa(unit, Q_cfa, DOsat_s20, T_air)
            if k_cfa > kLa: kLa = 0.
            else: kLa -= k_cfa
            
    d_submergence = (tank_depth - diffuser_height) * ft2m
    if aeration_type in aer_kwargs:
        kwargs.update(aer_kwargs[aeration_type])
    aer = DiffusedAeration('aer', DO_ID=unit.DO_ID, V=1000, KLa=120, DOsat_s20=DOsat_s20,
                           T_water=unit.outs[0].T, T_air=T_air+273.15,
                           d_submergence=d_submergence, **kwargs)
    return kLa, aer

def get_aeration_demand(unit, kLa=None, DOsat_s20=8.0, alpha=0.5, 
                        tank_depth=15, diffuser_height=1.0, T_air=13,
                        aeration_type='fine_bubble', **kwargs):
    '''
    Estimate aeration air flowrate [m3/d] from kLa.

    Parameters
    ----------
    unit : :class:`SanUnit`
        A sanunit with aeration demand.
    alpha : float, optional
        Alpha factor. Must have the same shape as kLa.
    tank_depth : float, optional
        Water depth in ft. The default is 15.
    diffuser_height : float, optional
        Height of the diffuser from the tank floor, in ft. The default is 1.0.
    T_air : float, optional
        Inlet air temperature, in degree Celcius. The default is 13.
    aeration_type : str, optional
        Provides default parameter values for certain types. 
        'fine_bubble' for suspended growth bioreactors, 
        'digester' for aerobic sludge digester.
    **kwargs : optional
        Any keyword arguments to pass onto the `DiffusedAeration` object.
    
    See Also
    --------
    :class:`qsdsan.processes.DiffusedAeration`

    '''
    
    kLa, aer = _setup_diffused_aeration(unit, kLa, DOsat_s20, T_air, tank_depth, 
                                        diffuser_height, aeration_type, **kwargs)
    
    if isa(unit, CSTR):
        aer.alpha = alpha
        aer.V = unit.V_max
        aer.kLa = kLa
        return aer.Q_air
    elif isa(unit, PFR):
        Q_air = 0
        for k, a, V in zip(kLa, alpha, unit.V_tanks):
            if k <= 0: continue
            aer.alpha = a
            aer.V = V
            aer.kLa = k
            Q_air += aer.Q_air
        return Q_air
    else:
        raise TypeError(f'unrecognized reactor type {type(unit)}')

def get_aeration_energy(Q_air, tank_depth=15, diffuser_height=1.0, T_air=13,
                        blowermotor_efficiency=0.7, P_inlet_loss=1.0,
                        P_outlet_loss=7.0, P_atm=101.325,
                        **kwargs):
    '''
    Estimate energy required for aeration [kW].

    Parameters
    ----------
    Q_air : float, optional
        Air flow rate [m3/d].
    blowermotor_efficiency : float, optional
        Combined blower and motor efficiency. The default is 0.7.
    P_inlet_loss : float, optional
        Pressure drop at the inlet [kPa]. The default is 1.0.
    P_outlet_loss : float, optional
        Head loss in piping and diffuser [kPa]. The default is 7.0.
    P_atm : float, optional
        Atmospheric pressure [kPa]. The default is 101.325.

    See Also
    --------
    `get_aeration_demand`

    '''
    d_submergence = (tank_depth - diffuser_height) * ft2m
    P_in = P_atm - P_inlet_loss
    P_out = P_atm + 9.81 * d_submergence + P_outlet_loss
    P_blower = 1.4161e-5 * (T_air+273.15) * Q_air * ((P_out/P_in)**0.283 - 1) / blowermotor_efficiency
    return P_blower
    
default_alpha = {
    'rBOD': [0.3, 0.3, 0.4, 0.4, 0.5, 0.5],
    'others': [0.6, 0.6, 0.7, 0.7, 0.8, 0.8],
    }

def plantwide_aeration_demand(system):
    if system.ID[0] in 'BC': asr_alpha = default_alpha['rBOD']
    else: asr_alpha = default_alpha['others']
    Q_air = {}
    for u in system.units:
        if isa(u, AerobicDigester):
            kwargs = dict(alpha=0.6, tank_depth=20, diffuser_height=1.0, 
                          aeration_type='digester')
        elif isa(u, PFR):
            kwargs = dict(alpha=asr_alpha[:u.N_tanks_in_series],
                          aeration_type='fine_bubble')
        elif isa(u, CompletelyMixedMBR):
            kwargs = dict(alpha=asr_alpha[-1],
                          aeration_type='fine_bubble')
            Q_air[u.ID+'_cfa'] = 20000 * cfm2cmd
        elif isa(u, CSTR) and u.ID[-1].isdigit():
            i = int(u.ID[-1]) - 1
            kwargs = dict(alpha=asr_alpha[i],
                          aeration_type='fine_bubble')
        else:
            continue
        Q_air[u.ID] = get_aeration_demand(u, **kwargs)

    return Q_air

def plantwide_aeration_energy(system, Q_air_dct=None):
    if Q_air_dct is None:
        Q_air_dct = plantwide_aeration_demand(system)
    Pb = {}
    for uid, Q_air in Q_air_dct.items():
        kwargs = dict(P_inlet_loss=1.724, P_outlet_loss=17.24, blowermotor_efficiency=0.7)
        if uid == 'AED':
            kwargs.update(dict(
                tank_depth=20, diffuser_height=1.0, 
                blowermotor_efficiency=0.65
                ))
            if system.ID in ('E2P', 'N2'):
                kwargs['P_outlet_loss'] = 10.34
        Pb[uid] = get_aeration_energy(Q_air, **kwargs)
    return Pb