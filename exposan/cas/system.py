# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <zoe.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from flexsolve import IQ_interpolation
from qsdsan import Component, Components, WasteStream, set_thermo, System
from qsdsan.sanunits import Screening, Mixer, ActivatedSludgeProcess, \
    BeltThickener, SludgeDigester, CHP as CHPunit

__all__ = ('cmps', 'sys',)


# =============================================================================
# Components and streams
# =============================================================================

X = Component('X', phase='s', measured_as='COD', i_COD=0, description='Biomass',
              organic=True, particle_size='Particulate', degradability='Readily')
X_inert = Component('X_inert', phase='s', description='Inert biomass', i_COD=0,
                    organic=True, particle_size='Particulate', degradability='Undegradable')
Substrate = Component('Substrate', phase='l', i_COD=1,
                      organic=True, particle_size='Soluble', degradability='Readily')
CH4 = Component('CH4', phase='g', organic=True, particle_size='Dissolved gas',
                degradability='Readily')
cmps = Components([X, X_inert, Substrate, CH4])

# Default-adding components needed for combustion
cmps = Components.append_combustion_components(cmps, lock_state_at='')

# Define component groups
cmps.define_group('active_biomass', IDs=('X',))
cmps.define_group('inert_biomass', IDs=('X_inert',))
cmps.define_group('substrates', IDs=('Substrate',))

set_thermo(cmps)


inf = WasteStream('inf')
inf.set_flow_by_concentration(flow_tot=20,
                              concentrations={
                                  'Substrate': 300,
                                  'X_inert': 100,
                                  },
                              units=('mgd', 'mg/L'))


# =============================================================================
# Units and the system
# =============================================================================

U1 = Screening('U1', ins=inf)
M1 = Mixer('M1', ins=(U1-0, ''))
ASP = ActivatedSludgeProcess('ASP', ins=(M1-0, 'ASP_air'),
                            outs=('treated', 'was', 'offgas'))
GBT = BeltThickener('GBT', ins=ASP-1, outs=(1-M1, 'thickened'))
thickened = GBT.outs[1]
# Update the biomass concentration to be at the set value
# There will be a warning when solving it, due to the unrealistic
# sludge moisture content (i.e., it exceeds in the moisture content of the feeds)
# but it is fine since that's only testing the bound
def X_inert_at_mc(mc):
    GBT.sludge_moisture = mc
    GBT._set_split_at_mc()
    return thickened.iconc['X_inert']-25000
def GBT_spec():
    lb = 1e-3
    mixed_F_mass = GBT._mixed.F_mass
    if mixed_F_mass == 0:
        ub = 1-1e-3
    else:
        ub = max(lb+1e-3, GBT._mixed.imass['Water']/mixed_F_mass)
    IQ_interpolation(f=X_inert_at_mc, x0=lb, x1=ub, xtol=1e-3, ytol=1,
                     checkbounds=False)
GBT.specification = GBT_spec

AD = SludgeDigester('AD', ins=GBT-1, outs=('disposed', 'biogas'))

CHP = CHPunit('CHP', ins=(AD-1, 'natural_gas', 'CHP_air'), outs=('emission', 'solids'))

sys = System('sys', path=(U1, M1, ASP, GBT, AD, CHP))
sys.simulate()