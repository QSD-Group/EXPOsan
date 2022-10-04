# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

'''
TODO:
    - Add cost
    - Some of the design may refer to Table 10-3 on page 1068 of Metcalf and Eddy
        - Process (a) low loaded anaerobic lagoon system or (n) plug flow anaerobic system
'''

import qsdsan as qs
from flexsolve import IQ_interpolation
from qsdsan import WasteStream, sanunits as su
from exposan.cas import create_components

__all__ = ('create_system')

def create_system_without_cmps(flowsheet=None):
    flowsheet = flowsheet or qs.Flowsheet('cas')
    qs.main_flowsheet.set_flowsheet(flowsheet)
    
    # Influent
    inf = WasteStream('inf')
    inf.set_flow_by_concentration(
        flow_tot=20,
        concentrations={'Substrate': 300, 'X_inert': 100,},
        units=('mgd', 'mg/L'))
    
    # Units and the system
    with qs.System('sys') as sys:
        U1 = su.Screening('U1', ins=inf)
        M1 = su.Mixer('M1', ins=(U1-0, ''))
        ASP = su.ActivatedSludgeProcess('ASP', ins=(M1-0, 'ASP_air'),
                                    outs=('treated', 'was', 'offgas'))
        GBT = su.BeltThickener('GBT', ins=ASP-1, outs=(1-M1, 'thickened'))
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
        GBT.add_specification(GBT_spec)
        
        AD = su.SludgeDigester('AD', ins=GBT-1, outs=('disposed', 'biogas'))
        
        su.CHP('CHP', ins=(AD-1, 'natural_gas', 'CHP_air'), outs=('emission', 'solids'))
    
    return sys


def create_system(flowsheet=None):
    try: sys = create_system_without_cmps(flowsheet=flowsheet)
    except: # mostly for testing where components from the previous system will be carried over
        create_components()
        sys = create_system_without_cmps(flowsheet=flowsheet)
    return sys


# %%

if __name__ == '__main__':
    sys = create_system()
    sys.simulate()