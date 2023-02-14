# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from qsdsan import (
    processes as pc, 
    WasteStream, System, TEA, LCA, PowerUtility,
    ImpactIndicator as IInd, 
    ImpactItem as IItm, 
    StreamImpactItem as SIItm,
    )
from qsdsan.utils import ospath, ExogenousDynamicVariable as EDV
from exposan.metab_mock import (
    rhos_adm1_ph_ctrl,
    DegassingMembrane as DM, 
    METAB_AnCSTR as AB,
    IronSpongeTreatment as IST,
    DoubleMembraneGasHolder as GH,
    data_path
    )

folder = ospath.dirname(__file__)

__all__ = (
    'create_systems', 
    'default_inf_concs',
    'default_R1_init_conds',
    'default_R2_init_conds',
    'default_R_init_conds',
    'R1_ss_conds',
    'R2_ss_conds',
    'yields_bl', 'mus_bl', 'Ks_bl',
    'fermenters', 'methanogens', 'biomass_IDs',
    'vfa_IDs',
    )

#%% default values
scale = 1
Q = 5*scale          # influent flowrate [m3/d]
T1 = 273.15+35  # temperature [K]
# T1 = 273.15+25
Vl1 = 5*scale         # liquid volume [m^3]
Vg1 = 0.556*(scale**0.5)     # headspace volume [m^3]
ph1 = 5.8

T2 = 273.15+25    
Vl2 = 75*scale
Vg2 = 5*(scale**0.5)
ph2 = 7.2

T3 = T2
Vl3 = 5*scale
Vg3 = 0.556*(scale**0.5)

bl = 1   # yr, bead lifetime
# bl = 10

Temp1 = EDV('T1', function=lambda t: T1)
Temp2 = EDV('T2', function=lambda t: T2)
pH1 = EDV('pH1', function=lambda t: ph1)
pH2 = EDV('pH2', function=lambda t: ph2)

fermenters = ('X_su', 'X_aa', 'X_fa', 'X_c4', 'X_pro')
methanogens = ('X_ac', 'X_h2')
biomass_IDs = (*fermenters, *methanogens)
vfa_IDs = ('S_va', 'S_bu', 'S_pro', 'S_ac')

C_mw = 12.0107
N_mw = 14.0067

default_inf_concs = {
    'S_su':3.0,
    'S_aa':0.6,
    'S_fa':0.4,
    'S_va':0.4,
    'S_bu':0.4,
    'S_pro':0.4,
    'S_ac':0.4,
    'S_h2':5e-9,
    'S_ch4':5e-6,
    'S_IC':0.04*C_mw,
    'S_IN':0.01*N_mw,
    'S_I':0.02,
    'X_c':0.1,
    'X_ch':0.3,
    'X_pr':0.5,
    'X_li':0.25,
    'X_aa':1e-3,
    'X_fa':1e-3,
    'X_c4':1e-3,
    'X_pro':1e-3, 
    'X_ac':1e-3, 
    'X_h2':1e-3, 
    'X_I':0.025, 
    'S_cat':0.04, 
    'S_an':0.02
    }

yields_bl = {
         'Y_su': 0.1,
         'Y_aa': 0.08,
         'Y_fa': 0.06,
         'Y_c4': 0.06,
         'Y_pro': 0.04,
         'Y_ac': 0.05,
         'Y_h2': 0.06
         }

mus_bl = np.array([5.0e-01, 1.0e+01, 1.0e+01, 1.0e+01, 3.0e+01, 5.0e+01, 6.0e+00,
                   2.0e+01, 2.0e+01, 1.3e+01, 8.0e+00, 3.5e+01, 2.0e-02, 2.0e-02,
                   2.0e-02, 2.0e-02, 2.0e-02, 2.0e-02, 2.0e-02])

Ks_bl = np.array([5.0e-01, 3.0e-01, 4.0e-01, 2.0e-01, 
                  2.0e-01, 1.0e-01, 1.5e-01, 7.0e-06])

default_R1_init_conds = {
    'S_su': 0.0124*1e3,
    'S_aa': 0.0055*1e3,
    'S_fa': 0.1074*1e3,
    'S_va': 0.0123*1e3,
    'S_bu': 0.0140*1e3,
    'S_pro': 0.0176*1e3,
    'S_ac': 0.0893*1e3,
    'S_h2': 2.5055e-7*1e3,
    'S_ch4': 0.0555*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.1309*1e3,
    'X_ch': 0.0205*1e3,
    'X_pr': 0.0842*1e3,
    'X_li': 0.0436*1e3,
    'X_su': 1.87*1e3,
    'X_aa': 5.58*1e3,
    'X_fa': 2.03*1e3,
    'X_c4': 2.15*1e3,
    'X_pro': 1.00*1e3,
    }

default_R2_init_conds = {
    'S_su': 0.0124*1e3,
    'S_aa': 0.0055*1e3,
    'S_fa': 0.1074*1e3,
    'S_va': 0.0123*1e3,
    'S_bu': 0.0140*1e3,
    'S_pro': 0.0176*1e3,
    'S_ac': 0.0893*1e3,
    'S_h2': 2.5055e-7*1e3,
    'S_ch4': 0.0555*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.1309*1e3,
    'X_ac': 8.80*1e3,
    'X_h2': 3.70*1e3,
    }

default_R_init_conds = {
 'S_su': 0.32679646003805314,
 'S_aa': 0.3801610819236484,
 'S_fa': 12.437222319633748,
 'S_va': 0.3719673543175543,
 'S_bu': 0.47283246583627593,
 'S_pro': 0.3946420365926535,
 'S_ac': 10.182894473261367,
 'S_h2': 1.1655700001506622e-05,
 'S_ch4': 67.17348627201263,
 'S_IC': 846.4879614661522,
 'S_IN': 222.13725282096297,
 'S_I': 110.71467278942289,
 'X_c': 107.43132381172228,
 'X_ch': 1.2600235711799973,
 'X_pr': 1.3804329631122664,
 'X_li': 1.7696259648387357,
 'X_su': 732.9760678333023,
 'X_aa': 224.81751931525334,
 'X_fa': 126.7301174776879,
 'X_c4': 227.8726398428066,
 'X_pro': 140.2738127019708,
 'X_ac': 669.4626559278454,
 'X_h2': 245.67774602566578,
 'X_I': 206.42934561053158,
 'S_cat': 40.0,
 'S_an': 20.0,
 }

R1_ss_conds = {
    'S_su': 0.0145871088552909*1e3,
    'S_aa': 0.00643308564144693*1e3,
    'S_fa': 0.634823005990967*1e3,
    'S_va': 0.624510322247682*1e3,
    'S_bu': 1.03793927591996*1e3,
    'S_pro': 1.24676871525373*1e3,
    'S_ac': 2.00250371674824*1e3,
    'S_h2': 0.00850364943532684*1e3,
    'S_ch4': 0.0000422133982597226*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.027310256066728*1e3,
    'X_c': 0.146203507058736*1e3,
    'X_ch': 0.0286018513139117*1e3,
    'X_pr': 0.0467836694957302*1e3,
    'X_li': 0.0247209587890493*1e3,
    'X_su': 4.69052782535406*1e3,
    'X_aa': 1.22829926704024*1e3,
    'X_fa': 0.0147446263753011*1e3,
    'X_c4': 0.0149933579422897*1e3,
    'X_pro': 0.0145343147735253*1e3,
    'X_ac': 0.00098041337766024*1e3,
    'X_h2': 0.00110808891184369*1e3,
    'X_I': 0.0396205121367899*1e3
    }

R2_ss_conds = {
    'S_su': 0.00106990968535691*1e3,
    'S_aa': 0.00125571416517827*1e3,
    'S_fa': 0.121097573221394*1e3,
    'S_va': 0.0132519103137696*1e3,
    'S_bu': 0.0172912281196732*1e3,
    'S_pro': 0.020032163173878*1e3,
    'S_ac': 0.00574002755366853*1e3,
    'S_h2': 3.76969944940856e-08*1e3,
    'S_ch4': 0.0499411746585487*1e3,
    'S_IC': 0.0951*C_mw*1e3,
    'S_IN': 0.0945*N_mw*1e3,
    'S_I': 0.105601391746794*1e3,
    'X_c': 0.0897520281015078*1e3,
    'X_ch': 0.00108163641708242*1e3,
    'X_pr': 0.00120204580901502*1e3,
    'X_li': 0.00150204523369107*1e3,
    'X_su': 0.195961987850137*1e3,
    'X_aa': 0.059723477130333*1e3,
    'X_fa': 0.0351858744892462*1e3,
    'X_c4': 0.0812315951844566*1e3,
    'X_pro': 0.0503466475437059*1e3,
    'X_ac': 1.1653549028287*1e3,
    'X_h2': 0.4352809013846*1e3,
    'X_I': 0.196117291164614*1e3
    }

#%% Systems
IInd.load_from_file(ospath.join(data_path, 'TRACI_indicators.xlsx'), sheet=0)
IItm.load_from_file(ospath.join(data_path, '_impact_items.xlsx'))
IItm('Stainless_steel', source='stainless_steel')
bg_offset_CFs = IItm.get_item('biogas_offset').CFs
NaOCl_item = IItm.get_item('NaOCl')
citric_acid_item = IItm.get_item('citric_acid')

def get_fug_ch4(ws):
    '''Returns kg/hr fugitive CH4.'''
    cmps = ws.components
    return ws.imass['S_ch4']*cmps.S_ch4.i_mass

def get_NG_eq(bg):
    '''Returns m3/hr natural gas equivalent.'''
    cmps = bg.components
    KJ_per_kg = cmps.i_mass/cmps.chem_MW*cmps.LHV
    return -sum(bg.mass*KJ_per_kg)*1e-3/39

def add_strm_iitm(sys):
    for ws in sys.products:
        if ws.phase == 'l':
            SIItm(ID=f'{ws.ID}_fugitive_ch4', linked_stream=ws, 
                  flow_getter=get_fug_ch4,
                  GWP100=28, MIR=0.0143794871794872)
        elif ws.phase == 'g':
            SIItm(ID=f'{ws.ID}_NG_offset', linked_stream=ws, functional_unit='m3',
                  flow_getter=get_NG_eq, # m3/hr natural-gas-equivalent
                  **bg_offset_CFs)
    for u in sys.units:
        if isinstance(u, DM):
            SIItm(linked_stream=u.NaOCl, **NaOCl_item.CFs)
            SIItm(linked_stream=u.citric_acid, **citric_acid_item.CFs)

# Operation items
kWh = lambda lca: lca.system.power_utility.rate*lca.lifetime_hr
def MJ(lca):
    sys = lca.system
    duties = [hu.duty for hu in sys.heat_utilities]
    MJ_per_hr = sum(d for d in duties if d > 0)*1e-3
    return MJ_per_hr*lca.lifetime_hr
    
def add_TEA_LCA(sys, irr, lt):
    TEA(sys, discount_rate=irr, lifetime=lt, simulate_system=False, CEPCI=708)   
    LCA(
        sys, lifetime=lt, simulate_system=False,
        electricity = kWh,
        heat_onsite = MJ
        )

def create_systems(lifetime=30, discount_rate=0.1, electric_price=0.0913,
                   flowsheet_A=None, flowsheet_B=None, flowsheet_C=None, flowsheet_D=None,
                   inf_concs={}, R1_init_conds={}, R2_init_conds={}, R_init_conds={}, 
                   which=None, selective=False):
    PowerUtility.price = electric_price
    which = which or ('A', 'B', 'C', 'D')
    if isinstance(which, str): which = (which,)
    
    if selective: 
        R1_retain = fermenters
        R2_retain = methanogens
    else: R1_retain = R2_retain = biomass_IDs
    pc.create_adm1_cmps()
    adm1 = pc.ADM1()
    adm1_phctrl = pc.ADM1()
    dyn_params = adm1_phctrl.rate_function.params.copy()
    adm1_phctrl.set_rate_function(rhos_adm1_ph_ctrl)
    adm1_phctrl.rate_function._params = dyn_params
    
    inf_concs = inf_concs or default_inf_concs.copy()
    brewery_ww = WasteStream('BreweryWW_A')
    brewery_ww.set_flow_by_concentration(Q, concentrations=inf_concs, 
                                         units=('m3/d', 'kg/m3'))
    R1_init_conds = R1_init_conds or R1_ss_conds
    R2_init_conds = R2_init_conds or R2_ss_conds
    R_init_conds = R_init_conds or default_R_init_conds
    systems = []
    
    if 'A' in which:
        flowsheet_A = flowsheet_A or qs.Flowsheet('A')
        qs.main_flowsheet.set_flowsheet(flowsheet_A)
        flowsheet_A.stream.register(brewery_ww.ID, brewery_ww)
        
        effA = WasteStream('Effluent_A', T=T2)
        bg1A = WasteStream('Biogas_1A', phase='g')
        bg2A = WasteStream('Biogas_2A', phase='g')
        
        ISTA = IST(ID='ISTA')
        GHA = GH(ID='GHA')
        
        R1A = AB('R1A', ins=brewery_ww, outs=(bg1A, ''), 
                 V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1_phctrl, 
                 retain_cmps=R1_retain, bead_lifetime=bl,
                 exogenous_vars=(Temp1, pH1),
                 F_BM_default=1)
        R2A = AB('R2A', ins=R1A-1, outs=(bg2A, effA), 
                 V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1_phctrl,
                 retain_cmps=R2_retain, bead_lifetime=bl,
                 exogenous_vars=(Temp2, pH2),
                 equipment=[ISTA, GHA],
                 F_BM_default=1)
    
        R1A.set_init_conc(**R1_init_conds)
        R2A.set_init_conc(**R2_init_conds)
    
        sysA = System('sysA', path=(R1A, R2A))
        sysA.set_dynamic_tracker(R1A, R2A, bg1A, bg2A, effA)
        add_strm_iitm(sysA)
        add_TEA_LCA(sysA, discount_rate, lifetime)
        systems.append(sysA)
        
    if 'B' in which:
        flowsheet_B = flowsheet_B or qs.Flowsheet('B')
        qs.main_flowsheet.set_flowsheet(flowsheet_B)
        
        infB = brewery_ww.copy('BreweryWW_B')
        effB = WasteStream('Effluent_B', T=T2)
        bg1B = WasteStream('Biogas_1B', phase='g')
        bgh2B = WasteStream('Biogas_hsp_2B', phase='g')
        bgm2B = WasteStream('Biogas_mem_2B', phase='g')
        
        ISTB = IST(ID='ISTB')
        GHB = GH(ID='GHB')
        
        R1B = AB('R1B', ins=infB, outs=(bg1B, ''), 
                 V_liq=Vl1, V_gas=Vg1, T=T1, model=adm1_phctrl,
                 retain_cmps=R1_retain, bead_lifetime=bl,
                 exogenous_vars=(Temp1, pH1),
                 F_BM_default=1)
        R2B = AB('R2B', ins=R1B-1, outs=(bgh2B, ''), 
                 V_liq=Vl2, V_gas=Vg2, T=T2, model=adm1_phctrl,
                 retain_cmps=R2_retain, bead_lifetime=bl,
                 exogenous_vars=(Temp2, pH2),
                 equipment=[ISTB, GHB],
                 F_BM_default=1)
        DM2B = DM('DM2B', ins=R2B-1, outs=(bgm2B, effB),
                  F_BM_default=1)
    
        R1B.set_init_conc(**R1_init_conds)
        R2B.set_init_conc(**R2_init_conds)
        
        sysB = System('sysB', path=(R1B, R2B, DM2B))
        sysB.set_dynamic_tracker(R1B, R2B, bg1B, bgh2B, bgm2B, effB)
        add_strm_iitm(sysB)
        add_TEA_LCA(sysB, discount_rate, lifetime)
        systems.append(sysB)
    
    if 'C' in which:
        flowsheet_C = flowsheet_C or qs.Flowsheet('C')
        qs.main_flowsheet.set_flowsheet(flowsheet_C)
        
        infC = brewery_ww.copy('BreweryWW_C')
        effC = WasteStream('Effluent_C', T=T2)
        bgC = WasteStream('Biogas_C', phase='g')

        ISTC = IST(ID='ISTC')
        GHC = GH(ID='GHC')
        
        RC = AB('RC', ins=infC, outs=(bgC, effC), 
                V_liq=Vl3, V_gas=Vg3, T=T3, model=adm1,
                retain_cmps=biomass_IDs, bead_lifetime=bl,
                equipment=[ISTC, GHC],
                F_BM_default=1)        
        RC.set_init_conc(**R_init_conds)
    
        sysC = System('sysC', path=(RC,))
        sysC.set_dynamic_tracker(RC, bgC, effC)
        add_strm_iitm(sysC)
        add_TEA_LCA(sysC, discount_rate, lifetime)
        systems.append(sysC)
        
    if 'D' in which:
        flowsheet_D = flowsheet_D or qs.Flowsheet('D')
        qs.main_flowsheet.set_flowsheet(flowsheet_D)
        
        infD = brewery_ww.copy('BreweryWW_D')
        effD = WasteStream('Effluent_D', T=T2)
        bghD = WasteStream('Biogas_hsp_D', phase='g')
        bgmD = WasteStream('Biogas_mem_D', phase='g')

        ISTD = IST(ID='ISTD')
        GHD = GH(ID='GHD')
        
        RD = AB('RD', ins=infD, outs=(bghD, ''), 
                V_liq=Vl3, V_gas=Vg3, T=T3, model=adm1,
                retain_cmps=biomass_IDs, bead_lifetime=bl,
                equipment=[ISTD, GHD],
                F_BM_default=1)
        DMD = DM('DMD', ins=RD-1, outs=(bgmD, effD),
                  F_BM_default=1)
        RD.set_init_conc(**R_init_conds)
    
        sysD = System('sysD', path=(RD, DMD))
        sysD.set_dynamic_tracker(RD, bghD, bgmD, effD)
        add_strm_iitm(sysD)
        add_TEA_LCA(sysD, discount_rate, lifetime)
        systems.append(sysD)
    
    return systems

#%%
if __name__ == '__main__':
    systems = create_systems()
    for sys in systems:
        sys.simulate(
            t_span=(0,200),
            state_reset_hook='reset_cache',
            method='BDF'
            )
        sys.diagram()