# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

Reference:
    Alex, J.; Benedetti, L.; Copp, J. B.; Gernaey, K. V.; Jeppsson, U.;
    Nopens, I.; Pons, M. N.; Rosen, C.; Steyer, J. P.; Vanrolleghem, P. A.
    Benchmark Simulation Model No. 2 (BSM2).
    http://iwa-mia.org/wp-content/uploads/2022/09/TR3_BSM_TG_Tech_Report_no_3_BSM2_General_Description.pdf.

'''
import numpy as np, qsdsan as qs
from qsdsan import WasteStream, sanunits as su, processes as pc

# Parameters for ASM No. 1 at 15 degC, Tables 2-3
default_asm_kwargs = dict(
    # Table 2 Stoichiometric parameters
    Y_H=0.67,
    Y_A=0.24,
    f_P=0.08,
    i_XB=0.08, #0.086;
    i_XP=0.06,
    fr_SS_COD=0.75, # X_I2TSS, X_S2TSS, X_BH2TSS, X_BA2TSS, X_P2TSS
    # Table 3 Kinetic parameters
    mu_H=4.0, #6.0;
    K_S=10.0, #20;
    K_O_H=0.2, # K_OH = 0.2;
    K_NO=0.5,
    b_H=0.3, #0.62;
    mu_A=0.5, #0.8;
    K_NH=1.0,
    K_O_A=0.4, # K_OA = 0.4;
    b_A=0.05, #0.2;
    eta_g=0.8, # ny_g
    k_a=0.05, #0.08;
    k_h=3.0,
    K_X=0.1, #0.03;
    eta_h=0.8, # ny_h #0.4;
    # path=os.path.join(data_path, '_asm1.tsv'),
    )

cmps_asm1 = pc.create_asm1_cmps()
thermo_asm1 = qs.get_thermo()

cmps_adm1 = pc.create_adm1_cmps()
thermo_adm1 = qs.get_thermo()
adm1 = pc.ADM1()
effluent = WasteStream('effluent', T=35+273.15)


# %%

qs.set_thermo(thermo_asm1)

# Anaerobic digester influent (pre ASM2ADM interface)
T_inf = 14.8581 + 273.15
default_inf_kwargs = {
    'flow_tot': 178.4674,
    'concentrations': {
        'S_I': 28.0665,
        'S_S': 48.9526,
        'X_I': 10361.7101,
        'X_S': 20375.0176,
        'X_BH': 10210.0698,
        'X_BA': 553.2808,
        'X_P': 3204.6601,
        'S_O': 0.25225,
        'S_NO': 1.6871,
        'S_NH': 28.9098,
        'S_ND': 4.6834,
        'X_ND': 906.0933,
        'S_ALK': 7.1549*12,
        },
    'units': ('m3/d', 'mg/L'),
    }
# TSS = 33528.5538 mg SS/l


influent = WasteStream('influent', T=14.8581+273.15)

influent.set_flow_by_concentration(**default_inf_kwargs)

J1 = su.ASMtoADM('J1', upstream=influent, downstream=effluent, 
                 thermo=thermo_adm1, isdynamic=True, adm1_model=adm1)

sys = qs.System('sys', path=(J1,))

t = 1
t_step = 1
# method = 'RK45'
# method = 'RK23'
# method = 'DOP853'
# method = 'Radau'
method = 'BDF'
# method = 'LSODA'

sys.simulate(
    state_reset_hook='reset_cache',
    t_span=(0,t),
    t_eval=np.arange(0, t+t_step, t_step),
    method=method,
    # rtol=1e-2,
    # atol=1e-3,
    # export_state_to=f'results/sol_{t}d_{method}.xlsx',
    )

effluent_conc = dict(zip(effluent.components.IDs, effluent.iconc.data))

# effluent_conc
# S_su 0.0
# S_aa 40.99187113202961
# S_fa 0.0
# S_va 0.0
# S_bu 0.0
# S_pro 0.0
# S_ac 0.0
# S_h2 0.0
# S_ch4 0.0
# S_IC 94.847730527175
# S_IN 28.711513621117057
# S_I 28.066499999999962
# X_c 0.0
# X_ch 3340.077084921927
# X_pr 16560.50572492688
# X_li 7793.513198151161
# X_su 0.0
# X_aa 0.0
# X_fa 0.0
# X_c4 0.0
# X_pro 0.0
# X_ac 0.0
# X_h2 0.0
# X_I 17010.64239199998
# S_cat 0.0
# S_an 5.210421371425576
# H2O 966179.9215060878

# MATLAB
# ADM1 influent (post ASM2ADM interface)
# ************************************* Ssu = monosacharides (kg COD/m3) = 0
# Saa = amino acids (kg COD/m3) = 0.04388
# Sfa = long chain fatty acids (LCFA) (kg COD/m3) = 0
# Sva = total valerate (kg COD/m3) = 0
# Sbu = total butyrate (kg COD/m3) = 0
# Spro = total propionate (kg COD/m3) = 0
# Sac = total acetate (kg COD/m3) = 0
# Sh2 = hydrogen gas (kg COD/m3) = 0
# Sch4 = methane gas (kg COD/m3) = 0
# Sic = inorganic carbon (kmole C/m3) = 0.0079326
# Sin = inorganic nitrogen (kmole N/m3) = 0.0019721 
# Si = soluble inerts (kg COD/m3) = 0.028067
# Xc = composites (kg COD/m3) = 0
# Xch = carbohydrates (kg COD/m3) = 3.7236
# Xpr = proteins (kg COD/m3) = 15.9235
# Xli = lipids (kg COD/m3) = 8.047
# Xsu = sugar degraders (kg COD/m3) = 0
# Xaa = amino acid degraders (kg COD/m3) = 0
# Xfa = LCFA degraders (kg COD/m3) = 0
# Xc4 = valerate and butyrate degraders (kg COD/m3) = 0
# Xpro = propionate degraders (kg COD/m3) = 0
# Xac = acetate degraders (kg COD/m3) = 0
# Xh2 = hydrogen degraders (kg COD/m3) = 0
# Xi = particulate inerts (kg COD/m3) = 17.0106
# Scat+ = cations (base) (kmole/m3) = 0
# San- = anions (acid) (kmole/m3) = 0.0052101
# Flow rate (m3/d) = 178.4674
# Temperature (degC) = 35
