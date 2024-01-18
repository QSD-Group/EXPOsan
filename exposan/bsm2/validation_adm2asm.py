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
    K_O_H=0.2,
    K_NO=0.5,
    b_H=0.3, #0.62;
    mu_A=0.5, #0.8;
    K_NH=1.0,
    K_O_A=0.4,
    b_A=0.05, #0.2;
    eta_g=0.8, # ny_g
    k_a=0.05, #0.08;
    k_h=3.0,
    K_X=0.1, #0.03;
    eta_h=0.8, # ny_h #0.4;
    # path=os.path.join(data_path, '_asm1.tsv'),
    )

cmps_asm1 = pc.create_asm1_cmps()
asm1 = pc.ASM1(components=cmps_asm1, **default_asm_kwargs)
thermo_asm1 = qs.get_thermo()
effluent = WasteStream('effluent', T=14.8581+273.15)

cmps_adm1 = pc.create_adm1_cmps()
thermo_adm1 = qs.get_thermo()
adm1 = pc.ADM1()


# %%

qs.set_thermo(thermo_adm1)

# MATLAB
# ADM1 effluent (pre ADM2ASM interface)
# ***************************************
# Ssu = monosacharides (kg COD/m3) = 0.012394
# Saa = amino acids (kg COD/m3) = 0.0055432
# Sfa = long chain fatty acids (LCFA) (kg COD/m3) = 0.10741
# Sva = total valerate (kg COD/m3) = 0.012333
# Sbu = total butyrate (kg COD/m3) = 0.014003
# Spro = total propionate (kg COD/m3) = 0.017584
# Sac = total acetate (kg COD/m3) = 0.089315
# Sh2 = hydrogen gas (kg COD/m3) = 2.5055e-07
# Sch4 = methane gas (kg COD/m3) = 0.05549
# Sic = inorganic carbon (kmole C/m3) = 0.095149
# Sin = inorganic nitrogen (kmole N/m3) = 0.094468 (= 1.3226 kg N/m3)
# Si = soluble inerts (kg COD/m3) = 0.13087
# Xc = composites (kg COD/m3) = 0.10792
# Xch = carbohydrates (kg COD/m3) = 0.020517
# Xpr = proteins (kg COD/m3) = 0.08422
# Xli = lipids (kg COD/m3) = 0.043629
# Xsu = sugar degraders (kg COD/m3) = 0.31222
# Xaa = amino acid degraders (kg COD/m3) = 0.93167
# Xfa = LCFA degraders (kg COD/m3) = 0.33839
# Xc4 = valerate and butyrate degraders (kg COD/m3) = 0.33577
# Xpro = propionate degraders (kg COD/m3) = 0.10112
# Xac = acetate degraders (kg COD/m3) = 0.67724
# Xh2 = hydrogen degraders (kg COD/m3) = 0.28484
# Xi = particulate inerts (kg COD/m3) = 17.2162
# Scat+ = cations (base) (kmole/m3) = -4.0789e-34
# San- = anions (acid) (kmole/m3) = 0.0052101
# Flow rate (m3/d) = 178.4674
# Temperature (degC) = 35

# pH = pH within AD system = 7.2631
# S_H+ = protons (kmole/m3) = 5.4562e-08
# Sva- = valerate (kg COD/m3) = 0.012284
# Sbu- = butyrate (kg COD/m3) = 0.013953
# Spro- = propionate (kg COD/m3) = 0.017511
# Sac- = acetate (kg COD/m3) = 0.089035
# Shco3- = bicarbonate (kmole C/m3) = 0.08568
# Sco2 = carbon dioxide (kmole C/m3) = 0.0094689
# Snh3 = ammonia (kmole N/m3) = 0.001884
# Snh4+ = ammonium (kmole N/m3) = 0.092584
# Sgas,h2 = hydrogen concentration in gas phase (kg COD/m3) = 1.1032e-05
# Sgas,ch4 = methane concentration in gas phase (kg COD/m3) = 1.6535
# Sgas,co2 = carbon dioxide concentration in gas phase (kmole C/m3) = 0.01354
# pgas,h2 = partial pressure of hydrogen gas (bar, true value i.e. not normalized) = 1.7666e-05
# pgas,ch4 = partial pressure of methane gas (bar, true value i.e. not normalized) = 0.66195
# pgas,co2 = partial pressure of carbon dioxide gas (bar, true value, i.e. not normalized) = 0.34691
# pgas,total = total head space pressure of H2+CO2+CH4+H2O (bar, true value, i.e. not normalized) = 1.0645
# qgas = gas flow rate normalized to atmospheric pressure (m3/d) = 2708.3431
#
# Extra calculated outputs
# ------------------------
# Produced hydrogen gas (kg H2/d) = 0.0035541
# Produced methane gas (kg CH4/d) = 1065.3523
# Produced carbon dioxide gas (kg CO2/d) = 1535.4118
# Energy content of methane gas (MJ/d) = 53282.5305
# Energy content of methane gas (kWh/d) = 14800.7029

T_inf = 14.8581 + 273.15
default_inf_kwargs = {
    'flow_tot': 178.4674,
    'concentrations': {
        'S_su': 0.012394,
        'S_aa': 0.0055432,
        'S_fa': 0.10741,
        'S_va': 0.012333,
        'S_bu': 0.014003,
        'S_pro': 0.017584,
        'S_ac': 0.089315,
        'S_h2': 2.5055e-07,
        'S_ch4': 0.05549,
        'S_IC': 0.095149*12,
        'S_IN': 1.3226,
        'S_I': 0.13087,
        'X_c': 0.10792,
        'X_ch': 0.020517,
        'X_pr': 0.08422,
        'X_li': 0.043629,
        'X_su': 0.31222,
        'X_aa': 0.93167,
        'X_fa': 0.33839,
        'X_c4': 0.33577,
        'X_pro': 0.10112,
        'X_ac': 0.67724,
        'X_h2': 0.28484,
        'X_I': 17.2162,
        # Scat+ = cations (base) (kmole/m3) = -4.0789e-34
        # San- = anions (acid) (kmole/m3) = 0.0052101
        },
    'units': ('m3/d', 'kg/m3'),
    }



influent = WasteStream('influent', T=T_inf)

influent.set_flow_by_concentration(**default_inf_kwargs)

J1 = su.ADMtoASM('J1', upstream=influent, downstream=effluent, 
                 thermo=thermo_asm1, isdynamic=True, adm1_model=adm1)
J1.bio_to_xs = 0.79
sys = qs.System('sys', path=(J1,))

t = 1 # simulation time shouldn't matter
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

xN_0 = influent.composite('N', particle_size='x')
sN_0 = influent.composite('N', particle_size='s')

xN_qs = effluent.composite('N', particle_size='x')
sN_qs = effluent.composite('N', particle_size='s')


# effluent_conc, #!!! to be updated
# {'S_I': 130.87,
#  'S_S': 258.5822000000001,
#  'X_I': 17216.200000000008,
#  'X_S': 2611.4735000000014,
#  'X_BH': 0.0,
#  'X_BA': 0.0,
#  'X_P': 626.0625000000002,
#  'S_O': 0.0,
#  'S_NO': 0.0,
#  'S_NH': 1443.3491250362472,
#  'S_ND': 0.5434935760800002,
#  'X_ND': 100.92453092319155,
#  'S_ALK': 1180.3134876913525,
#  'S_N2': 0.0,
#  'H2O': 982527.7939584954}

xN_mt = (17216.2434 + 626.0652) * 0.06 + 100.8668
sN_mt = 1442.7882 + 0.54323

# MATLAB
# Anaerobic digester output (post ADM2ASM interface)
# **************************************************
# SI = 130.867 mg COD/l
# SS = 258.5789 mg COD/l
# XI = 17216.2434 mg COD/l
# XS = 2611.4843 mg COD/l
# XBH = 0 mg COD/l
# XBA = 0 mg COD/l
# XP = 626.0652 mg COD/l
# SO = 0 mg -COD/l
# SNO = 0 mg N/l
# SNH = 1442.7882 mg N/l
# SND = 0.54323 mg N/l
# XND = 100.8668 mg N/l
# SALK = 97.8459 mol HCO3/m3 (* 12 = 1174.1508 gC/m3)
# TSS = 15340.3447 mg SS/l
# Flow rate = 178.4674 m3/d
# Temperature = 14.8581 degC