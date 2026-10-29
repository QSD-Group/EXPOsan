# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Python evaluation of the published PaCTEA active-composting equations.

Reference: Ravoahangy et al., Waste Management 212 (2026), 115337.
https://doi.org/10.1016/j.wasman.2026.115337
Equation/parameter reference:
https://github.com/nomenazo/PaCTEA/blob/316fa8f9cc809a7bc69d2c0f7d3c46b0121d9fd7/composting.m
State units are kg/kg initial mixture, temperature K, and time hours.
The reference equations include rounded/nonconservative stoichiometry. Report
the resulting mass/element residuals; do not silently normalize product yields.
"""
from dataclasses import dataclass
from math import isfinite
import numpy as np
from scipy.integrate import solve_ivp

SOURCE = 'https://doi.org/10.1016/j.wasman.2026.115337'
STATES = ('C','P','L','H','CE','LG','Xi','Sc','Sp','Sl','Sh','Slg',
          'Xmb','Xtb','Xma','Xta','Xmf','Xtf','Xdb','CO2','W','T',
          'CH4gen','CH4oxi','CH4','Xa','NO3','N2O','N2','NH3','NH4',
          'Wvap','O2diss','O2gaz','O2cons','O2out','NH3gaz')
SUBSTRATES = ('Carbohydrates','Proteins','Lipids','Hemicellulose','Cellulose','Lignin')
SEED = np.array([5e-4,5e-4,1e-4,1e-4,1e-6,1e-6,1e-4])
# PaCTEA reference substrate, biomass, and formed-inert elemental mass fractions.
CF = np.array([72/180,192/352,300/393,120/282,72/180,240/366,
               1.48*12/34.54,72/180,192/352,300/393,120/282,240/366,
               60/113,60/113,60/113,60/113,120/247,120/247,60/113])
NF = np.array([0,56/352,0,0,0,0,.12*14/34.54,0,56/352,0,0,0,
               14/113,14/113,14/113,14/113,14/247,14/247,14/113])


@dataclass(frozen=True)
class PaCTEAConfig:
    """Reference defaults; pile parameters are independent of plant throughput."""
    ambient_temperature_k: float = 298.
    aeration_l_min_kg: float = .3
    bulk_density_kg_m3: float = 750.
    free_air_space: float = .21
    pile_volume_m3: float = .000135
    heat_transfer_area_m2: float = .1
    heat_transfer_kj_m2_k_h: float = .09
    initial_dissolved_oxygen: float = .003
    initial_gas_oxygen: float = .05
    initial_ammonium: float = 0.
    relative_tolerance: float = 1e-7
    absolute_tolerance: float = 1e-10
    # https://github.com/nomenazo/PaCTEA/blob/main/composting.m h^-1
    growth_rates: tuple = (.2,.18,.1,.12,.1,.1,.003)
    death_rates: tuple = (.012,.008,.01,.015,.01,.01,.0083)

    def __post_init__(self):
        for key,value in self.__dict__.items():
            if isinstance(value, tuple):
                if len(value) != 7 or any(not isfinite(v) or v < 0 for v in value):
                    raise ValueError(f'{key} must contain seven nonnegative finite rates')
            elif not isfinite(value) or value < 0:
                raise ValueError(f'{key} must be nonnegative and finite')
        for key in ('ambient_temperature_k','bulk_density_kg_m3','pile_volume_m3',
                    'free_air_space','relative_tolerance','absolute_tolerance'):
            if getattr(self,key) <= 0:
                raise ValueError(f'{key} must be positive')
        if self.free_air_space >= 1:
            raise ValueError('free_air_space must be less than one')


_CO2 = np.array([.445717506,.26543624,.165284084]*3 + [.139609659] +
                [.445717506,.26543624,.165284084,.139609659] +
                [.41509283,.254264706,.160882526,.136456284,.19653097]*2)
_WATER = np.array([.716485507,1.791988467,.380836347]*3 + [.329440099] +
                  [.716485507,1.791988467,.380836347,.329440099] +
                  [.830451489,2.728498673,.410802056,.351627849,.625936213]*2)
_O2 = np.array([.61286157,.34845815,.156798785]*3 + [.191963281] +
               [.61286157,.34845815,.156798785,.191963281] +
               [.592664534,.341834717,.148720799,.189935871,.210088681]*2)
_O2MOL = np.array([.173553719,.098678414,.044403196]*3 + [.054361283] +
                  [.173553719,.098678414,.044403196,.054361283] +
                  [.07678245,.044286279,.019267472,.024607076,.027217967]*2)
_NH4 = np.array([6.64705882352941,2.490625447,6.647058824,
                 6.64705882352941,2.490625447,6.647058824,
                 .705882353,2.490625447,6.647058824,6.647058824,
                 .705882353,2.490625447,6.647058824,6.647058824] +
                [14.52941176,2.069869946,14.52941176,14.52941176,14.52941176]*2)
_GROUP = np.repeat(np.arange(6), [3,3,4,4,5,5])
_SOLUBLE = np.array([0,1,2,0,1,2,0,1,2,3,0,1,2,3,0,1,2,3,4,0,1,2,3,4])
_HYDRO_S = np.array([0,1,2,0,1,2,3,4,5,3,4,5])
_HYDRO_X = np.array([0,0,0,1,1,1,3,5,5,2,4,4])
_KH = np.array([.0293,.1508*.5,1.4563e-7,.0053,.1731*.5,.0182,
                .009,.007,.0068,.009,.007,.0068])


def resolve_carbohydrates(composition, fractions=None):
    """Resolve optional children of the carbohydrate umbrella for the ODE.

    Without fiber characterization, assign the umbrella to the reference C
    pool as an explicit lumped approximation. Structural-fiber kinetics and
    recalcitrance then cannot be inferred from the input data.
    """
    result = dict(composition)
    keys = ('Cellulose','Hemicellulose','Lignin')
    if fractions:
        if set(fractions) != set(keys):
            raise ValueError('Supply all three carbohydrate child fractions')
        if any(not isfinite(v) or v < 0 for v in fractions.values()):
            raise ValueError('Carbohydrate child fractions must be finite and nonnegative')
        remainder = result['Carbohydrates']-sum(fractions.values())
        if remainder < -1e-9:
            raise ValueError('Carbohydrate child fractions exceed their umbrella total')
        result['Carbohydrates'] = max(remainder,0.)
        result.update(fractions)
        return result, 'resolved carbohydrate children'
    if any(key in result for key in keys):
        if not all(key in result for key in keys):
            raise ValueError('Detailed PaCTEA input requires all three fiber pools')
        return result, 'resolved biochemical pools'
    result.update(dict.fromkeys(keys,0.))
    return result, 'lumped carbohydrate approximation; fiber kinetics unresolved'


def rhs(t, state, config):
    """Evaluate the 37 published state equations (nonnegative rate arguments)."""
    y = np.maximum(state, 0.)
    C,P,L,H,CE,LG,Xi,Sc,Sp,Sl,Sh,Slg,Xmb,Xtb,Xma,Xta,Xmf,Xtf,Xdb,CO2,W,T,CH4gen,CH4oxi,CH4,Xa,NO3,N2O,N2,NH3,NH4,Wvap,O2diss,O2gaz,O2cons,O2out,NH3gaz = y
    W = max(W,1e-12)
    T = max(T,200.)
    tm = config.pile_volume_m3*config.bulk_density_kg_m3
    vg = config.free_air_space*config.pile_volume_m3
    qair = config.aeration_l_min_kg*.001*tm*1.2*60
    x = y[12:18]
    vhydro = _KH*x[_HYDRO_X]*y[_HYDRO_S]/np.maximum(.0001*x[_HYDRO_X]+y[_HYDRO_S],1e-30)
    ft = ((T-355)*(T-273)**2)/(49*(49*(T-322)+33*(595-2*T)))
    ft = max(0.,ft) if 273 < T < 355 else 0.
    soluble = y[7:12]/W
    oxygen = O2diss/W
    growth = np.asarray(config.growth_rates)[:6][_GROUP]*x[_GROUP]*soluble[_SOLUBLE]/(6.2e-5+soluble[_SOLUBLE])*ft*oxygen/(.007e-6+oxygen)
    death = np.asarray(config.death_rates)[:6]*x
    lysis = .0025*Xdb
    kn = 10**(.051*(T-275)-7.158)
    nitrification = config.growth_rates[6]*Xa*(NH4/W)/(kn+NH4/W)*oxygen/(.0608e-5+oxygen)
    ftdenit = np.exp((((T-275)-(11 if T < 286 else 20))*np.log(89 if T < 286 else 2.1)-(9*np.log(2.1) if T < 286 else 0))/10)
    denitrification = .01*NO3*(NO3/W)/(.084+NO3/W)*ftdenit
    autodeath = config.death_rates[6]*Xa
    gas = max(CO2+CH4+N2+NH3gaz+N2O+O2gaz,1e-30)
    volatilization = qair*NH3gaz/(gas*tm)
    oxidized = 5.35e-4*(CH4gen/W)/(.72+CH4gen/W)*(O2diss/(W*.032))/(.033+O2diss/(W*.032))
    hydro_totals = np.bincount(_HYDRO_S, weights=vhydro,minlength=6)
    soluble_hydro = hydro_totals[[0,1,2,3,5]].copy()
    soluble_hydro[0] += hydro_totals[4]
    d = np.zeros(37)
    d[:6] = -hydro_totals
    d[1] += .7*(lysis+autodeath)
    d[6] = .31*(lysis+autodeath)
    d[7:12] = soluble_hydro-np.bincount(_SOLUBLE,weights=growth/.35,minlength=5)
    d[12:18] = np.bincount(_GROUP,weights=growth,minlength=6)-death
    d[18] = sum(death)-lysis
    d[19] = np.sum(growth/_CO2)+2.75*oxidized
    # Preserve reference evaporation and oxygen-transfer equations. Wvap is a
    # gas inventory, not cumulative exhaust; integrate exhaust separately.
    ps = 10**(22.443-2795/T-1.6798*np.log(T))
    evaporation = 1e-4/tm*(ps-Wvap/.018*8.134*T/vg)
    vapor_out = qair*Wvap/(1.2*vg*tm)
    d[31] = evaporation-vapor_out
    d[20] = np.sum(growth/_WATER)-d[31]+2.25*oxidized
    o2eq = 1.3e-5*np.exp(12000/8.314*(1/T-1/298))*(O2gaz*8.134*T/(vg*.032))*W*.001*.032*tm
    transfer = 180*1.02**(T-293)*(o2eq-O2diss)
    consumption = np.sum(growth/_O2)+nitrification/.040522541
    d[32] = transfer-consumption-4*oxidized
    d[34] = consumption
    d[35] = qair*O2gaz/(gas*tm*1.2)
    d[33] = qair*.23/tm-transfer-d[35]
    qbio = 440*tm*np.sum(growth/_O2MOL/np.where(_GROUP<4,.113,.247))
    d[21] = (qbio-(T-config.ambient_temperature_k)*(qair+config.heat_transfer_kj_m2_k_h*config.heat_transfer_area_m2))/(tm*(2+4*W))
    d[22] = np.dot([.267,.375,.707,.284,.535],soluble_hydro)/(1+4e5*(O2diss/(W*.032)))
    d[23] = oxidized
    d[24] = d[22]-d[23]
    d[25] = nitrification-autodeath
    d[26] = nitrification/.0390553-denitrification
    d[27] = .0710*denitrification
    d[28] = .181*denitrification
    d[29] = volatilization
    nh4eq = 6e-4*(NH3gaz*8.134*T/(vg*.017))*W*.017*tm
    nh4transfer = 100*np.exp(.03*(T-298))*(NH4-nh4eq)
    nh4_weights = -1/_NH4
    nh4_weights[_SOLUBLE==1] *= -1
    # Reference implementation uses CO2 coefficients for fungal protein terms.
    nh4_weights[15] = 1/_CO2[15]
    nh4_weights[20] = 1/_CO2[20]
    d[30] = -nitrification/.13+np.dot(nh4_weights,growth)-nh4transfer
    d[36] = nh4transfer-volatilization
    # MATLAB's published driver uses odeset('NonNegative',1:length(y0)).
    d[(state <= 0.) & (d < 0.)] = 0.
    return d


def simulate(composition, days, config=None):
    """Predict conversion for an explicitly fractionated wet feed/bulking mixture.

    Ash is conserved outside the organic inert pool. Inoculum is an explicit
    additional input at the reference kg/kg-mixture seed concentrations.
    """
    c = config or PaCTEAConfig()
    if not isfinite(days) or days <= 0:
        raise ValueError('Composting duration must be positive')
    missing = [key for key in SUBSTRATES if key not in composition]
    if missing:
        raise ValueError(f'PaCTEA requires explicit biochemical fractions: {missing}; supply fiber characterization in deploy/feedstocks.py or characterization=')
    if (any(not isfinite(v) or v < 0 for v in composition.values()) or
        abs(sum(composition.values())-1) > 1e-9 or
        set(composition)-set((*SUBSTRATES,'Water','Ash'))):
        raise ValueError('PaCTEA composition must contain only six substrates, Water and Ash, summing to one')
    y0 = np.zeros(37)
    y0[:6] = [composition[k] for k in SUBSTRATES]
    y0[12:18], y0[25] = SEED[:6], SEED[6]
    y0[20],y0[21] = composition['Water'],c.ambient_temperature_k
    y0[30],y0[32],y0[33] = c.initial_ammonium,c.initial_dissolved_oxygen,c.initial_gas_oxygen
    sol = solve_ivp(lambda t,y:rhs(t,y,c),(0.,days*24),y0,method='BDF',
                    rtol=c.relative_tolerance,atol=c.absolute_tolerance,
                    max_step=12.)
    if not sol.success or not np.isfinite(sol.y).all():
        raise RuntimeError(f'PaCTEA integration failed: {sol.message}')
    if np.min(sol.y) < -max(1e-7,100*c.absolute_tolerance):
        raise RuntimeError('PaCTEA produced negative inventories; revise characterization/operating conditions')
    end = np.maximum(sol.y[:,-1],0.)
    dry = sum(end[:19])+end[25]+end[26]+end[30]+composition['Ash']
    if dry <= 0:
        raise RuntimeError('PaCTEA predicts no dry product')
    initial_c = float(CF@y0[:19]+60/113*y0[25])
    initial_n = float(NF@y0[:19]+14/113*y0[25]+14/18*y0[30])
    retained_c = float(CF@end[:19]+60/113*end[25])
    organic_n = float(NF@end[:19]+14/113*end[25])
    retained_n = organic_n+14/18*end[30]+14/62*end[26]
    c_out = retained_c+12/44*end[19]+12/16*end[24]
    n_out = retained_n+14/17*(end[29]+end[36])+28/44*end[27]+end[28]
    return dict(model='PaCTEA',source=SOURCE,days=days,
                source_revision='316fa8f9cc809a7bc69d2c0f7d3c46b0121d9fd7',
                validation_status='reference_equations; not calibrated to this feedstock',
                dry_product_kg_kg_mixture=float(dry),
                model_water_kg_kg_mixture=float(end[20]),
                co2_kg_kg_mixture=float(end[19]),methane_kg_kg_mixture=float(end[24]),
                ammonia_kg_kg_mixture=float(end[29]+end[36]),
                n2o_kg_kg_mixture=float(end[27]),nitrogen_gas_kg_kg_mixture=float(end[28]),
                oxygen_consumed_kg_kg_mixture=float(end[34]+4*end[23]),
                inoculum_kg_kg_mixture=float(sum(SEED)),
                initial_carbon_kg_kg_mixture=initial_c,
                initial_nitrogen_kg_kg_mixture=initial_n,
                carbon_kg_kg_mixture=retained_c,
                organic_nitrogen_kg_kg_mixture=organic_n,
                ammonium_n_kg_kg_mixture=float(14/18*end[30]),
                nitrate_n_kg_kg_mixture=float(14/62*end[26]),
                nitrogen_kg_kg_mixture=retained_n,
                carbon_balance_residual_kg_kg_mixture=c_out-initial_c,
                nitrogen_balance_residual_kg_kg_mixture=n_out-initial_n,
                final_temperature_k=float(end[21]),
                final_states=dict(zip(STATES,map(float,end))))
