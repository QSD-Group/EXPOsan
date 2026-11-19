# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""CHONS-conserving scenario chemistry and explicit product characterization."""
import numpy as np

ELEMENTS = ('C', 'H', 'O', 'N', 'S')
from chemicals.elements import periodic_table
ATOMIC = np.array([periodic_table[e].MW for e in ELEMENTS])
# https://github.com/QSD-Group/EXPOsan/blob/main/exposan/saf/_components.py
# Preserve HTL biochemical formulas and use deterministic CAS identities.
# Formula-only protein lookup is ambiguous across database/cache versions;
# resolved CAS/identity is exported. Fibers partition carbohydrate mass.
FEED_FORMULAS = dict(Lipids='C16H32O2', Proteins='C3H7NO2', Carbohydrates='C6H12O6',
                     Cellulose='C6H10O5', Hemicellulose='C5H8O4', Lignin='C20H22O7')


def atom_vector(chemical):
    return np.array([chemical.atoms.get(k, 0.) for k in ELEMENTS])


def empirical_formula(atoms):
    scale = atoms[0] if atoms[0] > 1e-12 else max(atoms)
    return ''.join(f'{e}{n/scale:.12g}' for e, n in zip(ELEMENTS, atoms) if n > 1e-12)


def pseudo(ID, formula, phase, density, cp=1.5):
    import thermosteam as tmo
    chemical = tmo.Chemical.blank(ID, formula=formula, phase=phase)
    chemical.default()
    # https://thermosteam.readthedocs.io/en/latest/Chemical.html : explicit
    # pseudocomponent property models. Density and 1.5 kJ/(kg K) Cp are engineering
    # scenarios, not measured oil/char properties; no molecular vapor pressure fit.
    chemical.V.add_model(chemical.MW/density/1000, top_priority=True)
    chemical.Cn.add_model(cp*chemical.MW, top_priority=True)
    # https://doi.org/10.1016/S0016-2361(01)00186-7 : Channiwala-Parikh HHV,
    # MJ/kg; C,H,O,N,S are dry mass percent. Ash is a separate inert component.
    fractions = atom_vector(chemical)*ATOMIC/chemical.MW
    C, H, O, N, S = fractions*100
    hhv = max(.3491*C+1.1783*H+.1005*S-.1034*O-.0151*N, 0.)
    a = atom_vector(chemical)
    chemical.Hf = a[0]*(-393522.) + a[1]/2*(-285830.) + a[4]*(-296842.) + hhv*1000*chemical.MW
    chemical.HHV = hhv*1000*chemical.MW
    chemical.LHV = max(chemical.HHV-a[1]/2*44010., 0.)
    return chemical


def create_thermo(feedstock, config, technology):
    """Return fresh local thermo; never replace the caller's global thermo."""
    import thermosteam as tmo
    if not feedstock.composition:
        raise ValueError('Use Feedstock.from_name or a consistent wet composition')
    composition = dict(feedstock.composition)
    fibers = feedstock.carbohydrate_fractions
    if fibers:
        if sum(fibers.values()) > composition['Carbohydrates']+1e-10:
            raise ValueError('Fiber partition exceeds carbohydrate umbrella')
        composition['Carbohydrates'] -= sum(fibers.values())
        composition.update(fibers)
    unknown = set(composition)-set(FEED_FORMULAS)-{'Water','Ash'}
    if unknown:
        raise ValueError(f'Unsupported feed fractions: {sorted(unknown)}')
    chemicals = [tmo.Chemical('Water')]
    # https://thermosteam.readthedocs.io/en/latest/Chemical.html : phase locks
    # prevent fictitious liquid-gas property extrapolation for noncondensables.
    # NH3/H2S dissolution is an empirical split; no electrolyte EOS is claimed.
    chemicals.extend(tmo.Chemical(i,phase='g') for i in ('CO2','CO','CH4','H2','NH3','H2S','O2','N2','SO2'))
    formulas = dict(FEED_FORMULAS)
    if config.ultimate_daf:
        formula = empirical_formula(np.array([config.ultimate_daf.get(e,0.) for e in ELEMENTS])/ATOMIC)
        formulas = {name:formula for name in formulas}
    from chemicals.identifiers import search_chemical
    # https://webbook.nist.gov/cgi/cbook.cgi?ID=C56417&Mask=1 : L-alanine,
    # C3H7NO2. Explicit CAS avoids formula-isomer lookup ambiguity in the parent
    # HTL source. Palmitic acid and D-glucose retain their native feed identities.
    identities = dict(Lipids='57-10-3', Proteins='56-41-7', Carbohydrates='50-99-7')
    for name, formula in formulas.items():
        if name in identities and not config.ultimate_daf:
            identity = search_chemical(identities[name])
            chemical = tmo.Chemical(name, search_ID=identity.CASs, phase='s')
            chemical.default()
            # Effective feed lumps decompose during heating; extrapolating the
            # stable molecular solid Cp to 500-650 C is outside its domain.
            # https://thermosteam.readthedocs.io/en/latest/Chemical.html
            # 1.5 kJ/(kg K) is an explicit shared screening heat-capacity scenario.
            chemical.Cn.add_model(1.5*chemical.MW,top_priority=True)
            chemicals.append(chemical)
        else:
            chemicals.append(pseudo(name, formula, 's', 1400.))
    ash = tmo.Chemical.blank('Ash', phase='s', MW=1.)
    ash.default()
    # https://github.com/QSD-Group/EXPOsan/blob/main/exposan/saf/_components.py
    # Ash is an inert lump; 2500 kg/m3 and 0.8 kJ/(kg K) are uncalibrated
    # mineral-mixture scenarios, not a claim about silica or measured sludge ash.
    ash.V.add_model(1/2500000., top_priority=True)
    ash.Cn.add_model(.8, top_priority=True)
    ash.Hf = ash.HHV = ash.LHV = 0.
    chemicals.append(ash)
    by_name = {i.ID:i for i in chemicals}
    atoms = sum((mass/by_name[name].MW*atom_vector(by_name[name]) for name,mass in composition.items()
                 if name in FEED_FORMULAS), np.zeros(5))
    products = {}
    if technology in ('fast_pyrolysis','htc'):
        for ID, retention, phase, rho in (
                ('PyrolysisOil', config.oil_retention, 'l', 1200.),
                ('BiocharMatrix' if technology == 'fast_pyrolysis' else 'HydrocharMatrix', config.char_retention, 's', 1500.),
                ('AqueousOrganics', config.aqueous_retention, 'l', 1100.)):
            retained = atoms*np.array([retention.get(e,0.) for e in ELEMENTS])
            if retained.sum() > 1e-14:
                chemical = pseudo(ID, empirical_formula(retained), phase, rho)
                chemicals.append(chemical)
                products[ID] = dict(atoms_per_kg_wet=retained, mass_per_kg_wet=float(retained@ATOMIC))
                if ID == 'PyrolysisOil':
                    # https://chemicals.readthedocs.io/chemicals.phase_change.html
                    # Native database Hvap at reference normal boiling point;
                    # J/mol divided by g/mol is numerically kJ/kg. This is an
                    # explicit property surrogate, not a molecular oil analysis.
                    reference = tmo.Chemical(config.oil_latent_reference)
                    if not reference.Tb or not reference.Hvap(reference.Tb):
                        raise ValueError('Oil latent reference needs a database Tb and Hvap')
                    products[ID].update(latent_kj_kg=float(reference.Hvap(reference.Tb)/reference.MW),
                        latent_reference=config.oil_latent_reference,reference_CAS=reference.CAS,
                        reference_boiling_k=reference.Tb)
    else:
        for name in FEED_FORMULAS:
            if composition.get(name,0.):
                products[name] = dict(mass_per_kg_wet=composition[name]*(1-config.gasification_conversion))
    chemicals = tmo.Chemicals(chemicals)
    chemicals.compile()
    return tmo.Thermo(chemicals), composition, products, atoms


def gasify(atoms, completion, co_fraction):
    """Stoichiometric reforming plus an explicit methanation extent, kmol/h.

    This is a constrained reaction scenario, not Gibbs equilibrium or kinetics.
    Positive reaction water is consumed from the incoming liquid inventory.
    """
    # https://www.psi.ch/en/cpe/hydropilot : gas species and aqueous N boundary.
    # CHONS balance derives reforming; completion/CO partition are explicit inputs.
    C,H,O,N,S = atoms
    co = C*co_fraction
    hydrogen = H/2-1.5*N-S+2*C-O-co
    if hydrogen < -1e-9:
        raise ValueError('Selected elemental feed cannot support this reforming scenario')
    methane = min(max(hydrogen,0.)/4, C-co)*completion
    return dict(CO=co, CO2=C-co-methane, CH4=methane,
                H2=max(hydrogen-4*methane,0.), NH3=N, H2S=S,
                reaction_water=2*C-O-co-2*methane)


def complete_gas_balance(atoms, water_kmol):
    """Nonnegative atom balance with a stated carbon-gas allocation objective.

    Minimize departure from CO/CO2/CH4 = 0.4/0.5/0.1 of residual carbon.
    The result completes a scenario balance; it is not measured gas speciation.
    """
    from scipy.optimize import linprog
    C,H,O,N,S = atoms
    species = ('CO','CO2','CH4','H2','Water','NH3','H2S')
    A = np.array(((1,1,1,0,0,0,0), (0,0,4,2,2,3,2),
                  (1,2,0,0,1,0,0), (0,0,0,0,0,1,0), (0,0,0,0,0,0,1)), dtype=float)
    b = np.array((C,H+2*water_kmol,O+water_kmol,N,S))
    # https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html
    # 0.4/0.5/0.1 is a numerical allocation objective, not measured composition.
    target = np.array((.4*C,.5*C,.1*C))
    eq = np.zeros((8,13))
    eq[:5,:7] = A
    for j in range(3):
        eq[5+j,j] = 1
        eq[5+j,7+j] = -1
        eq[5+j,10+j] = 1
    objective = np.r_[np.zeros(7),np.ones(6)]
    objective[3] = 1e-6
    solution = linprog(objective, A_eq=eq, b_eq=np.r_[b,target], bounds=(0,None), method='highs')
    if not solution.success:
        raise ValueError('Product retentions cannot close CHONS with the available feed water; revise characterization/scenario')
    values = np.maximum(solution.x[:7],0.)
    if np.max(np.abs(A@values-b)) > max(1e-9, np.max(np.abs(b))*1e-8):
        raise ValueError('Gas balance solver did not close the represented elements')
    return dict(zip(species,values))


def characterize_product(stream, target=None):
    """Whole-product water, dry ultimate/ash, HHV and dry gas mole fractions."""
    mass = float(stream.F_mass)
    water = float(stream.imass['Water'])
    dry = mass-water
    ash = float(stream.imass['Ash'])
    atoms = sum((float(stream.imol[c.ID])*atom_vector(c) for c in stream.chemicals if c.ID != 'Water'), np.zeros(5))
    element_mass = atoms*ATOMIC
    gas_ids = ('CH4','H2','CO2','CO','NH3','H2S','N2','O2','SO2')
    gas_moles = {k:float(stream.imol[k]) for k in gas_ids}
    gas_total = sum(gas_moles.values())
    hhv = sum(float(stream.imol[c.ID])*float(c.HHV or 0.) for c in stream.chemicals)/1000
    dry_mol = float(stream.F_mol)-float(stream.imol['Water'])
    return dict(mass_kg_h=mass, water_wt=water/mass if mass else 0.,
        ash_dry_wt=ash/dry if dry else None,
        ultimate_dry_wt={e:float(v/dry) if dry else None for e,v in zip(ELEMENTS,element_mass)},
        hhv_mj_kg=hhv/mass if mass else 0., hhv_mj_h=hhv,
        dry_gas_mol={k:v/gas_total for k,v in gas_moles.items()} if gas_total else {},
        purity_mol=float(stream.imol[target])/float(stream.F_mol) if target and stream.F_mol else None,
        normal_volume_nm3_h=gas_total*22.414,
        normal_reference='0 C, 101325 Pa, ideal gas',
        proximate_analysis='not inferred; volatile matter and fixed carbon require measurement',
        molecular_speciation='explicit gas species; condensed pseudocomponents are elemental surrogates')
