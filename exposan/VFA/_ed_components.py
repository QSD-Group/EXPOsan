"""Ionic components for the isolated VFA ED example; no pH/speciation model."""
from qsdsan import Component, Components, set_thermo
from exposan.utils import add_V_from_rho

__all__ = ('create_ed_components', 'ION_CHARGES')
ION_CHARGES = dict(Na=1, K=1, Cl=-1, Propionate=-1, Butyrate=-1, Hexanoate=-1)


def create_ed_components(set_thermo_flag=True):
    """Return water and explicit ions on a chemical-mass basis.

    VFA formulae exclude the acidic proton. Charges are supplied separately
    in ION_CHARGES for electrochemistry (no equilibrium dissociation implied).
    Volume/heat-capacity models are illustrative aqueous pseudo-properties;
    TODO: replace with partial molar volumes and electrolyte property data.
    This separate thermo must not be mixed with the original neutral-acid thermo.
    """
    water = Component('H2O', particle_size='Soluble',
                      degradability='Undegradable', organic=False)
    components = [water]
    formulas = dict(Na='Na', K='K', Cl='Cl', Propionate='C3H5O2',
                    Butyrate='C4H7O2', Hexanoate='C6H11O2')
    for ID, formula in formulas.items():
        c = Component(ID, search_db=False, formula=formula, phase='l',
                      particle_size='Soluble', degradability='Undegradable',
                      organic=ID in ('Propionate', 'Butyrate', 'Hexanoate'))
        c.i_charge = ION_CHARGES[ID] / c.MW  # mol positive charge/g chemical
        # ASSUMPTION: 1000 kg/m3 apparent solute density, not ionic solid density.
        add_V_from_rho(c, 1000.)
        c.Cn.add_model(75.)  # J/mol/K; illustrative, no thermal sizing here.
        c.mu.add_model(0.001)
        c.Hf = c.HHV = c.LHV = 0.
        components.append(c)
    cmps = Components(components)
    cmps.compile()
    cmps.set_alias('H2O', 'Water')
    if set_thermo_flag:
        set_thermo(cmps)
    return cmps
