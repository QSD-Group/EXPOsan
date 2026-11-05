"""Runnable isolated ED example: python -m exposan.VFA.ed_example.

Every coefficient below is ILLUSTRATIVE, not a validated VFA membrane dataset.
Original VFA systems/components/RedoxED are deliberately left unchanged.
"""
import numpy as np
from qsdsan import WasteStream, Flowsheet
from thermosteam import Thermo
from ._ed_components import create_ed_components, ION_CHARGES
from ._electrodialysis import Electrodialysis


def create_ed_feeds(thermo, diluate_water=100., concentrate_water=100.,
                    vfa_kmol_hr=None, chloride_kmol_hr=0.01,
                    concentrate_nacl_kmol_hr=0.01, potassium_kmol_hr=0.):
    """Construct already-ionized, charge-balanced liquid feeds [kg/hr, kmol/hr].

    Default VFAs: 0.003 kmol/hr each; all values are example assumptions.
    Sodium supplies the remaining positive equivalents after specified K.
    No neutral acid is converted and no base consumption is calculated.
    TODO: connect a separate, mass-balanced pH adjustment/deprotonation model.
    """
    vf = dict(Propionate=0.003, Butyrate=0.003, Hexanoate=0.003) if vfa_kmol_hr is None else dict(vfa_kmol_hr)
    if set(vf) - {'Propionate', 'Butyrate', 'Hexanoate'}:
        raise ValueError('Only the three VFA anions are accepted.')
    values = [*vf.values(), chloride_kmol_hr, concentrate_nacl_kmol_hr, potassium_kmol_hr]
    if not np.all(np.isfinite(values)) or min(values) < 0:
        raise ValueError('Feed molar flows must be finite and nonnegative.')
    if not np.isfinite(diluate_water) or not np.isfinite(concentrate_water) or min(diluate_water, concentrate_water) <= 0:
        raise ValueError('Water flows must be positive and finite.')
    sodium = chloride_kmol_hr + sum(vf.values()) - potassium_kmol_hr
    if sodium < 0:
        raise ValueError('Potassium exceeds available anionic equivalents.')
    d = WasteStream('ed_diluate_in', thermo=thermo, H2O=diluate_water, units='kg/hr')
    c = WasteStream('ed_concentrate_in', thermo=thermo, H2O=concentrate_water, units='kg/hr')
    for ID, flow in vf.items():
        d.imol[ID] = flow
    d.imol['Na'], d.imol['K'], d.imol['Cl'] = sodium, potassium_kmol_hr, chloride_kmol_hr
    c.imol['Na'] = c.imol['Cl'] = concentrate_nacl_kmol_hr
    return d, c


def create_example(simulate=True):
    """Return an isolated ED unit using provisional inputs, without global thermo changes.

    TODO: replace ALL numerical transport/electrical/geometry inputs below
    with measured or cited values for the chosen membrane and feed.
    """
    with Flowsheet('vfa_ed_example'):
        cmps = create_ed_components(set_thermo_flag=False)
        thermo = Thermo(cmps)
        feeds = create_ed_feeds(thermo)
        # CEM: Na carries current; AEM: equal VFA fractions, plus chloride.
        # K is tracked even though absent in these illustrative feeds.
        t = [[1., 0., 0., 0., 0., 0.], [0., 0., 0.25, 0.25, 0.25, 0.25]]
        ed = Electrodialysis('electrodialysis', ins=feeds,
            outs=(WasteStream('ed_diluate_out', thermo=thermo),
                  WasteStream('ed_concentrate_out', thermo=thermo)), thermo=thermo,
            charges=ION_CHARGES, target='Propionate', recovery=0.5,
            current_density=10., width=0.1, depth=0.001,
            transport_numbers=t, permeabilities=np.full((2, 6), 1e-8),
            ionic_conductivities=np.full(6, 0.005), limiting_coefficient=1.,
            membrane_resistance=[2e-4, 2e-4],
            membrane_resistance_coefficient=[0.01, 0.01])
        # Equal summed permeabilities preserve zero-charge diffusion for this
        # electroneutral example; arbitrary ion values need a coupled closure.
        if simulate:
            ed.simulate()
        return ed


if __name__ == '__main__':
    unit = create_example()
    print(unit.design_results)
    print('Signed ion transfers [mol/s]:', unit.transfers)
