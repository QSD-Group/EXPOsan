# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

Standalone RWW -> SandFilter -> effluent system for ESAP validation.

The wastewater profiles and component properties are shared directly with
TricklingFilter_validation so both unit validations use identical feeds.

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import qsdsan as qs
from qsdsan.unit_operations import SandFiltration

from .TricklingFilter_validation import (
    _NONDEFAULT_COMPONENT_PROPERTIES as _NONDEFAULT_COMPONENT_PROPERTIES,
    _WASTEWATER_ORGANIC_FRACTIONS as _WASTEWATER_ORGANIC_FRACTIONS,
    _WASTEWATER_STRENGTHS as _WASTEWATER_STRENGTHS,
    _create_validation_components,
    _create_validation_wastewater,
    _normalize_wastewater_strength,
    )

__all__ = ('create_sand_filter_validation_system',)

ID = 'SandFilter_validation'


# %% Standalone filtration topology (feed profiles are shared with TF)
def create_sand_filter_validation_system(
        flowsheet=None, wastewater_strength='default',
        model='constant_splitting', sf_kwargs=None):
    '''
    Create a standalone raw-wastewater sand-filtration validation system.

    Parameters
    ----------
    flowsheet : qsdsan.Flowsheet, optional
        Flowsheet to populate. A new validation flowsheet is created by default.
    wastewater_strength : {'default', 'low', 'mid', 'high'}, optional
        Influent profile shared with TricklingFilter_validation. Nondefault
        profiles use its organic fractions, strength values, and Table S1
        component N/P and COD-to-mass properties.
    model : {'constant_splitting', 'empirical_fouling'}, optional
        Sand-filtration model. The default is ``'constant_splitting'``.
        Accepted aliases are those supported by :class:`SandFiltration`.
        A ``model`` entry in ``sf_kwargs`` takes precedence.
    sf_kwargs : dict, optional
        Additional :class:`SandFiltration` arguments, including split or
        fouling parameters and ``energy_consumption`` [kWh/m3].

    Returns
    -------
    qsdsan.System
        For constant splitting, the unit has RWW as its inlet and ``effluent``
        and ``backwash`` as its outlets. For empirical fouling, the
        ``filtered_water_recycle`` outlet is connected to the backwash inlet,
        and ``backwash`` remains the solids-bearing discharge.
    '''
    flowsheet = flowsheet or qs.Flowsheet(ID)
    qs.main_flowsheet.set_flowsheet(flowsheet)

    wastewater_strength = _normalize_wastewater_strength(wastewater_strength)
    _create_validation_components(wastewater_strength)
    rww = _create_validation_wastewater(wastewater_strength)

    sf_kwargs = dict(sf_kwargs or ())
    model = sf_kwargs.pop('model', model)
    model = SandFiltration._normalize_model(model)

    if model == 'constant_splitting':
        SF = SandFiltration(
            'SF', ins=rww, outs=('effluent', 'backwash'),
            model=model, **sf_kwargs)
        system = qs.System(ID, path=(SF,))
    else:
        recycle = qs.WasteStream('filtered_water_recycle')
        SF = SandFiltration(
            'SF', ins=(rww, recycle),
            outs=('effluent', 'backwash', recycle),
            model=model, **sf_kwargs)
        system = qs.System(ID, path=(SF,), recycle=recycle)

    if SF.isdynamic:
        system.set_dynamic_tracker(rww, SF, *SF.outs)
    return system


# %% Interactive construction only
if __name__ == '__main__':
    sys = create_sand_filter_validation_system()
    RWW = sys.flowsheet.stream.RWW
    SF = sys.flowsheet.unit.SF
    effluent = sys.flowsheet.stream.effluent
