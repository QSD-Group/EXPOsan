# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Distinct methane and hydrogen supercritical-water gasification scenarios."""
from dataclasses import dataclass, field
from ._thermal_basis import ThermalConfig


@dataclass(frozen=True)
class HTGConfig(ThermalConfig):
    route: str = 'methane'
    temperature_k: float | None = None
    pressure_pa: float = 28e6
    # https://www.psi.ch/en/cpe/hydropilot : wet-feed pilot process basis.
    # 15 min and 10 wt% slurry are engineering scenarios, not fitted pilot kinetics.
    residence_time_h: float = .25
    slurry_solids: float = .10
    methanation_completion: float | None = None
    oil_retention: dict = field(default_factory=dict)
    char_retention: dict = field(default_factory=dict)
    aqueous_retention: dict = field(default_factory=dict)

    def __post_init__(self):
        if self.route not in ('methane','hydrogen'):
            raise ValueError('HTG route must be methane or hydrogen')
        if self.temperature_k is None:
            # https://www.psi.ch/en/cpe/hydropilot : methane 450 C/280 bar.
            # https://doi.org/10.1016/j.biombioe.2022.106508 : food-waste H2
            # experiments at650 C included14 wt% K2CO3. Only the temperature is
            # adopted as a scenario here; that catalyst/yield is NOT transferred.
            object.__setattr__(self,'temperature_k',723.15 if self.route == 'methane' else 923.15)
        if self.methanation_completion is None:
            # https://www.psi.ch/en/cpe/hydropilot : route distinction only.
            # Extents .97/.05 are explicit scenarios, not fitted equilibrium.
            object.__setattr__(self,'methanation_completion',.97 if self.route == 'methane' else .05)
        super().__post_init__()


def create_system(capacity_tpd, *, route=None, feedstock='food', config=None,
                  economics=None, characterization=None, markets=None, ID=None, **kwargs):
    """Return a native HTG ThermalResult for the selected gas-product route."""
    from ._thermal_system import create_system as build
    c = config or HTGConfig(route=route or 'methane')
    if route is not None and route != c.route:
        raise ValueError('route and HTGConfig.route disagree')
    return build('htg_'+c.route,capacity_tpd,feedstock=feedstock,config=c,economics=economics,
        characterization=characterization,markets=markets,ID=ID,**kwargs)
