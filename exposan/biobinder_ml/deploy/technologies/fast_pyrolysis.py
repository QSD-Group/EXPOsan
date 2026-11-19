# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Native fast-pyrolysis equipment network with explicit product-basis scenarios."""
from dataclasses import dataclass
from ._thermal_basis import ThermalConfig


@dataclass(frozen=True)
class FastPyrolysisConfig(ThermalConfig):
    # https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-23053.pdf
    # 500 C design condition; 2 s is a vapor-residence scenario, not solids holdup.
    residence_time_h: float = 2/3600


def create_system(capacity_tpd, *, feedstock='food', config=None, economics=None,
                  characterization=None, markets=None, ID=None, **kwargs):
    """Return a ThermalResult with native drying, conversion and recovery units."""
    from ._thermal_system import create_system as build
    return build('fast_pyrolysis',capacity_tpd,feedstock=feedstock,
        config=config or FastPyrolysisConfig(),economics=economics,
        characterization=characterization,markets=markets,ID=ID,**kwargs)
