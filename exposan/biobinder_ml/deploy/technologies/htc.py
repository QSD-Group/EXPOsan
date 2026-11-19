# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Hydrothermal carbonization with hydrochar matrix, ash and moisture reporting."""
from dataclasses import dataclass, field
from ._thermal_basis import ThermalConfig


@dataclass(frozen=True)
class HTCConfig(ThermalConfig):
    # https://www.mdpi.com/1996-1073/10/2/211 : 220 C, 1 h process design.
    temperature_k: float = 493.15
    # https://www.mdpi.com/1996-1073/10/2/211 AppendixA.3: Psat +10 bar.
    # 3.4 MPa rounds the 220 C saturation pressure plus gas-head allowance.
    pressure_pa: float = 3.4e6
    residence_time_h: float = 1.
    slurry_solids: float = .20  # specified slurry concentration, not grape-marc calibration.
    oil_retention: dict = field(default_factory=dict)
    # https://www.mdpi.com/1996-1073/16/5/2483 Table3: measured hydrochar is
    # CHONS+ash, not carbon. These retentions are uncalibrated scenarios; reported
    # HTC6 C/S retentions exceed feed availability and are not imposed as balances.
    char_retention: dict = field(default_factory=lambda:dict(C=.75,H=.55,O=.45,N=.65,S=.75))
    aqueous_retention: dict = field(default_factory=lambda:dict(C=.20,H=.40,O=.50,N=.30,S=.20))


def create_system(capacity_tpd, *, feedstock='food', config=None, economics=None,
                  characterization=None, markets=None, ID=None, **kwargs):
    """Return a ThermalResult with heat recovery, dewatering and final drying."""
    from ._thermal_system import create_system as build
    return build('htc',capacity_tpd,feedstock=feedstock,config=config or HTCConfig(),
        economics=economics,characterization=characterization,markets=markets,ID=ID,**kwargs)
