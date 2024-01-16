# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References
----------
.. [1] Seider, W. D., Lewin, D. R., Seader, J. D., Widagdo, S., Gani, R. & Ng, K. M. (2017). 
    Product and Process Design Principles: Synthesis, Analysis and Evaluation (4th ed.). 
    John Wiley; Sons Inc. pp 441, 458-459.

"""
import biosteam as bst
from biosteam import Unit
from biosteam.utils import list_available_names
from biosteam.units.design_tools.mechanical import (
    motor_efficiency,
    nearest_NEMA_motor_size,
    electric_motor_cost
)
from math import log, exp, ceil
from warnings import warn
from typing import NamedTuple, Tuple, Callable


class BlowerCostAlgorithm(NamedTuple): 
    #: Defines preliminary correlation algorithm for a compressor type
    compression_ratio_bounds: Tuple[float, float] #: Compression ratio (recommended range of operation, to autodetermine blower type and/or issue warning)
    icfm_bounds: Tuple[float, float] #: Actual cubic feet per minute (hard limit for parallel units)
    hp_bounds: Tuple[float, float] #: Horse power per machine (not a hard limit for costing, but included here for completion)
    cost: Callable #: function(horse_power) -> Baseline purchase cost
    efficiencies: Tuple[float, float] #: Heuristic efficiencies within the range of compression ratios.
    CE: float #: Chemical engineering price cost index.


class Blower(Unit, isabstract=True):
    """
    Abstract class for blowers that includes design and costing. Child classes
    should implement the `_run` method for mass and energy balances. Preliminary 
    design and costing is estimated according to [1]_.
    
    """
    _N_ins = 1
    _N_outs = 1
    
    _specific_heat_ratio = 1.4   # air
    
    _units = {
        'Brake horsepower': 'hp',
        'Power': 'hp',
        'Motor size': 'hp',
    }

    material_factors = {
        'Carbon steel': 1.0,
        'Stainless steel': 2.5,
        # 'Nickel alloy': 5.0,
        'Cast aluminum': 0.6,
    }
    _F_BM_default = {
        'Blower(s)': 2.15,
    }

     #: dict[str, BlowerCostAlgorithm] Cost algorithms by compressor type.
    baseline_cost_algorithms = {
        'Centrifugal': BlowerCostAlgorithm(
                compression_ratio_bounds=(1.1, 2),
                icfm_bounds=(100, 50000),
                hp_bounds=(5, 1000),
                cost=lambda Pc: exp(7.0187 + 0.7900 * log(Pc)),
                efficiencies=(0.70, 0.80),
                CE=567,
            ),
        'Rotary straight-lobe': BlowerCostAlgorithm(
                compression_ratio_bounds=(1.2, 1.3),
                icfm_bounds=(20, 50000),
                hp_bounds=(1, 1000),
                cost=lambda Pc: exp(7.71751 + 0.7932 * log(Pc) - 0.0129 * log(Pc)**2),
                efficiencies=(0.50, 0.70),
                CE=567,
            ),
    }

    def _init(self, P, eta=None, blower_type=None, material=None, dP_design=30000):
        self.P = P  #: Outlet pressure [Pa]
        self.eta = eta  #: Blower efficiency, if None given, determined based on compressor type
        self.blower_type = 'Default' if blower_type is None else blower_type
        self.material = 'Carbon steel' if material is None else material
        self.dP_design = dP_design
    
    @property
    def eta(self):
        return self._eta
    @eta.setter
    def eta(self, eta):
        """[float] Blower mechanical efficiency. If None, the efficiency 
        will be determined based on blower type."""
        if eta is not None and (eta <= 0 or eta > 1):
            raise ValueError(f'blower efficiency must be within  (0, 1] or None, not {eta}.')
        self._eta = eta
    
    @property
    def blower_type(self):
        return self._blower_type
    @blower_type.setter
    def blower_type(self, blower_type):
        """[str] Type of blower. If 'Default', the type will be determined based on the outlet pressure."""
        blower_type = blower_type.capitalize()
        if blower_type not in self.baseline_cost_algorithms and blower_type != 'Default':
            raise ValueError(
                f"blower type {repr(blower_type)} not available; "
                f"only {list_available_names(self.baseline_cost_algorithms)} are available"
                )
        self._blower_type = blower_type
        
    @property
    def material(self):
        """[str]. Construction material. Defaults to 'Carbon steel'."""
        return self._material
    @material.setter
    def material(self, material):
        try:
            self.F_M['Blower(s)'] = self.material_factors[material]
        except KeyError:
            raise AttributeError("material must be one of the following: "
                                 f"{list_available_names(self.material_factors)}")
        self._material = material

    @property
    def compression_ratio(self):
        Pi = self.ins[0].P
        rc = self.P / Pi
        if rc < 1: 
            warn('compression ratio specified by inlet and outlet pressures is'
                 f'{rc:.2f}, invalid for blowers. Defaulting pressure difference to'
                 f'{self.dP_design} Pa.')
            rc = (Pi + self.dP_design) / Pi
        return rc

    def _determine_blower_type(self):
        cost_algorithms = self.baseline_cost_algorithms
        rc = self.compression_ratio
        icfm = self.ins[0].get_total_flow('cfm')
        N = 1
        if icfm > 50000: 
            N = ceil(icfm/50000)
            icfm /= N
        for name, alg in cost_algorithms.items():
            if icfm >= alg.icfm_bounds[0]: 
                rc_lb, rc_ub = alg.compression_ratio_bounds
                if rc_lb <= rc <= rc_ub: 
                    return name
        name = 'Centrifugal'
        warn('no blower available that is recommended for a compression ratio of '
            f'{rc:.2f} and an inlet cubic feet per minute of {icfm:.2f}; '
            f'defaulting to {name.lower()} blower')
        return name

    def _design(self):
        design_results = self.design_results
        blower_type = self.blower_type 
        if blower_type == 'Default': blower_type = self._determine_blower_type()
        design_results['Type'] = blower_type
        alg = self.baseline_cost_algorithms[blower_type]
        icfm_lb, icfm_ub = alg.icfm_bounds
        icfm = self.ins[0].get_total_flow('cfm')
        N = ceil(icfm / icfm_ub) if icfm > icfm_ub else 1
        
        k = self._specific_heat_ratio
        Pi = self.ins[0].P
        eta_B = self.eta or sum(alg.efficiencies)/2
        rc = self.compression_ratio
        design_results['Brake horsepower'] = Pb = 0.00436 * k/(k-1) * (icfm/N*Pi/eta_B) * (rc**(1-1/k)-1) # braker horsepower per blower
        eta_M = motor_efficiency(Pb)
        design_results['Power'] = Pc = Pb / eta_M
        design_results['Motor size'] = motor_rating = nearest_NEMA_motor_size(Pc)
        self.add_power_utility(N * Pc / 1.341) # Add power for blowers in kW
        design_results['Blowers in parallel'] = N 
        design_results['Motors in parallel'] = N * ceil(Pc / motor_rating) # Multiple motors per pump is posible
        
    def _cost(self):
        D = self.design_results
        C = self.baseline_purchase_costs
        alg = self.baseline_cost_algorithms[D['Type']]
        Pc = D['Power']
        Pm = D['Motor size']
        N = D['Blowers in parallel']
        Nm = D['Motors in parallel']
        C['Blower(s)'] = N * bst.CE/ alg.CE * alg.cost(Pc)
        C['Motor(s)'] = Nm * electric_motor_cost(Pm)
        self.F_M['Blower(s)'] = self.material_factors[self.material]