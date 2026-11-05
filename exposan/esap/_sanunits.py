#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Project-specific unit operations for the ESAP fertilizer systems.
'''

from math import ceil
from warnings import warn

import biosteam as bst
from biosteam.exceptions import CostWarning
from biosteam.units import BatchCrystallizer
from biosteam.units.design_tools import size_batch

__all__ = ('ESAPBatchCrystallizer',)

_m3_to_gal = 1 / 0.003785411784
_kg_to_short_ton = 0.00110231131


class ESAPBatchCrystallizer(BatchCrystallizer):
    '''
    Black-box batch crystallizer with wet-cake separation.

    This project-specific child class adds a material-balance model to
    :class:`biosteam.units.BatchCrystallizer`. A specified fraction of a
    dissolved salt is converted to an equal molar amount of a solid crystal
    component. A specified fraction of those crystals is captured in a wet
    cake; uncrystallized salt, uncaptured crystals, and unused water leave in
    the mother liquor.

    Parameters
    ----------
    ins : stream
        Concentrated salt solution.
    outs : iterable[stream]
        Wet crystal cake and mother liquor, respectively.
    solute_ID : str
        Component ID for dissolved salt.
    crystal_ID : str
        Component ID for the corresponding solid crystal. The dissolved and
        crystal components must have equal molecular weights.
    crystallization_recovery : float, optional
        Fraction of dissolved salt converted to crystals. Default is 0.85.
    crystal_capture_efficiency : float, optional
        Fraction of formed crystals captured in the wet cake. Default is 0.98.
    wet_cake_moisture : float, optional
        Water mass fraction of the wet cake. Default is 0.20. Dissolved salt is
        conservatively assigned to the mother liquor rather than entrained in
        the cake.
    separator_type : str, optional
        Currently only ``'reciprocating_pusher'`` is supported. The purchase
        cost correlation is from BioSTEAM's ``SolidsCentrifuge``: ``68040 *
        loading**0.5`` at CEPCI 567, with loading in short ton crystal/hr and a
        bare-module factor of 2.03. Its published range is 1--20 short ton/hr;
        use below that range is an explicit extrapolation. Separator count is
        ``N_sep = ceil(loading / separator_capacity)``,
        with zero separators for zero loading. Separators share the loading
        equally and are independent of the number of crystallizer vessels.
        Total baseline purchase cost [USD] is
        ``N_sep * 68040 * (loading / N_sep)**0.5 * bst.CE / 567``.
        The bare-module factor is applied separately to obtain installed cost.
        This is an average-throughput sizing assumption, not a check of peak
        batch discharge capacity. Source: ``SolidsCentrifuge`` in BioSTEAM's
        ``units/solids_separation.py``.
    heat_of_crystallization : float, optional
        Heat released per kg of salt crystallized, [kJ/kg]. Default is zero
        because product-specific values have not yet been selected. Positive
        values increase the required cooling duty.
    separator_kWhr_per_m3 : float, optional
        Separator electricity intensity, [kWh/m3 feed]. Default is 1.40,
        matching BioSTEAM's ``SolidsCentrifuge`` default.
    separator_capacity : float, optional
        Maximum captured-crystal loading per separator, [short ton/hr]. Must
        be positive. Default is 20, the upper loading bound of BioSTEAM's
        reciprocating-pusher centrifuge cost correlation. The separator count
        is the total captured-crystal loading divided by this capacity,
        rounded up; zero loading requires zero separators. Values above 20
        may extrapolate the cost correlation beyond its published range.
    tau : float, optional
        Crystallization residence time, [hr]. Passed to the parent class.
    N : int, optional
        Number of crystallizer vessels. Passed to the parent class.
    V : float, optional
        Maximum vessel volume, [m3]. Specify either ``N`` or ``V``.
    T : float, optional
        Crystallization temperature, [K]. Default is 305.15 K.
    kW : float, optional
        Agitation power intensity, [kW/gal working volume]. The parent default
        is 0.00746 kW/gal; this class explicitly converts m3 to US gallons.

    Notes
    -----
    The model is intentionally non-equilibrium. Recovery, capture efficiency,
    and cake moisture must be supported by experiments, literature, or
    uncertainty analysis. Heat duty includes sensible cooling of the feed to
    ``T`` plus the specified heat of crystallization. The default heat of
    crystallization is zero.
    '''

    _N_ins = 1
    _N_outs = 2
    _units = {
        **BatchCrystallizer._units,
        'Crystallization recovery': '%',
        'Crystal capture efficiency': '%',
        'Crystal production': 'kg/hr',
        'Captured crystals': 'kg/hr',
        'Separator solids loading': 'ton/hr',
        'Number of separators': '',
        'Agitator power': 'kW',
        'Separator power': 'kW',
        }
    _F_BM_default = {
        **BatchCrystallizer._F_BM_default,
        'Crystal separator': 2.03,
        }

    def _init(self, solute_ID, crystal_ID, crystallization_recovery=0.85,
              crystal_capture_efficiency=0.98, wet_cake_moisture=0.20,
              separator_type='reciprocating_pusher',
              heat_of_crystallization=0, separator_kWhr_per_m3=1.40,
              separator_capacity=20,
              **kwargs):
        super()._init(**kwargs)
        self.solute_ID = solute_ID
        self.crystal_ID = crystal_ID
        self.crystallization_recovery = crystallization_recovery
        self.crystal_capture_efficiency = crystal_capture_efficiency
        self.wet_cake_moisture = wet_cake_moisture
        if heat_of_crystallization < 0:
            raise ValueError('`heat_of_crystallization` cannot be negative.')
        self.heat_of_crystallization = heat_of_crystallization
        if separator_kWhr_per_m3 < 0:
            raise ValueError('`separator_kWhr_per_m3` cannot be negative.')
        self.separator_kWhr_per_m3 = separator_kWhr_per_m3
        self.separator_capacity = separator_capacity
        if separator_type != 'reciprocating_pusher':
            raise ValueError(
                '`separator_type` must currently be "reciprocating_pusher".'
                )
        self.separator_type = separator_type
        dissolved = self.chemicals[solute_ID]
        crystal = self.chemicals[crystal_ID]
        if abs(dissolved.MW - crystal.MW) > 1e-9:
            raise ValueError(
                'Dissolved and crystal components must have equal molecular '
                'weights.'
                )

    @staticmethod
    def _validate_fraction(name, value, upper_open=False):
        valid = 0 <= value < 1 if upper_open else 0 <= value <= 1
        if not valid:
            upper = 'less than 1' if upper_open else 'at most 1'
            raise ValueError(f'`{name}` must be nonnegative and {upper}.')

    @property
    def separator_capacity(self):
        '''[float] Positive captured-crystal capacity per separator [short ton/hr].'''
        return self._separator_capacity

    @separator_capacity.setter
    def separator_capacity(self, value):
        if not 0 < value < float('inf'):
            raise ValueError('`separator_capacity` must be positive and finite.')
        self._separator_capacity = value

    @property
    def crystallization_recovery(self):
        return self._crystallization_recovery

    @crystallization_recovery.setter
    def crystallization_recovery(self, value):
        self._validate_fraction('crystallization_recovery', value)
        self._crystallization_recovery = value

    @property
    def crystal_capture_efficiency(self):
        return self._crystal_capture_efficiency

    @crystal_capture_efficiency.setter
    def crystal_capture_efficiency(self, value):
        self._validate_fraction('crystal_capture_efficiency', value)
        self._crystal_capture_efficiency = value

    @property
    def wet_cake_moisture(self):
        return self._wet_cake_moisture

    @wet_cake_moisture.setter
    def wet_cake_moisture(self, value):
        self._validate_fraction('wet_cake_moisture', value, upper_open=True)
        self._wet_cake_moisture = value

    def _setup(self):
        bst.Unit._setup(self)
        for stream in self.outs:
            stream.T = self.T
            stream.P = 101325

    def _run(self):
        feed = self.ins[0]
        wet_cake, mother_liquor = self.outs
        wet_cake.empty()
        mother_liquor.copy_like(feed)

        solute_ID = self.solute_ID
        crystal_ID = self.crystal_ID
        crystallized = feed.imol[solute_ID] * self.crystallization_recovery
        captured = crystallized * self.crystal_capture_efficiency

        mother_liquor.imol[solute_ID] -= crystallized
        mother_liquor.imol[crystal_ID] = crystallized - captured
        wet_cake.imol[crystal_ID] = captured

        crystal_mass = wet_cake.imass[crystal_ID]
        target_water = (
            crystal_mass * self.wet_cake_moisture
            / (1 - self.wet_cake_moisture)
            )
        available_water = mother_liquor.imass['H2O']
        cake_water = min(target_water, available_water)
        wet_cake.imass['H2O'] = cake_water
        mother_liquor.imass['H2O'] -= cake_water
        if cake_water + 1e-12 < target_water:
            warn(
                f'{self.ID} has insufficient water to meet the requested wet '
                'cake moisture.',
                RuntimeWarning,
                )

        wet_cake.phase = mother_liquor.phase = 'l'
        wet_cake.T = mother_liquor.T = self.T
        wet_cake.P = mother_liquor.P = 101325

    def _design(self):
        flow = sum(stream.F_vol for stream in self.outs)
        tau = self.tau
        tau_0 = self.tau_0
        working_fraction = self.V_wf
        if self.V:
            design = size_batch(
                flow, tau, tau_0, working_fraction, V_max=self.V,
                )
        else:
            design = size_batch(
                flow, tau, tau_0, working_fraction, N_reactors=self.N,
                )

        results = self.design_results
        results['Crystallizer volume'] = volume = design.pop('Reactor volume')
        results.update(design)
        results['Number of crystallizers'] = N = design['Number of reactors']
        results['Crystallization recovery'] = 100 * self.crystallization_recovery
        results['Crystal capture efficiency'] = 100 * self.crystal_capture_efficiency
        results['Crystal production'] = (
            self.ins[0].imol[self.solute_ID]
            * self.crystallization_recovery
            * self.chemicals[self.solute_ID].MW
            )
        results['Captured crystals'] = self.outs[0].imass[self.crystal_ID]

        working_volume_gal = working_fraction * volume * N * _m3_to_gal
        agitator_power = self.kW * working_volume_gal
        separator_power = self.separator_kWhr_per_m3 * flow
        results['Agitator power'] = agitator_power
        results['Separator power'] = separator_power
        self.add_power_utility(agitator_power + separator_power)
        feed = self.ins[0]
        cooled_feed = feed.copy(None)
        cooled_feed.T = self.T
        cooled_feed.P = 101325
        sensible_duty = cooled_feed.H - feed.H
        crystallization_duty = (
            -self.heat_of_crystallization * results['Crystal production']
            )
        self.add_heat_utility(sensible_duty + crystallization_duty, self.T)

        loading = results['Captured crystals'] * _kg_to_short_ton
        results['Separator solids loading'] = loading
        results['Number of separators'] = (
            max(1, ceil(loading / self.separator_capacity)) if loading else 0
            )
        # Parallel machines share average captured-crystal loading,
        # independently of vessel count.
        
    def _cost(self):
        super()._cost()
        loading = self.design_results['Separator solids loading']
        if not loading:
            self.baseline_purchase_costs['Crystal separator'] = 0
            return
        if loading < 1:
            warn(
                f'{self.ID} crystal separator solids loading '
                f'({loading:.4g} ton/hr) is below the 1 ton/hr lower bound '
                'of the reciprocating-pusher centrifuge cost correlation.',
                CostWarning,
                )
        N = self.design_results['Number of separators']
        loading_per_separator = loading / N
        # BioSTEAM SolidsCentrifuge reciprocating-pusher purchase cost:
        # 68040 * L**0.5 USD per machine at CEPCI 567, L in short ton solids/hr.
        # Escalate to bst.CE and sum N equally loaded parallel machines.
        # The 2.03 bare-module factor is applied separately by unit costing.
        self.baseline_purchase_costs['Crystal separator'] = (
            N * 68040 * loading_per_separator**0.5 * bst.CE / 567
            )
