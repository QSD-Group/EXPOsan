# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

Energy accounting for the ESAP trickling-filter validation system.

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import math

from qsdsan.utils import auom
from exposan.werf._units import TricklingFilter

__all__ = (
    'add_trickling_filter_energy_metrics',
    'get_influent_pumping_power',
    'get_trickling_filter_recirculation_pumping_power',
    'get_trickling_filter_rotating_arm_power',
    )


# %% Energy-accounting assumptions (not biological calibration)
TF_rotating_arm_power = 3e-3             # kW/m3 media volume
TF_distribution_headloss = 5.5           # ft
TF_influent_site_lift = 50.              # ft, site-level lift
TF_influent_headloss = 3.28              # ft, influent pump loss
TF_influent_pump_efficiency = 0.60       # -
TF_recirculation_pump_headloss = 3.28    # ft, recycle pump loss
TF_recirculation_pump_efficiency = 0.60  # -


# %% Hydraulic validation and power calculations
def _validate_efficiency(efficiency):
    """Return a finite hydraulic efficiency in the interval (0, 1]."""
    efficiency = float(efficiency)
    if not math.isfinite(efficiency) or not 0. < efficiency <= 1.:
        raise ValueError('efficiency must be finite and in the interval (0, 1]')
    return efficiency


def get_trickling_filter_rotating_arm_power(unit):
    """Return validation rotary-distributor power in kW.

    The 0.003 kW/m3 media-volume demand follows Tarallo et al. (2015).
    """
    if not isinstance(unit, TricklingFilter):
        raise TypeError('unit must be a TricklingFilter')
    media_volume = unit.filter_bed_surface * unit.filter_bed_depth
    return media_volume * TF_rotating_arm_power


def get_influent_pumping_power(
        stream, site_lift_ft=TF_influent_site_lift,
        headloss_ft=TF_influent_headloss,
        efficiency=TF_influent_pump_efficiency):
    """Return ESAP validation influent pumping power in kW.

    Power is based on raw influent flow, a 50-ft site lift plus 3.28-ft head
    loss by default, water density of 1000 kg/m3, and 60% efficiency.
    """
    efficiency = _validate_efficiency(efficiency)
    if stream.F_vol <= 0.:
        return 0.
    head_m = (site_lift_ft + headloss_ft) \
        * auom('ft').conversion_factor('m')
    flow_m3_s = stream.F_vol / 3600.
    return 1000. * 9.81 * flow_m3_s * head_m / efficiency / 1000.


def get_trickling_filter_recirculation_pumping_power(
        unit, efficiency=TF_recirculation_pump_efficiency):
    """Return ESAP validation recycle pumping power in kW.

    Recycle flow is ``Q_R = R * Q_in``. Head includes bed depth, the WERF
    5.5-ft distribution/collection loss, and a 3.28-ft recycle-pump loss.
    Hydraulic efficiency defaults to 60%.
    """
    if not isinstance(unit, TricklingFilter):
        raise TypeError('unit must be a TricklingFilter')
    feed = unit.ins[0]
    ratio = unit.recirculation_ratio
    if feed.F_vol <= 0. or ratio <= 0.:
        return 0.
    efficiency = _validate_efficiency(efficiency)
    head_m = unit.filter_bed_depth + (
        TF_distribution_headloss + TF_recirculation_pump_headloss
        ) * auom('ft').conversion_factor('m')
    flow_m3_s = feed.F_vol * ratio / 3600.
    return 1000. * 9.81 * flow_m3_s * head_m / efficiency / 1000.


# %% Model metric registration (reads the simulated state; no new solve)
def add_trickling_filter_energy_metrics(model, include_rotating_arm=True,
                                        include_cost=True):
    """Add ESAP validation trickling-filter electricity metrics to a model."""
    metric = model.metric
    system = model.system
    units = tuple(i for i in system.units if isinstance(i, TricklingFilter))
    if not units:
        raise ValueError('model system does not contain a TricklingFilter')

    @metric(name='influent pumping energy', units='kW', element='Pumping')
    def get_influent_pump_power():
        return get_influent_pumping_power(system.flowsheet.stream.RWW)

    @metric(name='recirculation pumping energy', units='kW', element='Pumping')
    def get_recirculation_pump_power():
        return sum(get_trickling_filter_recirculation_pumping_power(i)
                   for i in units)

    if include_rotating_arm:
        @metric(name='rotating arm energy consumption', units='kW',
                element='Pumping')
        def get_tf_rotating_arm_power():
            return sum(get_trickling_filter_rotating_arm_power(i) for i in units)

    if include_cost:
        @metric(name='trickling filter electricity cost', units='USD/d',
                element='OPEX')
        def get_tf_energy_cost():
            power = get_influent_pump_power() + get_recirculation_pump_power()
            if include_rotating_arm:
                power += get_tf_rotating_arm_power()
            return power * 24. * system.power_utility.price
