# -*- coding: utf-8 -*-
"""
Created on Thu Oct  9 17:45:54 2025

@author: aliah
"""

"""
Esterification + Washing + MeOH Recovery Train
==============================================

Hydrogen-free upgrading route for mid-cut (230–280 °C) HTL distillates
rich in fatty acids (C16–C18).  The sequence includes:

1. EsterificationUnit  (acid-catalyzed: FA + MeOH → FAME + H2O)
2. Acid quench and base wash
3. Water wash and drying (organic flash)
4. Methanol recovery and recycle loop

Designed for smooth integration in BioSTEAM/QSDsan systems.
"""

# from __future__ import annotations
import biosteam as bst
from biosteam.units.decorators import cost
from thermosteam import Reaction, ReactionSet

try:
    from qsdsan import System
except Exception:
    System = None


# --------------------------------------------------------------------------------------
# 1) Core esterification reactor
# --------------------------------------------------------------------------------------
@cost(basis='Wet mass flowrate', ID='Esterification reactor',
          units='kg/hr', CE=733.5, cost=35000, S=1000, n=0.6, BM=3.3)
class EsterificationUnit(bst.units.Reactor):
    """
    Acid-catalyzed esterification of free fatty acids (C16:0, C16:1, C18)
    with methanol to produce methyl esters (FAME) and water.
    """

    _N_ins = 3   # (distillate feed, MeOH, acid catalyst)
    _N_outs = 1

    def __init__(self, ID='', ins=(), outs=(), T=323.15, P=101325,
                 tau=2.0, conv_C16s=0.95, conv_C16u=0.95, conv_C18=0.95):
        super().__init__(ID, ins, outs, T=T, P=P, tau=tau)
        self.conv_C16s = conv_C16s
        self.conv_C16u = conv_C16u
        self.conv_C18  = conv_C18
        self._rxns = None

    def _setup_reactions(self):
        rxns = [
            Reaction('C16:0FA + MeOH -> C17FAME + H2O',  'C16:0FA',  self.conv_C16s),
            Reaction('C16:1FA + MeOH -> C17UFAME + H2O', 'C16:1FA', self.conv_C16u),
            Reaction('C18FACID + MeOH -> C19FAME + H2O', 'C18FACID', self.conv_C18),
        ]
        self._rxns = ReactionSet(rxns)

    def _run(self):
        feed, meoh, acid = self.ins
        out = self.outs[0]
        out.mix_from(self.ins)
        if self._rxns is None:
            self._setup_reactions()
        self._rxns.force_reaction(out.mol)
        out.T, out.P = self.T, self.P

    def _design(self):
        Fm = self.ins[0].F_mass
        self.design_results['Feed mass flow'] = Fm
        V = sum(s.F_vol for s in self.ins) * self.tau * 3600.0
        self.design_results['Reactor volume (m3)'] = max(V, 0.1)
        self.add_power_utility(0.2)  # agitator load


# --------------------------------------------------------------------------------------
# 2) Supporting wash & recovery units (simple placeholders)
# --------------------------------------------------------------------------------------
class AcidQuench(bst.Unit):
    """Neutralize residual acid catalyst with Na2CO3 (simplified)."""
    _N_ins = 2
    _N_outs = 1

    def _run(self):
        org, base = self.ins
        out = self.outs[0]
        out.mix_from(self.ins)
        out.T = org.T


class MixerSettlerWash(bst.units.MixerSettler):
    """Generic mixer–settler stage for acid/base or water washing."""
    pass


class OrganicDryer(bst.units.Flash):
    """Mild flash to strip MeOH/H2O from organic product."""
    pass


class RecoveryFlash(bst.units.Flash):
    """Flash to separate MeOH-rich vapor and aqueous bottoms."""
    pass


class MeOHRecycleSplitter(bst.units.Splitter):
    """Split condensed MeOH to recycle and purge."""
    pass


# --------------------------------------------------------------------------------------
# 3) Composite factory that wires everything together
# --------------------------------------------------------------------------------------
def build_esterification_wash_recovery(
    ID_prefix='EST',
    distillate=None,
    meoh_fresh=None,
    acid_catalyst=None,
    wash_water=None,
    base_solution=None,
    recycle_meoh_target=0.90,
    dryer_T=338.15,
    rec_T=351.15,
    rec_P=101325,
):
    """
    Build the esterification + washing + MeOH recycle sequence.

    Returns
    -------
    units : dict
    streams : dict
    sys : bst.System
    _qsdsys : optional QSDsan.System
    """
    # Mix fresh and recycled MeOH
    meoh_recycle = bst.Stream(f'{ID_prefix}_meoh_recycle')
    meoh_mix = bst.Mixer(f'{ID_prefix}_meoh_mix', ins=(meoh_fresh, meoh_recycle))

    # Reactor
    R_est = EsterificationUnit(f'{ID_prefix}_R',
                               ins=(distillate, meoh_mix-0, acid_catalyst))

    # Quench & washes
    Q = AcidQuench(f'{ID_prefix}_quench', ins=(R_est-0, base_solution))
    MS1 = MixerSettlerWash(f'{ID_prefix}_MS1', ins=(Q-0, wash_water))
    org1, aq1 = MS1.outs
    MS2 = MixerSettlerWash(f'{ID_prefix}_MS2', ins=(org1, wash_water.copy()))
    org2, aq2 = MS2.outs

    # Drying and recovery
    Dryer = OrganicDryer(f'{ID_prefix}_Dry', ins=org2, T=dryer_T, P=101325)
    org_vap, org_liq = Dryer.outs
    aq_pool = bst.Mixer(f'{ID_prefix}_aq_pool', ins=(aq1, aq2, org_vap))
    Rec = RecoveryFlash(f'{ID_prefix}_Rec', ins=aq_pool-0, T=rec_T, P=rec_P)
    meoh_vap, water_btms = Rec.outs
    Cond = bst.units.Condenser(f'{ID_prefix}_Cond', ins=meoh_vap)
    Split = MeOHRecycleSplitter(f'{ID_prefix}_Split', ins=Cond-0,
                                split=recycle_meoh_target)
    meoh_recycle_stream, meoh_purge = Split.outs
    meoh_recycle.copy_like(meoh_recycle_stream)

    product_fuel = org_liq
    wastewater = water_btms

    units = dict(R_est=R_est, Q=Q, MS1=MS1, MS2=MS2, Dryer=Dryer,
                 aq_pool=aq_pool, Rec=Rec, Cond=Cond, Split=Split)
    streams = dict(product_fuel=product_fuel, meoh_recycle=meoh_recycle,
                   meoh_purge=meoh_purge, wastewater=wastewater)

    sys = bst.System(f'{ID_prefix}_sys',
                     path=[meoh_mix, R_est, Q, MS1, MS2, Dryer,
                           aq_pool, Rec, Cond, Split])

    _qsdsys = System(ID=f'{ID_prefix}_qsdsan', system=bst.main_flowsheet.system) if System else None
    return units, streams, sys, _qsdsys


# --------------------------------------------------------------------------------------
# 4) Helper: one-call builder for easy system inclusion
# --------------------------------------------------------------------------------------
def add_esterification_section(
    upstream_stream,
    wash_water,
    base_solution,
    meoh_fresh,
    acid_catalyst=None,
    ID_prefix='EST',
    recycle_meoh_target=0.9,
):
    """
    High-level wrapper to attach the esterification train to an existing system.

    Example:
        EST_section = add_esterification_section(BiofuelFlash-1, wash, base, MeOH)
    """
    units, streams, sys, _ = build_esterification_wash_recovery(
        ID_prefix=ID_prefix,
        distillate=upstream_stream,
        meoh_fresh=meoh_fresh,
        acid_catalyst=acid_catalyst,
        wash_water=wash_water,
        base_solution=base_solution,
        recycle_meoh_target=recycle_meoh_target,
    )
    return dict(units=units, streams=streams, sys=sys)


# --------------------------------------------------------------------------------------
# 5) Minimal local test
# --------------------------------------------------------------------------------------
if __name__ == '__main__':
    dist = bst.Stream('dist', **{'C16:0FA': 18, 'C16:1FA': 206, 'C18FACID': 24},
                      units='kg/hr', T=323.15)
    meoh = bst.Stream('meoh_fresh', MeOH=50, units='kg/hr')
    acid = bst.Stream('acid_cat', H2SO4=1, units='kg/hr')
    wash = bst.Stream('wash', H2O=200, units='kg/hr')
    base = bst.Stream('base', Na2CO3=2, H2O=50, units='kg/hr')

    u, s, sys, _ = build_esterification_wash_recovery(
        distillate=dist, meoh_fresh=meoh,
        acid_catalyst=acid, wash_water=wash, base_solution=base)

    # simple placeholder splits
    for unit in (u['MS1'], u['MS2']):
        unit.split = 0.2
    u['Dryer'].V = 0.05
    u['Rec'].V = 0.6

    sys.simulate()
    print('\n--- Product fuel ---')
    print(s['product_fuel'].show())
    print('\n--- Wastewater ---')
    print(s['wastewater'].show())
