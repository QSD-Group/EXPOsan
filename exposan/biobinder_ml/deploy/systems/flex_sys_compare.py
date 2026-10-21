# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References
[1] Snowden-Swan et al., Wet Waste Hydrothermal Liquefaction and Biocrude Upgrading to Hydrocarbon Fuels:
    2021 State of Technology; PNNL-32731; Pacific Northwest National Lab. (PNNL), Richland, WA (United States), 2022.
    https://doi.org/10.2172/1863608.
[2] Feng et al, Characterizing the Opportunity Space for Sustainable
    Hydrothermal Valorization of Wet Organic Wastes.
    Environ. Sci. Technol. 2024, 58 (5), 2528–2541.
    https://doi.org/10.1021/acs.est.3c07394.
[3] Nordahl et al., Life-Cycle Greenhouse Gas Emissions and Human Health Trade-Offs
    of Organic Waste Management Strategies.
    Environ. Sci. Technol. 2020, 54 (15), 9200–9209.
    https://doi.org/10.1021/acs.est.0c00364.
[4] Lopez-Ruiz et al., Electrocatalytic Valorization into H2 and Hydrocarbons
    of an Aqueous Stream Derived from Hydrothermal Liquefaction.
    J Appl Electrochem 2021, 51 (1), 107–118.
    https://doi.org/10.1007/s10800-020-01452-x.
[5] Lopez-Ruiz et al., Low-Temperature Electrochemical Wastewater Oxidation;
    PNNL-35535; 2023. https://doi.org/10.2172/2332860.


'''

# !!! Temporarily ignoring warnings
import warnings
# Keep simulation warnings visible for cut feasibility checks.

import os
import numpy as np
import biosteam as bst
import qsdsan as qs
from biosteam import IsenthalpicValve
from qsdsan import sanunits as qsu
from qsdsan.utils import clear_lca_registries
from exposan.htl import create_tea

# ============================
# UPSTREAM: use biobinder_ml
# ============================
from exposan.saf import (
    _HHV_per_GGE,
    _load_components,
    _load_process_settings,
    _units as u_saf,
    data_path,
    dry_flowrate,
    find_Lr_Hr,
    get_mass_energy_balance,
    gwp_dct,
    # HTL_yields,
    price_dct,
    results_path,
    tea_kwargs,
    uptime_ratio,
)

from exposan.biobinder_ml import (
    central_dry_flowrate as default_central,
    pilot_dry_flowrate as default_pilot,
    feedprice_dct as feedprice_dct,
    HTL_yields,
    limit_internal_hx,
    configure_utilities,
    set_eff_hx_temperature,
    BiobinderTEA,
)
from exposan.biobinder_ml._feedstocks import get_feedstock_composition
from exposan.biobinder_ml import _units as u_bb
from exposan.biobinder_ml.Dist_flex import gwp_dict
# ============================
# DOWNSTREAM: SAF helpers
# ============================
from exposan.saf import price_dct as _base_price_dct, gwp_dct

# Shared runtime key search and trial registry.
from exposan.biobinder_ml._distill_utils_flex import add_runtime_LHK_LrHr_spec, get_distillation_trials


_psi_to_Pa = 6894.76
_m3_to_gal = 264.172
SAF_ECON_MODE = "market"  # "market" or "mfsp"
configure_utilities()

__all__ = ("create_system", "get_MFSP", "simulate_and_print", "fuel_blend_ratios")

class _DummyEffHX:
    """
    Minimal object that satisfies exposan.saf._units.Hydroprocessing._design()
    as 'eff_hx'. It prevents thermo/chemical mismatch during mix_from(self.outs).
    """
    def __init__(self, ID, thermo):
        self.ID = ID
        # IMPORTANT: both ins and outs must be real streams (not None)
        self.ins = [bst.Stream(f"{ID}_in", thermo=thermo)]
        self.outs = [bst.Stream(f"{ID}_out", thermo=thermo)]
        # Hydroprocessing._design uses eff_hx.Hnet
        self.Hnet = 0.0

    def simulate_as_auxiliary_exchanger(self, ins, outs, duty=None):
        self._last_duty = duty
        return


_RF_MODELS = {}
_INTEGRATION_ASSUMPTIONS = {}


def _get_rf_model(model_path=None):
    """Load/cache the existing trained model; never silently use default yields."""
    from pathlib import Path
    from joblib import load
    if model_path is None:
        model_path = Path(r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\results\rf_yield_model.joblib")
    model_path = Path(model_path).expanduser().resolve()
    if not model_path.is_file():
        raise FileNotFoundError(
            f"RF model not found: {model_path}. Supply rf_model_path to create_system.")
    key = str(model_path)
    if key not in _RF_MODELS:
        print(f"Loading RF model: {model_path}")
        _RF_MODELS[key] = load(model_path)
    return _RF_MODELS[key]


def _predict_rf_yields(feedstock_composition, rf_model_path=None):
    """Reuse Dist_flex's feature mapping and predictor at the actual HTL baseline."""
    from exposan.biobinder_ml._rf_htl_predictor import predict_htl_yields_from_rf
    names = ("Carbohydrates", "Proteins", "Lipids", "Ash")
    weights = [float(feedstock_composition[name]) for name in names]
    if any(not np.isfinite(x) or x < 0 for x in weights) or sum(weights) <= 0:
        raise ValueError("RF feedstock dry-component masses must be finite, nonnegative, and have positive total")
    carb, protein, lipid, ash = [100*x/sum(weights) for x in weights]
    predicted = predict_htl_yields_from_rf(
        rf_model=_get_rf_model(rf_model_path),
        carb_wt=carb, protein_wt=protein, lipids_wt=lipid, ash_wt=ash,
        process={"Temperature (C)": 280, "Residence Time": 15,
                 "Solid content (w/w) %": 20})
    predicted = dict(predicted)
    if not predicted or any(not np.isfinite(float(x)) or float(x) < 0
                            for x in predicted.values()):
        raise ValueError("RF predictor returned empty, negative, or nonfinite HTL yields")
    print("RF-predicted HTL yields:", predicted)
    return predicted


class JointBiocrudeSplitter(u_bb.BiocrudeSplitter):
    """Local splitter with indexed bins, empty-bin validation and mass closure."""

    def _update_component_ratios(self):
        if not hasattr(self, '_biocrude_ratios') or not hasattr(self, '_cutoff_fracs'):
            return
        boundaries = tuple(self.cutoff_Tbs)
        fractions = tuple(self.cutoff_fracs)
        if len(fractions) != len(boundaries) + 1:
            raise ValueError("Each boiling boundary requires one additional fraction")
        if any(a >= b for a, b in zip(boundaries, boundaries[1:])):
            raise ValueError("Boiling boundaries must increase strictly")
        bins = [{} for _ in fractions]
        for ID, weight in self.biocrude_ratios.items():
            if weight <= 0:
                continue
            Tb = float(self.components[ID].Tb)
            if not np.isfinite(Tb):
                raise ValueError(f"Missing/nonfinite boiling point for {ID}")
            idx = int(np.searchsorted(boundaries, Tb, side="left"))
            bins[idx][ID] = float(weight)
        for idx, group in enumerate(bins):
            total = sum(group.values())
            if total <= 0:
                raise ValueError(f"Cut bin {idx} has no represented components; cannot assign fraction {fractions[idx]}")
            bins[idx] = {ID: weight/total for ID, weight in group.items()}
        self._joint_bins = tuple(bins)
        self._keys = [(max(bins[i], key=lambda ID: self.components[ID].Tb),
                       min(bins[i+1], key=lambda ID: self.components[ID].Tb))
                      for i in range(len(boundaries))]

    def _run(self):
        inlet, outlet = self.ins[0], self.outs[0]
        total = float(inlet.imass[self.biocrude_IDs].sum())
        outlet.copy_like(inlet)
        # Do not discard pre-existing represented compounds.
        for fraction, group in zip(self.cutoff_fracs, self._joint_bins):
            for ID, ratio in group.items():
                outlet.imass[ID] += total*fraction*ratio
        outlet.imass[self.biocrude_IDs] = 0
        if not np.isclose(outlet.F_mass, inlet.F_mass, rtol=1e-8, atol=1e-8):
            raise RuntimeError(f"Biocrude split mass imbalance: in={inlet.F_mass}, out={outlet.F_mass}")


def _joint_biocrude_ratios(intermediate_surrogate_ratios):
    """Preserve original bin compositions; represent empty intermediate explicitly."""
    import inspect
    base = inspect.signature(u_bb.BiocrudeSplitter.__init__).parameters['biocrude_ratios'].default
    ratios = dict(base)
    if intermediate_surrogate_ratios is None:
        raise ValueError("Supply intermediate_surrogate_ratios explicitly; e.g. {'C18H38': 1.0} is only a surrogate sensitivity, not an HTL assay")
    components = bst.settings.thermo.chemicals
    for ID, weight in intermediate_surrogate_ratios.items():
        weight = float(weight)
        if not np.isfinite(weight) or weight <= 0:
            raise ValueError("Intermediate surrogate weights must be finite and positive")
        Tb = float(components[ID].Tb)
        if not 300 + 273.15 < Tb <= 343 + 273.15:
            raise ValueError(f"{ID} Tb={Tb-273.15:.2f} C is outside the intermediate bin")
        ratios[ID] = weight
    return ratios


def validate_cut_fractions(values):
    """Require explicit fractions of represented dry organic biocrude.

    Order: <150, 150-300, 300-343, >=343 C (normal boiling points).
    Water/ash/removal losses must be treated separately. No invented default assay.
    Positive fractions are required because all three shortcut columns are present.
    """
    if values is None:
        raise ValueError("Supply four assay fractions or explicitly labeled sensitivity fractions: <150,150-300,300-343,>=343 C")
    values = tuple(float(x) for x in values)
    if len(values) != 4 or any(not np.isfinite(x) or x <= 0 for x in values):
        raise ValueError("Four finite positive fractions are required")
    if not np.isclose(sum(values), 1.0, rtol=0, atol=1e-8):
        raise ValueError("cutoff_fracs must sum to one; do not normalize missing assay mass")
    return values

def audit_cut_stream(stream, label, show_zero=False):
    """Audit normal-Tb bins; this is composition bookkeeping, not VLE."""
    boundaries = (150.0, 300.0, 343.0)
    names = ('<=150', '150-300', '300-343', '>343')
    cuts = np.zeros(4)
    unknown = 0.0
    rows = []
    for cmp in stream.chemicals:
        mass = float(stream.imass[cmp.ID])
        if not show_zero and mass == 0:
            continue
        try:
            tb = float(cmp.Tb) - 273.15
        except (TypeError, ValueError):
            tb = float('nan')
        if np.isfinite(tb):
            idx = int(np.searchsorted(boundaries, tb, side='left'))
            cuts[idx] += mass
            name = names[idx]
        else:
            unknown += mass
            name = 'UNKNOWN'
        rows.append((cmp.ID, tb, name, mass))
    rows.sort(key=lambda r: (r[1] if np.isfinite(r[1]) else float('inf'), r[0]))
    print(f"\n[AUDIT] {label} | stream={stream.ID}")
    print(f"Total={stream.F_mass:.8f} kg/h | T={stream.T:.3f} K | P={stream.P:.3f} Pa")
    print(f"{'Component':<22}{'Tb C':>12}{'Cut C':>14}{'kg/h':>18}")
    for ID, tb, name, mass in rows:
        tb_text = f'{tb:.3f}' if np.isfinite(tb) else 'UNKNOWN'
        print(f'{ID:<22}{tb_text:>12}{name:>14}{mass:>18.8f}')
    print('Cut masses kg/h:', np.round(cuts, 8))
    print('Cut fractions:', np.round(cuts / stream.F_mass, 8) if stream.F_mass > 0 else np.zeros(4))
    print('Unknown-Tb mass kg/h:', unknown)
    return cuts


def audit_cut_column(col, label):
    feed = col.ins[0]
    top, bottom = col.outs
    print(f"\n[COLUMN AUDIT] {col.ID} | {label}")
    print(f'LHK={col.LHK}, Lr={col.Lr}, Hr={col.Hr}, P={col.P} Pa')
    f = audit_cut_stream(feed, f'{col.ID} inlet')
    d = audit_cut_stream(top, f'{col.ID} overhead')
    b = audit_cut_stream(bottom, f'{col.ID} bottoms')
    print('Mass closure error kg/h:', top.F_mass + bottom.F_mass - feed.F_mass)
    print('Cut closure errors kg/h:', np.round(d + b - f, 8))
    if feed.F_mass > 0:
        print('Actual overhead/feed mass ratio:', top.F_mass / feed.F_mass)
    for ID in col.LHK:
        m = float(feed.imass[ID])
        mt, mb = float(top.imass[ID]), float(bottom.imass[ID])
        print(f'Key {ID}: feed={m:.8f}, top={mt:.8f}, bottom={mb:.8f} kg/h')
        if m > 0:
            print(f'  Top recovery={mt/m:.8f}; bottom recovery={mb/m:.8f}')


def audit_joint_system(system, failed=False):
    """Failed/unreached unit outlets may contain trial, stale or empty states."""
    units = system.flowsheet.unit
    label = 'FAILED attempt: outputs may be trial/stale/unreached' if failed else 'completed simulation'
    audit_cut_stream(units.CrudeSplitter.outs[0], 'CrudeSplitter outlet')
    # Include the entire thermo inventory once, with zero-flow components.
    audit_cut_stream(units.CrudeSplitter.outs[0], 'Full thermo component inventory', show_zero=True)
    for ID in ('CrudeLightDis', 'CrudeHeavyDis', 'BinderCutDis'):
        audit_cut_column(getattr(units, ID), label)


def export_joint_audit_excel(system, fractions, status, error_text='',
                             run_started=None, output_dir=None):
    """Snapshot all streams reachable downstream of HTL biocrude.

    Also includes every input of reached units (H2, catalysts, utilities/feed).
    Failed snapshots may contain trial/stale/empty states; no convergence claim.
    Component mass fractions use each stream's total mass, including water.
    """
    from datetime import datetime
    from pathlib import Path
    import pandas as pd

    now = datetime.now().astimezone()
    if output_dir is None:
        output_dir = Path(__file__).resolve().parent / 'joint_cut_audits'
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    tag = '_'.join(f'{v:.2f}' for v in fractions)
    path = output_dir / f'joint_audit_{now:%Y%m%d_%H%M%S_%f}_{tag}_{status}.xlsx'
    units = system.flowsheet.unit
    # Begin at the HTL biocrude outlet, preserving scaler/pump connections.
    root = units.HTL.outs[2]  # gas, aqueous, biocrude, char
    pending = [root]
    seen_streams, seen_units = {}, {}
    while pending:
        stream = pending.pop()
        if id(stream) in seen_streams:
            continue
        seen_streams[id(stream)] = stream
        sink = getattr(stream, 'sink', None)
        if sink is not None and id(sink) not in seen_units:
            seen_units[id(sink)] = sink
            pending.extend(sink.outs)
    # Side inputs are recorded without traversing backwards to unrelated trains.
    for unit in seen_units.values():
        for stream in unit.ins:
            seen_streams.setdefault(id(stream), stream)

    summaries, compositions, connections, balances = [], [], [], []
    cut_names = ('<=150', '150-300', '300-343', '>343')
    state_note = ('Completed simulation snapshot; qualification not established'
                  if status == 'PASS' else
                  'Failed simulation snapshot: trial, stale or unreached states possible')
    for number, stream in enumerate(seen_streams.values(), 1):
        source = getattr(stream, 'source', None)
        sink = getattr(stream, 'sink', None)
        source_id = getattr(source, 'ID', '')
        sink_id = getattr(sink, 'ID', '')
        total = float(stream.F_mass)
        cut_masses = np.zeros(4)
        unknown = 0.0
        for cmp in stream.chemicals:
            mass = float(stream.imass[cmp.ID])
            try:
                tb = float(cmp.Tb) - 273.15
            except (TypeError, ValueError):
                tb = float('nan')
            cut = 'UNKNOWN'
            if np.isfinite(tb):
                idx = int(np.searchsorted((150., 300., 343.), tb, side='left'))
                cut = cut_names[idx]
                cut_masses[idx] += mass
            else:
                unknown += mass
            # Include zero-flow components so the full thermo list is retained.
            compositions.append(dict(
                Stream_number=number, Stream_ID=stream.ID,
                Source_unit=source_id, Sink_unit=sink_id,
                Component=cmp.ID, Mass_kg_h=mass,
                Mass_fraction=mass/total if total > 0 else float('nan'),
                Tb_C=tb, Normal_Tb_cut_C=cut,
                Nonzero=mass != 0.0, Run_status=status))
        row = dict(Stream_number=number, Stream_ID=stream.ID,
                   Source_unit=source_id, Sink_unit=sink_id,
                   Terminal_stream=sink is None,
                   Mass_kg_h=total, Temperature_K=float(stream.T),
                   Pressure_Pa=float(stream.P), Phase=str(stream.phase),
                   Unknown_Tb_mass_kg_h=unknown, Run_status=status)
        for name, mass in zip(cut_names, cut_masses):
            row[f'Cut_{name}_kg_h'] = mass
            row[f'Cut_{name}_mass_fraction'] = mass/total if total > 0 else float('nan')
        summaries.append(row)
    for unit in seen_units.values():
        inlet = sum(float(st.F_mass) for st in unit.ins)
        outlet = sum(float(st.F_mass) for st in unit.outs)
        balances.append(dict(Unit_ID=unit.ID, Unit_class=type(unit).__name__,
                             Inlet_kg_h=inlet, Outlet_kg_h=outlet,
                             Outlet_minus_inlet_kg_h=outlet-inlet))
        for direction, streams in (('inlet', unit.ins), ('outlet', unit.outs)):
            for port, stream in enumerate(streams):
                connections.append(dict(Unit_ID=unit.ID, Direction=direction,
                                        Port=port, Stream_ID=stream.ID))
    trials = []
    for unit in seen_units.values():
        trials.extend(get_distillation_trials(unit))
    metadata = [
        ('Run_started', run_started or ''),
        ('Snapshot_time', now.isoformat()), ('Status', status),
        ('Error', error_text), ('Fractions', str(tuple(fractions))),
        ('Fraction_order', '<=150,150-300,300-343,>343 C'),
        ('Snapshot_interpretation', state_note),
        ('Mass_fraction_basis', 'Component kg/h divided by total stream kg/h; includes water'),
        ('Tb_interpretation', 'Normal component Tb; not operating temperature or product qualification'),
        ('Scope', 'HTL biocrude outlet and reachable downstream streams, plus side inputs'),
        ('Balance_interpretation', 'Scalers intentionally change throughput; failed states may not close'),
    ]
    # ThermoSTEAM exception strings may contain terminal control characters.
    # Strip XML-illegal controls from text only; keep numeric audit values intact.
    import re
    illegal_controls = re.compile(r'[\x00-\x08\x0b\x0c\x0e-\x1f]')
    def clean_text(value):
        return illegal_controls.sub('', value) if isinstance(value, str) else value
    metadata = [(clean_text(k), clean_text(v)) for k, v in metadata]
    for records in (summaries, compositions, connections, balances, trials):
        for record in records:
            for key, value in record.items():
                record[key] = clean_text(value)
    with pd.ExcelWriter(path, engine='openpyxl') as writer:
        pd.DataFrame(metadata, columns=['Field', 'Value']).to_excel(writer, sheet_name='Run_info', index=False)
        pd.DataFrame(summaries).to_excel(writer, sheet_name='Streams', index=False)
        pd.DataFrame(compositions).to_excel(writer, sheet_name='Components', index=False)
        pd.DataFrame(connections).to_excel(writer, sheet_name='Connections', index=False)
        pd.DataFrame(balances).to_excel(writer, sheet_name='Unit_balances', index=False)
        pd.DataFrame(trials if trials else [{'Status': 'NO_SEARCH_RECORDED'}]).to_excel(writer, sheet_name='Distillation_trials', index=False)
        for ws in writer.book.worksheets:
            ws.freeze_panes = 'A2'
            ws.auto_filter.ref = ws.dimensions
            for col in ws.columns:
                letter = col[0].column_letter
                ws.column_dimensions[letter].width = min(55, max(14, max(len(str(cell.value or '')) for cell in col[:100]) + 2))
    print(f'Excel audit saved: {path}')
    return path


def create_system(
    flowsheet=None,
    feedstock_id="fog",

    # biobinder-style knobs
    central_dry_flowrate=None,
    pilot_dry_flowrate=None,
    decentralized_HTL=False,
    decentralized_upgrading=False,

    # SAF toggles
    include_PSA=True,
    include_EC=False,

    # electricity overrides
    electricity_price=None,
    electricity_GHG=None,

    # optional override of HTL yields (e.g., RF predicted)
    HTL_yields_override=None,
    use_RF=True,
    rf_model_path=None,

    # runtime distillation tuning
    prefer="min_purchase_cost",
    tol=0.05,
    Lr_grid=None,
    Hr_grid=None,

    # Four dry-organic assay/sensitivity fractions; no inferred 300-343 yield.
    cutoff_fracs=None,
    intermediate_surrogate_ratios=None,
    intermediate_route="hydroprocess",  # hydroprocess | chp | sale
    intermediate_price=0.5,             # $/kg, not a fuel-specification claim
    biobinder_price=0.10,               # inherited binder scenario, $/kg
    binder_storage_T=298.15,            # inherited cooling; heated tank not modeled
    tea_overrides=None,
    # Retain legacy kwargs for compatibility.
    **kwargs,
):
    # ---------------------------------------------------------------------
    # 0) Clear registries like biobinder does (prevents stale items)
    # ---------------------------------------------------------------------
    fractions = validate_cut_fractions(cutoff_fracs)
    if intermediate_route not in ("hydroprocess", "chp", "sale", "marine"):
        raise ValueError("intermediate_route must be hydroprocess, chp, sale, or marine")
    if not np.isfinite(intermediate_price) or intermediate_price < 0:
        raise ValueError("intermediate_price must be a finite nonnegative $/kg")
    if not np.isfinite(biobinder_price) or biobinder_price < 0:
        raise ValueError("biobinder_price must be a finite nonnegative $/kg")
    price_dct = dict(_base_price_dct)
    qs.main_flowsheet.clear()
    clear_lca_registries()

    if central_dry_flowrate is None:
        central_dry_flowrate = default_central
    if pilot_dry_flowrate is None:
        pilot_dry_flowrate = default_pilot

    # ---------------------------------------------------------------------
    # 1) Match biobinder CHCU/DHCU/DHDU logic exactly
    # ---------------------------------------------------------------------
    if decentralized_HTL is False:
        if decentralized_upgrading is False:
            flowsheet_ID = "joint_CHCU"
            N_HTL = N_upgrading = 1
        else:
            raise ValueError("Centralized HTL + decentralized upgrading is not a valid configuration.")
    else:
        if decentralized_upgrading is False:
            flowsheet_ID = "joint_DHCU"
            N_HTL = max(1, round(central_dry_flowrate / pilot_dry_flowrate))
            N_upgrading = 1
            pilot_dry_flowrate = central_dry_flowrate / N_HTL
        else:
            flowsheet_ID = "joint_DHDU"
            N_HTL = N_upgrading = 1
            central_dry_flowrate = pilot_dry_flowrate

    # ---------------------------------------------------------------------
    # 2) Flowsheet init
    # ---------------------------------------------------------------------
    if flowsheet is None:
        flowsheet = qs.Flowsheet(flowsheet_ID)
        qs.main_flowsheet.set_flowsheet(flowsheet)
    else:
        qs.main_flowsheet.set_flowsheet(flowsheet)

    _load_process_settings()
    _load_components()

    # -----------------------------------------
    # 3) Feedstock scaling + pricing
    # -----------------------------------------
    scaled_feedstock = qs.WasteStream("scaled_feedstock", price=0.0)
    try:
        scaled_feedstock.price = float(feedprice_dct[feedstock_id])
    except KeyError:
        raise KeyError(f"No feedstock price defined for '{feedstock_id}' in feedprice_dct")

    FeedstockScaler = u_bb.Scaler(
        "FeedstockScaler", ins=scaled_feedstock, outs="feedstock",
        scaling_factor=N_HTL, reverse=True,
    )

    FeedstockTrans = u_bb.Transportation(
        "FeedstockTrans",
        ins=(FeedstockScaler - 0, "feedstock_trans_surrogate"),
        outs=("transported_feedstock",),
        N_unit=N_HTL,
        copy_ins_from_outs=True,
        transportation_unit_cost=0,   # set later if you want
        transportation_distance=1,
    )

    # Process water (biobinder style)
    scaled_process_water = qs.WasteStream("scaled_process_water")
    ProcessWaterScaler = u_bb.Scaler(
        "ProcessWaterScaler", ins=scaled_process_water, outs="htl_process_water",
        scaling_factor=N_HTL, reverse=True,
    )

    feedstock_composition = get_feedstock_composition(feedstock_id)
    print(f"Feedstock used is {feedstock_id}")
    print("Corresponding composition:")
    for k, v in feedstock_composition.items():
        print(f"  {k}: {v:.4f}")

    FeedstockCond = u_bb.Conditioning(
        "FeedstockCond",
        ins=(FeedstockTrans - 0, ProcessWaterScaler.outs[0]),
        outs="conditioned_feedstock",
        feedstock_composition=feedstock_composition,
        feedstock_dry_flowrate=central_dry_flowrate if decentralized_HTL is False else pilot_dry_flowrate,
        target_HTL_solid_loading=0.2,
    )
    FeedstockCond.N_unit = N_HTL

    # ---------------------------------------------------------------------
    # 4) HTL (biobinder CentralizedHTL / PilotHTL)
    # ---------------------------------------------------------------------
    aqueous_composition = {"N": 0.48 / 100}
    aqueous_composition["HTLaqueous"] = 1 - sum(aqueous_composition.values())

    if HTL_yields_override is not None:
        dw_yields = dict(HTL_yields_override)
        yield_source = "explicit override"
    elif use_RF:
        dw_yields = _predict_rf_yields(feedstock_composition, rf_model_path)
        yield_source = "RF predictor"
    else:
        dw_yields = dict(HTL_yields)
        yield_source = "inherited defaults (use_RF=False)"

    HTL_kwargs = dict(
        ID="HTL",
        ins=FeedstockCond.outs[0],
        outs=("gas", "HTL_aqueous", "biocrude", "hydrochar"),
        T=280 + 273.15,
        P=12.4e6,
        tau=15 / 60,
        dw_yields=dw_yields,
        gas_composition={"CO2": 1},
        aqueous_composition=aqueous_composition,
        biocrude_composition={"HTLbiocrude": 1},
        char_composition={"HTLchar": 1},
        internal_heat_exchanging=True,
        eff_T=60 + 273.15,
        eff_P=30 * _psi_to_Pa,
    )

    if decentralized_HTL is False:
        HTL_unit = u_bb.CentralizedHTL
    else:
        HTL_unit = u_bb.PilotHTL
        HTL_kwargs["N_unit"] = N_HTL

    HTL = HTL_unit(**HTL_kwargs)

    # scale gas/char to central totals (biobinder pattern)
    HTLgasScaler = u_bb.Scaler("HTLgasScaler", ins=HTL - 0, outs="scaled_HTLgas", scaling_factor=N_HTL, reverse=False)
    HTLcharScaler = u_bb.Scaler("HTLcharScaler", ins=HTL - 1 + 2, outs="scaled_HTLchar", scaling_factor=N_HTL, reverse=False) \
        if False else u_bb.Scaler("HTLcharScaler", ins=HTL.outs[-1], outs="scaled_HTLchar", scaling_factor=N_HTL, reverse=False)

    # internal HX guardrails (same helpers you use in biobinder)
    limit_internal_hx(HTL, coldest_id="chilled_water", cap_margin_K=2.0, cooler_margin_K=3.0)
    set_eff_hx_temperature(HTL, clamp_to="cooling_water")

    # ---------------------------------------------------------------------
    # 5) Biocrude transport + scaling to central (biobinder pattern)
    #     - This ensures downstream upgrading sees ONE stream representing the total system
    # ---------------------------------------------------------------------
    BiocrudeTrans = u_bb.Transportation(
        "BiocrudeTrans",
        ins=(HTL - 2, "biocrude_trans_surrogate"),
        outs=("transported_biocrude",),
        N_unit=N_HTL,
        transportation_unit_cost=0,
        transportation_distance=1,
    )

    BiocrudeScaler = u_bb.Scaler(
        "BiocrudeScaler",
        ins=BiocrudeTrans - 0,
        outs="scaled_biocrude",
        scaling_factor=N_HTL,
        reverse=False,
    )

    # =========================================================================
    # ================== DOWNSTREAM SAF SECTION STARTS HERE ====================
    # =========================================================================

    # --- crude split ---
    cutoff_fracs = fractions

    CrudeSplitter = JointBiocrudeSplitter(
        "CrudeSplitter",
        ins=BiocrudeScaler - 0,
        outs="splitted_crude",
        biocrude_IDs=("HTLbiocrude",),
        cutoff_fracs=cutoff_fracs,
        biocrude_ratios=_joint_biocrude_ratios(intermediate_surrogate_ratios),
        cutoff_Tbs=(150 + 273.15, 300 + 273.15, 343 + 273.15),
    )

    if len(CrudeSplitter.keys) != 3:
        raise RuntimeError("Local BiocrudeSplitter must support four fractions and three cutoff keys; inspect _units.py before running.")

    CrudePump = qsu.Pump("CrudePump", ins=CrudeSplitter - 0, outs="crude_to_dist", init_with="Stream")

    CrudeLightDis = qsu.ShortcutColumn(
        "CrudeLightDis",
        ins=CrudePump - 0,
        outs=("crude_light", "crude_medium_heavy"),
        LHK=CrudeSplitter.keys[0],
        P=50 * _psi_to_Pa,
        Lr=0.87,
        Hr=0.98,
        k=2,
        is_divided=True,
    )

    CrudeLightFlash = qsu.Flash(
        "CrudeLightFlash",
        ins=CrudeLightDis - 0,
        outs=("crude_light_gas", "crude_light_liq"),
        T=298.15,
        P=101325,
    )
    
    HTLaqScaler = u_bb.Scaler("HTLaqScaler", ins=HTL-1,
        outs="scaled_HTLaq", scaling_factor=N_HTL, reverse=False)
    HTLaqMixer = qsu.Mixer('HTLaqMixer', ins=(HTLaqScaler-0, CrudeLightFlash-1), outs='HTL_aq')

    # --- separate medium from char ---
    CrudeHeavyDis = qsu.ShortcutColumn(
        "CrudeHeavyDis",
        ins=CrudeLightDis - 1,
        outs=("crude_medium", "crude_above300"),
        LHK=("4M-PHYNO", "INDOLE"),
        P=50 * _psi_to_Pa,
        Lr=0.8,
        Hr=0.85,
        k=2,
        is_divided=True,
    )

    ratio0 = CrudeSplitter.cutoff_fracs[1] / sum(CrudeSplitter.cutoff_fracs[1:])
    add_runtime_LHK_LrHr_spec(
        splitter=CrudeSplitter,
        col=CrudeHeavyDis,
        idx=1,
        fallback_LHK=CrudeSplitter.keys[1],
        target_ratio=ratio0,
        tol=float(tol),
        Lr_grid=Lr_grid, Hr_grid=Hr_grid,
        prefer=str(prefer),
        require_design=True,
        require_cost=True, max_purchase_cost=1e7, cost_band=0.0, cost_tiebreak="min_err",
        vle_screen=True,
        design_sanity=True,
    )

    # =========================================================================
    # Hydrocracking (SAF downstream)
    # =========================================================================
    # Added 343 C normal-boiling-point boundary. Column pressure stays inherited.
    # These are shortcut separations with finite recoveries, not exact cut bins.
    BinderCutDis = qsu.ShortcutColumn(
        "BinderCutDis", ins=CrudeHeavyDis-1,
        outs=("intermediate_300_343", "hot_biobinder"),
        LHK=CrudeSplitter.keys[2], P=50*_psi_to_Pa,
        Lr=0.8, Hr=0.85, k=2, is_divided=True)
    def audit_binder_feed():
        feed = BinderCutDis.ins[0]
        cuts = audit_cut_stream(feed, 'BinderCutDis inlet BEFORE recovery search')
        pair = CrudeSplitter.keys[2]
        print('Primary LHK:', pair)
        print('Fallback LHK:', ('C18H38', 'C18FACID'))
        print('Specified top/feed target:', fractions[2] / sum(fractions[2:]))
        print('Actual feed Tb <=343 C fraction:', sum(cuts[:3]) / feed.F_mass if feed.F_mass > 0 else float('nan'))
        for ID in pair:
            print(f'Key {ID}: Tb={float(feed.chemicals[ID].Tb)-273.15:.3f} C; mass={float(feed.imass[ID]):.8f} kg/h')

    BinderCutDis.check_LHK = False
    BinderCutDis.add_specification(audit_binder_feed)
    add_runtime_LHK_LrHr_spec(
        splitter=CrudeSplitter, col=BinderCutDis, idx=2,
        fallback_LHK=('C18H38', 'C18FACID'),
        target_ratio=fractions[2]/sum(fractions[2:]),
        tol=0.10, x_min=0.0,
        Lr_grid=[0.995, 0.99, 0.98, 0.95, 0.90, 0.85, 0.80, 0.70],
        Hr_grid=[0.99, 0.98, 0.95, 0.90, 0.85, 0.80, 0.70],
        prefer='min_purchase_cost', require_design=True, require_cost=True, max_purchase_cost=1e7, cost_band=0.0, cost_tiebreak="min_err",
        vle_screen=True, design_sanity=True)
    BiobinderHX = qsu.HXutility("BiobinderHX", ins=BinderCutDis-1,
        outs="cooled_biobinder", T=binder_storage_T)
    biobinder = qs.WasteStream("biobinder", price=biobinder_price)
    BiobinderStorage = qsu.StorageTank("BiobinderStorage", ins=BiobinderHX-0,
        outs=biobinder, tau=72, vessel_material="Stainless steel",
        include_construction=False)
    hc_feeds = [CrudeHeavyDis-0]
    intermediate_chp = []
    if intermediate_route == "hydroprocess":
        hc_feeds.append(BinderCutDis-0)
    elif intermediate_route == "chp":
        intermediate_chp.append(BinderCutDis-0)
    else:
        IntermediateHX = qsu.HXutility("IntermediateHX", ins=BinderCutDis-0,
            outs="cooled_intermediate", T=298.15, rigorous=True)
        intermediate = qs.WasteStream("intermediate", price=intermediate_price)
        IntermediateTank = qsu.StorageTank("IntermediateTank", ins=IntermediateHX-0,
            outs=intermediate, tau=72, vessel_material="Stainless steel",
            include_construction=False)
    if intermediate_route == 'marine':
        from scenario_pricing import marine_price
        def set_marine_price():
            intermediate.price = marine_price(IntermediateTank.ins[0])
        IntermediateTank.add_specification(set_marine_price)
    SAFfeedMixer = qsu.Mixer("SAFfeedMixer", ins=hc_feeds, outs="SAF_precursor")

    HCcatalyst_in = qs.WasteStream("HCcatalyst_in", HCcatalyst=1, price=price_dct["HCcatalyst"])

    HC = u_saf.Hydroprocessing(
        "HC",
        ins=(SAFfeedMixer - 0, "H2_HC", HCcatalyst_in),
        outs=("HC_out", "HCcatalyst_out"),
        T=400 + 273.15,
        P=1500 * _psi_to_Pa,
        WHSV=0.625,
        catalyst_ID="HCcatalyst",
        catalyst_lifetime=5 * uptime_ratio * 365 * 24,
        hydrogen_rxned_to_inf_oil=0.0111,
        hydrogen_ratio=5.556,
        include_PSA=include_PSA,
        gas_yield=0.2665,
        oil_yield=0.7335,
        gas_composition={
            "CO2": 1 - 0.08809,
            "CH4": 0.02280, "C2H6": 0.02923,
            "C3H8": 0.01650, "C4H10": 0.00870,
            "TWOMBUTAN": 0.00408, "NPENTAN": 0.00678,
        },
        oil_composition={
            'TWOMPENTA': 0.00408, 'HEXANE': 0.00408,
            'TWOMHEXAN': 0.00408, 'HEPTANE': 0.00408,
            'CC6METH': 0.01020, 'PIPERDIN': 0.00408,
            'TOLUENE': 0.01020, 'THREEMHEPTA': 0.01020,
            'OCTANE': 0.01020, 'ETHCYC6': 0.00408,
            'ETHYLBEN': 0.02040, 'OXYLENE': 0.01020,
            'C9H20': 0.00408, 'PROCYC6': 0.00408,
            'C3BENZ': 0.01020, 'FOURMONAN': 0,
            'C10H22': 0.00203, 'C4BENZ': 0.01223,
            'C11H24': 0.02040, 'C10H12': 0.02040,
            'C12H26': 0.02040, 'OTTFNA': 0.01020,
            'C6BENZ': 0.02040, 'OTTFSN': 0.02040,
            'C7BENZ': 0.02040, 'C8BENZ': 0.02040,
            'C10H16O4': 0.01837, 'C15H32': 0.06120,
            'C16H34': 0.18360, 'C17H36': 0.08160,
            'C18H38': 0.04080, 'C19H40': 0.04080,
            'C20H42': 0.10200, 'C21H44': 0.04080,
            'TRICOSANE': 0.04080, 'C24H38O4': 0.00817,
            'C26H42O4': 0.01020, 'C30H62': 0.00203,
        },
        aqueous_composition={"Water": 1},
        internal_heat_exchanging=False,
        use_decorated_cost="Hydrocracker",
        tau=15/60,
        V_wf=0.4,
        length_to_diameter=2,
        vessel_material="Stainless steel 316",
        vessel_type="Vertical",
    )

    HC_HX = qsu.HXutility("HC_HX", ins=HC - 0, outs="cooled_HC_eff", T=60 + 273.15, init_with="Stream", rigorous=True)
    HC.eff_hx = _DummyEffHX("HC_eff_hx", bst.settings.thermo)
    HC_IV = IsenthalpicValve("HC_IV", ins=HC_HX - 0, outs="cooled_depressed_HC_eff", P=30 * _psi_to_Pa, vle=True)
    HCflash = qsu.Flash("HCflash", ins=HC_IV - 0, outs=("HC_fuel_gas", "HC_liquid"), T=60.2 + 273.15, P=30 * _psi_to_Pa)
    HCpump = qsu.Pump("HCpump", ins=HCflash - 1, init_with="Stream")

    HCliquidSplitter = qsu.Splitter(
        "HCliquidSplitter",
        ins=HCpump - 0,
        outs=("HC_ww", "HC_oil"),
        split={"H2O": 1},
        init_with="Stream",
    )

    # =========================================================================
    # Hydrotreating
    # =========================================================================
    HTcatalyst_in = qs.WasteStream("HTcatalyst_in", HTcatalyst=1, price=price_dct["HTcatalyst"])
    oil_fracs = [0.3455, 0.4479, 0.2066]

    HT = u_saf.Hydroprocessing(
        "HT",
        ins=(HCliquidSplitter - 1, "H2_HT", HTcatalyst_in),
        outs=("HTout", "HTcatalyst_out"),
        WHSV=0.625,
        catalyst_lifetime=2 * uptime_ratio * 365 * 24,
        catalyst_ID="HTcatalyst",
        T=300 + 273.15,
        P=1500 * _psi_to_Pa,
        hydrogen_rxned_to_inf_oil=0.0207,
        hydrogen_ratio=3,
        include_PSA=include_PSA,
        gas_yield=0.2143,
        oil_yield=0.8637,
        gas_composition={"CO2": 0.03880, "CH4": 0.00630},
        oil_composition={"Gasoline": oil_fracs[0], "Jet": oil_fracs[1], "Diesel": oil_fracs[2]},
        aqueous_composition={"Water": 1},
        internal_heat_exchanging=False,
        use_decorated_cost="Hydrotreater",
        tau=0.5,
        V_wf=0.4,
        length_to_diameter=2,
        vessel_material="Stainless steel 316",
        vessel_type="Vertical",
    )

    HT_HX = qsu.HXutility("HT_HX", ins=HT - 0, outs="cooled_HT_eff", T=60 + 273.15, init_with="Stream", rigorous=True)
    HT.eff_hx = _DummyEffHX("HT_eff_hx", bst.settings.thermo)
    HT_IV = IsenthalpicValve("HT_IV", ins=HT_HX - 0, outs="cooled_depressed_HT_eff", P=717.4 * _psi_to_Pa, vle=True)
    HTflash = qsu.Flash("HTflash", ins=HT_IV - 0, outs=("HT_fuel_gas", "HT_liquid"), T=43 + 273.15, P=55 * _psi_to_Pa)
    HTpump = qsu.Pump("HTpump", ins=HTflash - 1, init_with="Stream")

    HTliquidSplitter = qsu.Splitter(
        "HTliquidSplitter",
        ins=HTpump - 0,
        outs=("HT_ww", "HT_oil"),
        split={"H2O": 1},
        init_with="Stream",
    )

    GasolineDis = qsu.ShortcutColumn(
        "OilLightDis",
        ins=HTliquidSplitter - 1,
        outs=("hot_gasoline", "jet_diesel"),
        LHK=("Gasoline", "Jet"),
        Lr=0.99,
        Hr=0.99,
        k=2,
        is_divided=True,
    )
    GasolineFlash = qsu.Flash("GasolineFlash", ins=GasolineDis - 0, outs=("", "cooled_gasoline"), T=298.15, P=101325)

    JetDis = qsu.ShortcutColumn(
        "JetDis",
        ins=GasolineDis - 1,
        outs=("hot_jet", "hot_diesel"),
        LHK=("Jet", "Diesel"),
        Lr=0.99,
        Hr=0.99,
        k=2,
        is_divided=True,
    )
    JetFlash = qsu.Flash("JetFlash", ins=JetDis - 0, outs=("", "cooled_jet"), T=298.15, P=101325)
    DieselHX = qsu.HXutility("DieselHX", ins=JetDis - 1, outs="cooled_diesel", T=298.15, init_with="Stream", rigorous=True)

    # =========================================================================
    # Products
    # =========================================================================
    def _set_price_from_gal(stream, price_per_gal):
        """Convert $/gal to $/kg using the stream density."""
        rho = getattr(stream, "rho", None)
        if rho is None or (not np.isfinite(rho)) or rho <= 0:
            fallback_rho = {
            "gasoline": 740.0,   # kg/m3
            "jet": 800.0,
            "diesel": 832.0,
             }
            rho = fallback_rho.get(stream.ID, 800.0)
        rho_kg_per_gal = rho / 264.172
        stream.price = float(price_per_gal) / float(rho_kg_per_gal)

        if rho_kg_per_gal <= 0:
           raise ValueError(f"Non-positive density for stream {stream.ID}")
        stream.price = float(price_per_gal) / float(rho_kg_per_gal)
    GasolinePC = qsu.PhaseChanger("GasolinePC", ins=GasolineFlash - 1)
    gasoline = qs.WasteStream("gasoline", Gasoline=1)
    GasolineTank = qsu.StorageTank("GasolineTank", ins=GasolinePC - 0, outs=(gasoline,),
                                   tau=3*24, init_with="WasteStream", vessel_material="Carbon steel")

    JetPC = qsu.PhaseChanger("JetPC", ins=JetFlash - 1)
    jet = qs.WasteStream("jet", Jet=1)
    JetTank = qsu.StorageTank("JetTank", ins=JetPC - 0, outs=(jet,),
                              tau=3*24, init_with="WasteStream", vessel_material="Carbon steel")

    DieselPC = qsu.PhaseChanger("DieselPC", ins=DieselHX - 0)
    diesel = qs.WasteStream("diesel", Diesel=1)
    DieselTank = qsu.StorageTank("DieselTank", ins=DieselPC - 0, outs=(diesel,),
                                 tau=3*24, init_with="WasteStream", vessel_material="Carbon steel")

    mixed_fuel = qs.WasteStream("mixed_fuel", price=0.0)
    GasolineTank._run()
    JetTank._run()
    DieselTank._run()
    def _update_product_prices_and_mixed_fuel():        
      
        # MARKET / CUTS MODE: gasoline, jet, diesel are sold products
        if SAF_ECON_MODE == "market":
           _set_price_from_gal(gasoline, price_dct["gasoline"])
           _set_price_from_gal(jet,      price_dct["jet"])
           _set_price_from_gal(diesel,   price_dct["diesel"])
           mixed_fuel.price = 0.0

        # MFSP MODE: cuts are bookkeeping only, mixed fuel is the economic product
        elif SAF_ECON_MODE == "mfsp":
           gasoline.price = 0.0
           jet.price = 0.0
           diesel.price = 0.0
           mixed_fuel.price = 0.0

        else:
           raise ValueError(f"Unknown SAF_ECON_MODE: {SAF_ECON_MODE}")    

        # Always keep a bookkeeping mixed-fuel stream for reporting
        mixed_fuel.mix_from((gasoline, jet, diesel))

    # attach to last product unit so it updates each simulation
    DieselTank.add_specification(_update_product_prices_and_mixed_fuel)
    DieselTank.run_after_specifications = True
     

    # =========================================================================
    # Electrochemical Unit
    # =========================================================================
    # All wastewater streams
    ww_streams = [HTLaqMixer-0, HCliquidSplitter-0, HTliquidSplitter-0]
    # Wastewater sent to municipal wastewater treatment plant
    ww_to_disposal = qs.WasteStream('ww_to_disposal')

    WWmixer = qsu.Mixer('WWmixer', ins=ww_streams)
    
    fuel_gases = [
        HTLgasScaler-0, CrudeLightFlash-0, # aggregate HTL gases
        HCflash-0, HTflash-0, # post-hydroprocessing gases
        GasolineFlash-0, JetFlash-0, # final distillation fuel gases
        ]
       
    recovered_N = qs.WasteStream('recovered_N', price=price_dct['N'])
    recovered_P = qs.WasteStream('recovered_P', price=price_dct['P'])
    recovered_K = qs.WasteStream('recovered_K', price=price_dct['K'])

    EC = u_saf.SAFElectrochemical(
        'EC',
        ins=(WWmixer-0, 'EC_replacement_surrogate'),
        outs=('EC_gas', 'EC_H2', recovered_N, recovered_P, recovered_K, ww_to_disposal),
        COD_removal=0.95, # assumed
        N_recovery=0.8,
        P_recovery=0.99,
        K_recovery=0.8,
        include_PSA=include_PSA,
        PSA_efficiency=0.95,
        )
    EC.register_alias('Electrochemical')
    EC.include_PSA_cost = False # HC/HT has PSA
    fuel_gases.append(EC-0)
    if include_EC is False:
        EC.skip = True
    else:
        EC.skip = False
        if type(include_EC) is dict:
            for attr, val in include_EC.items():
                setattr(EC, attr, val)

    def adjust_prices():
        FeedstockTrans._run()        
        # Centralized HTL and upgrading, transport feedstock
        if decentralized_HTL is False:
            dw_price = price_dct['trans_feedstock'] # $/dry mass
            factor = 1 - FeedstockTrans.ins[0].imass['Water']/FeedstockTrans.ins[0].F_mass
            FeedstockTrans.transportation_unit_cost = dw_price * factor
            BiocrudeTrans.transportation_unit_cost = 0
        # Decentralized HTL, centralized upgrading, transport biocrude
        elif decentralized_upgrading is False:
            FeedstockTrans.transportation_unit_cost = 0
            GGE_price = price_dct['trans_biocrude'] # $/GGE
            # 1e3 to convert from kJ/hr to MJ/hr, 264.172 is m3/hr to gal/hr
            factor = BiocrudeTrans.ins[0].HHV/1e3/(BiocrudeTrans.ins[0].F_vol*264.172)/_HHV_per_GGE
            BiocrudeTrans.transportation_unit_cost = GGE_price * factor
        # Decentralized HTL and upgrading, no transportation needed
        else:
            FeedstockTrans.transportation_unit_cost = BiocrudeTrans.transportation_unit_cost = 0
        # Wastewater
        ww_to_disposal.source._run()
        COD_mass_content = ww_to_disposal.COD*ww_to_disposal.F_vol/1e3 # mg/L*m3/hr to kg/hr
        factor = COD_mass_content/ww_to_disposal.F_mass
        ww_to_disposal.price = price_dct['COD']*factor
        ww_to_disposal_item = qs.ImpactItem.get_item('ww_to_disposal_item')
        try: ww_to_disposal_item.CFs['GWP'] = gwp_dct['COD']*factor
        except: pass
    ww_to_disposal.source.add_specification(adjust_prices)

    GasMixer = qsu.Mixer('GasMixer', ins=fuel_gases, outs=('waste_gases'))

    # =========================================================================
    # Facilities
    # =========================================================================
    
    # Adding HXN only saves cents/GGE with HTL internal HX, eliminate for simpler system
    # HXN = qsu.HeatExchangerNetwork('HXN', T_min_app=86, force_ideal_thermo=True)
    # 86 K: Jones et al. PNNL, 2014
    
    natural_gas = qs.WasteStream('natural_gas', CH4=1, price=price_dct['natural_gas'])
    solids_to_disposal = qs.WasteStream('solids_to_disposal', price=price_dct['solids'])
    CHPmixer = qsu.Mixer('CHPmixer', ins=[GasMixer-0, HTLcharScaler-0] + intermediate_chp)
    CHP = qsu.CombinedHeatPower('CHP', 
                                ins=(CHPmixer-0, natural_gas, 'air'),
                                outs=('gas_emissions', solids_to_disposal),
                                init_with='WasteStream',
                                supplement_power_utility=False)
    
    H2C = u_saf.HydrogenCenter(
        'H2C',
        process_H2_streams=(HC.ins[1], HT.ins[1]),
        recycled_H2_streams=EC-1,
        )
    H2C.register_alias('HydrogenCenter')
    H2C.makeup_H2_price = H2C.excess_H2_price = price_dct['H2'] # expected H2 price
    # H2C.makeup_H2_price = H2C.excess_H2_price = 33.4 # current H2 price
    PWC = u_bb.ProcessWaterCenter(
        'ProcessWaterCenter',
        process_water_streams=scaled_process_water,
        process_water_price=price_dct['process_water']
        )
    PWC.register_alias('PWC')
    @PWC.add_specification
    def run_scalers():
        FeedstockScaler._run()
        ProcessWaterScaler._run()
        PWC._run()
    
    
    # =========================================================================
    # System, TEA, LCA
    # =========================================================================   
    sys = qs.System.from_units(
        'sys',
        units=list(flowsheet.unit),
        operating_hours=365*24*uptime_ratio,
        )
    for unit in sys.units: unit.include_construction = False
    joint_tea_kwargs = dict(tea_kwargs)
    joint_tea_kwargs.update(tea_overrides or {})
    tea = create_tea(sys, cls=BiobinderTEA, **joint_tea_kwargs)
    land_factor = 115000/24 #115000 gal/day, PNNL 32731
    tea.land = lambda: 90000+(tea.system.flowsheet.unit.BiocrudeScaler.outs[0].F_vol*264.172)/land_factor*4.692*1e6 # 6% of 78.2M$, PNNL 32731
    

    # Add characterization factors for each impact item
    clear_lca_registries()
    GWP = qs.ImpactIndicator('GWP',
                             alias='GlobalWarmingPotential',
                             method='GREET',
                             category='environmental impact',
                             unit='kg CO2-eq',)

    feedstock_item = qs.StreamImpactItem(
        ID='feedstock_item',
        linked_stream=scaled_feedstock,
        # feedstock, landfill, composting, anaerobic_digestion
        # may or may not be good to assume landfill offsetting
        GWP=-gwp_dct['landfill'],
        )
    trans_feedstock_item = qs.StreamImpactItem(
        ID='trans_feedstock_item',
        linked_stream=FeedstockTrans.ins[1],
        GWP=gwp_dct['trans_feedstock'],
        )
    makeup_H2_item = qs.StreamImpactItem(
        ID='makeup_H2_item',
        linked_stream=H2C.ins[0],
        GWP=gwp_dct['H2'],
        )
    excess_H2_item = qs.StreamImpactItem(
        ID='excess_H2_item',
        linked_stream=H2C.outs[1],
        GWP=-gwp_dct['H2'],
        )
    trans_biocrude_item=qs.StreamImpactItem(
        ID='biocrude_trans_surrogate_item',
        linked_stream=BiocrudeTrans.ins[1],
        GWP=gwp_dict['trans_biocrude'],
        )
    HCcatalyst_item = qs.StreamImpactItem(
        ID='HCcatalyst_item',
        linked_stream=HC.ins[-1],
        GWP=gwp_dct['HCcatalyst'],
        )
    HTcatalyst_item = qs.StreamImpactItem(
        ID='HTcatalyst_item',
        linked_stream=HT.ins[-1],
        GWP=gwp_dct['HTcatalyst'],
        )
    natural_gas_item = qs.StreamImpactItem(
        ID='natural_gas_item',
        linked_stream=natural_gas,
        GWP=gwp_dct['natural_gas'],
        )
    # Assume no impacts from process water
    # process_water_item = qs.StreamImpactItem(
    #     ID='process_water_item',
    #     linked_stream=PWC.ins[-1],
    #     GWP=gwp_dct['process_water'],
    #     )
    ww_to_disposal_item = qs.StreamImpactItem(
        ID='ww_to_disposal_item',
        linked_stream=ww_to_disposal,
        GWP=gwp_dct['COD'], # will be updated based on COD content
        )
    solids_to_disposal_item = qs.StreamImpactItem(
        ID='solids_to_disposal_item',
        linked_stream=CHP.outs[1],
        GWP=gwp_dct['solids'],
        )
    qs.PowerUtility.price = electricity_price if electricity_price is not None else qs.PowerUtility.price
    e_item = qs.ImpactItem(
        ID='e_item',
        GWP=electricity_GHG if electricity_GHG is not None else gwp_dct['electricity'],
        )
    steam_item = qs.ImpactItem(
        ID='steam_item',
        GWP=gwp_dct['steam'],
        )
    cooling_item = qs.ImpactItem(
        ID='cooling_item',
        GWP=gwp_dct['cooling'],
        )
    recovered_N_item = qs.StreamImpactItem(
        ID='recovered_N_item',
        linked_stream=recovered_N,
        GWP=gwp_dct['N'],
        )
    recovered_P_item = qs.StreamImpactItem(
        ID='recovered_P_item',
        linked_stream=recovered_P,
        GWP=gwp_dct['P'],
        )
    recovered_K_item = qs.StreamImpactItem(
        ID='recovered_K_item',
        linked_stream=recovered_K,
        GWP=gwp_dct['K'],
        )
    # LCA adjustment based on system configuration
    fake_stream= qs.SanStream('Nothing', price=0)
    nothing_item = qs.StreamImpactItem(
        ID='nothing_item',
        linked_stream=fake_stream,
        GWP=0,
        )
    if decentralized_HTL is False:
    # Centralized HTL, Centralized upgrading
        trans_feedstock_item.linked_stream = FeedstockTrans.ins[1]  # feedstock transportation stream
        trans_biocrude_item.linked_stream = fake_stream             # No biocrude transportation
    elif decentralized_upgrading is False:
    # Decentralized HTL, centralized upgrading
        trans_feedstock_item.linked_stream = fake_stream            # No feedstock transportation
        trans_biocrude_item.linked_stream = BiocrudeTrans.ins[1]    #biocrude tranportation stream
    else:
    # Fully decentralized (no transportation needed)
        trans_feedstock_item.linked_stream = fake_stream            # No feedstock transportation
        trans_biocrude_item.linked_stream = fake_stream             # No biocrude transportation

    lifetime = tea.duration[1]-tea.duration[0]
    lca = qs.LCA(
        system=sys,
        lifetime=lifetime,
        uptime_ratio=uptime_ratio,
        simulate_system=False,
        e_item=lambda:(sys.get_electricity_consumption()-sys.get_electricity_production())*lifetime,
        steam_item=lambda:sys.get_heating_duty()/1000*lifetime, # kJ/yr to MJ/yr, include natural gas, but all offset in CHP
        cooling_item=lambda:sys.get_cooling_duty()/1000*lifetime, # kJ/yr to MJ/yr
        )
    
    _INTEGRATION_ASSUMPTIONS[sys] = dict(cutoff_fracs=fractions,
        HTL_yield_source=yield_source, HTL_yields=dict(dw_yields),
        distillation_setup="Shared adaptive helper + qsu.ShortcutColumn",
        prefer=prefer, tolerance=tol, Lr_grid=Lr_grid, Hr_grid=Hr_grid,
        intermediate_surrogate_ratios=dict(intermediate_surrogate_ratios),
        intermediate_route=intermediate_route, intermediate_price=intermediate_price,
        biobinder_price=biobinder_price, binder_storage_T=binder_storage_T,
        hydroprocessing_yields="Inherited SAF fixed yields; cut-specific validation needed",
        chp_boundary="Aggregate central gas/char/aqueous collection; logistics not added")
    return sys




def simulate_and_report(sys):
    """Return joint-plant IRR at fixed sale prices; never solve MSP first.

    TEA.solve_IRR() is the model's unlevered project IRR. Financing assumptions
    affect NPV/cash flow differently; do not label this equity IRR.
    Exceptions (e.g. no IRR root) propagate rather than reporting a fabricated IRR.
    """
    sys.simulate()
    products = {}
    for name in ("gasoline", "jet", "diesel", "biobinder", "intermediate"):
        stream = getattr(sys.flowsheet.stream, name, None)
        if stream is not None:
            products[name] = dict(kg_per_hr=stream.F_mass,
                metric_tonnes_per_year=stream.F_mass*sys.operating_hours/1000,
                price_USD_per_kg=stream.price,
                sales_USD_per_year=stream.cost*sys.operating_hours)
    tea = sys.TEA
    result = dict(project_IRR=float(tea.solve_IRR()), NPV=float(tea.NPV),
        annual_sales=float(tea.sales), annual_operating_cost=float(tea.AOC),
        products=products, assumptions=_INTEGRATION_ASSUMPTIONS.get(sys, {}))
    print(result)
    return result

INTEGRATION_NOTES = """
Distillation uses qsu.ShortcutColumn and the shared adaptive helper in
exposan.biobinder_ml._distill_utils_flex. Candidate pairs are generated from
actual column inlet composition at runtime. Minimum purchase-cost selection
retains the configured pressure, recovery grids, tolerance and feasibility
checks. Distillation trials are available through get_distillation_trials.

Four fractions: <150, 150-300, 300-343, >343 C, using <= cutoff conventions from
original splitter. Example (0.01,0.59,0.10,0.30) is an assumed assay scenario.
Existing C18H38 represents the intermediate ONLY as a surrogate sensitivity.
Its normal Tb is in 300-343 C; no evidence of pure C18H38 HTL cut composition is
claimed. Local JointBiocrudeSplitter fixes empty-bin mass loss, duplicate-fraction
collisions, and verifies splitter mass closure without changing original _units.py.

150-300 stream -> inherited HC -> HT -> gasoline/jet/diesel separation/storage.
300-343 -> hydroprocess (default), chp, or sale (default intermediate price=0).
>343 -> binder cooling/storage at inherited 0.10 USD/kg scenario price.
One shared HTL/TEA/utilities/wastewater section. RF output dictionaries are local.
Aggregate DHCU gas/char/aqueous routing and central land sizing remain corrected;
central gas/char/aqueous collection logistics are not separately costed.
Bookkeeping mixed_fuel price remains zero to avoid duplicate fuel revenue.
IRR reporter keeps sale prices rather than solving MSP first.

Limitations: fixed inherited hydroprocessing yields, H2 demand/catalyst assumptions,
HC->HT order, zero-duty internal HC/HT heat exchangers, ambient binder storage and
SAF-inherited LCA. No new policy/product qualification credits. Changes to cut
chemistry change separation/energy properties; modeled jet is a candidate stream.
PNNL pathway references support upgrading concepts, not this exact joint layout:
PNNL-32731, https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-32731.pdf
DOE Thorson 2023, https://www.energy.gov/sites/default/files/2023-05/beto-13-project-peer-review-waste-apr-2023-thorson.pdf
Marine blend experiments: doi:10.1016/j.fuel.2025.137221.
Binder blend experiments: doi:10.1021/acssuschemeng.5c03860.
NIST C18 surrogate physical properties: https://webbook.nist.gov/cgi/cbook.cgi?ID=C593453&Mask=4&Units=SI
Syntax and isolated wiring checks completed; full simulation not run in preparation
environment lacking EXPOsan/QSDsan/BioSTEAM/custom modules. No IRR result claimed.
"""

if __name__ == '__main__':
    # Previous single-case run -- commented out.
    # sys = create_system(
    #     feedstock_id="food", decentralized_HTL=True,
    #     decentralized_upgrading=False, include_PSA=True, include_EC=False,
    #     cutoff_fracs=(0.01, 0.50, 0.10, 0.39),
    #     intermediate_surrogate_ratios={"C18H38": 1.0},
    #     intermediate_route="hydroprocess", use_RF=True)
    # result = simulate_and_report(sys)

    # Light, 150-300 C, 300-343 C, >343 C; explicit sensitivity assumptions.
    fraction_cases = [(0.01, 0.59, 0.10, 0.30)]
    from datetime import datetime
    quick_results = []
    for fractions in fraction_cases:
        print(f"\n{'=' * 75}\nTesting cutoff_fracs={fractions}")
        case_system = None
        failed = True
        error_text = ''
        run_started = datetime.now().astimezone().isoformat()
        try:
            case_system = create_system(
                feedstock_id="food", decentralized_HTL=True,
                decentralized_upgrading=False,
                include_PSA=True, include_EC=False,
                cutoff_fracs=fractions,
                intermediate_surrogate_ratios={"C18H38": 1.0},
                intermediate_route="chp", use_RF=True)
            result = simulate_and_report(case_system)
            print(f"\nProject IRR: {result['project_IRR']:.2%}")
            print(f"NPV: ${result['NPV']:,.0f}")
            print(f"Annual sales: ${result['annual_sales']:,.0f}")
            print(f"Annual operating cost: ${result['annual_operating_cost']:,.0f}")
            failed = False
            quick_results.append((fractions, 'PASS', result))
            print(f'PASS: {fractions}')
        except Exception as error:
            message = f'{type(error).__name__}: {error}'
            error_text = message
            quick_results.append((fractions, 'FAIL', message))
            print(f'FAIL: {fractions}\n{message}')
        finally:
            if case_system is not None:
                try:
                    audit_joint_system(case_system, failed=failed)
                except Exception as audit_error:
                    # Preserve the original simulation failure and continue cases.
                    print(f'Audit could not complete: {type(audit_error).__name__}: {audit_error}')
                try:
                    export_joint_audit_excel(
                        case_system, fractions,
                        status='FAIL' if failed else 'PASS',
                        error_text=error_text, run_started=run_started)
                except Exception as export_error:
                    print(f'Excel export failed: {type(export_error).__name__}: {export_error}')

    print(f"\n{'=' * 75}\nQUICK CHECK SUMMARY")
    for fractions, status, _ in quick_results:
        print(f'{status:4s} | {fractions}')
