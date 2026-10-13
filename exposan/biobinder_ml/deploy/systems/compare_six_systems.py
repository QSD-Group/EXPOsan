# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from pathlib import Path
from datetime import datetime
import contextlib
import io
import json
import traceback
import warnings
import numpy as np
import flex_sys_compare as joint
import SAF_FLEX_compare as saf
import Dist_flex_compare as binder
from scenario_pricing import apply_price
from comparison_reporting import collect,failed_snapshot,export

# ONE EDITABLE COMPARISON BASIS. Rates are dry kg/h, prices USD/kg wet.
FEEDSTOCK='food'
DECENTRALIZED_HTL=True
CENTRAL_DRY_KG_HR=None  # None uses joint model default; exact realized rates logged
PILOT_DRY_KG_HR=None
FOUR_CUTS=(0.01,0.59,0.10,0.30) # <=150,150-300,300-343,>343 C
INTERMEDIATE_ROUTE='chp' # hydroprocess | chp | sale
INTERMEDIATE_PRICE=0.50 # explicit illustrative assumption, USD/kg; only sale uses it
BINDER_PRICE=0.10
FEEDSTOCK_PRICE_WET_TONNE=None # None inherits joint feedprice; negative=tipping income
COMMON_TEA_OVERRIDES={} # e.g. {'IRR':0.08,'duration':(2020,2050)}
INCLUDE_SECOND_REFINERY_ROUTE=True # False produces six cases
OUTPUT_DIR=Path(__file__).resolve().parent/'comparison_results'


def configure(system,basis):
    """Set identical shared price and TEA assumptions; leave pathway costs native."""
    fs=system.flowsheet.stream
    feed=fs.scaled_feedstock
    feed.price=basis['feed_price_USD_kg']
    for name in ('biobinder',):
        stream=getattr(fs,name,None)
        if stream is not None:stream.price=BINDER_PRICE
    tea=system.TEA
    for key,value in basis['tea_inputs'].items():setattr(tea,key,value)
    system.operating_hours=basis['hours']
    system.flowsheet.stream.natural_gas.price = 0.14
    # Land uses total biocrude volume for all three systems, not a pilot volume.
    scaler=system.flowsheet.unit.BiocrudeScaler
    tea.land=lambda:90000+scaler.outs[0].F_vol*264.172/(115000/24)*4.692e6
    # Verify actual feed basis after run; do not assume matching factory defaults.


def run_comparison():
    stamp=datetime.now().astimezone().strftime('%Y%m%d_%H%M%S_%f')
    run_dir=OUTPUT_DIR/stamp;run_dir.mkdir(parents=True,exist_ok=True)
    common=dict(feedstock_id=FEEDSTOCK,decentralized_HTL=DECENTRALIZED_HTL,
        decentralized_upgrading=False,central_dry_flowrate=CENTRAL_DRY_KG_HR,
        pilot_dry_flowrate=PILOT_DRY_KG_HR)
    a,b,c,d=joint.validate_cut_fractions(FOUR_CUTS)
    cases=[('Joint_CHP',joint,dict(common,include_PSA=True,include_EC=False,
        cutoff_fracs=(a,b,c,d),intermediate_surrogate_ratios={'C18H38':1.0},
        intermediate_route=INTERMEDIATE_ROUTE,intermediate_price=INTERMEDIATE_PRICE,
        biobinder_price=BINDER_PRICE,use_RF=True)),
        ('SAF',saf,dict(common,include_PSA=True,include_EC=False,cutoff_fracs=(a,b,c+d))),
        ('Binder_diesel',binder,dict(common,skip_EC=True,generate_H2=False,cutoff_fracs=(a,b+c,d)))]
    cases.extend([
        ('Joint_hydroprocess',joint,dict(common,include_PSA=True,include_EC=False,cutoff_fracs=(a,b,c,d),intermediate_surrogate_ratios={'C18H38':1.0},intermediate_route='hydroprocess',biobinder_price=BINDER_PRICE,use_RF=True)),
        ('Binder_esterification',binder,dict(common,skip_EC=True,generate_H2=False,cutoff_fracs=(a,b+c,d))),
        ('Joint_marine',joint,dict(common,include_PSA=True,include_EC=False,cutoff_fracs=(a,b,c,d),intermediate_surrogate_ratios={'C18H38':1.0},intermediate_route='marine',biobinder_price=BINDER_PRICE,use_RF=True)),
        ('Refinery_CHP',binder,dict(common,skip_EC=True,generate_H2=False,cutoff_fracs=(a,b+c,d),heavy_route='chp')),
    ])
    if INCLUDE_SECOND_REFINERY_ROUTE:
        cases.append(('Refinery_binder',binder,dict(common,skip_EC=True,generate_H2=False,cutoff_fracs=(a,b+c,d),heavy_route='binder')))
    snapshots=[];basis=None
    for name,module,kwargs in cases:
        capture=io.StringIO();system=None;error='';settings=dict(Feedstock=FEEDSTOCK,
            Decentralized_HTL=DECENTRALIZED_HTL,Four_cut_basis=str(FOUR_CUTS),
            Native_cut_fractions=str(kwargs['cutoff_fracs']),Intermediate_route=kwargs.get('intermediate_route',kwargs.get('heavy_route','not applicable')))
        print(f'Running {name} ...')
        try:
            with contextlib.redirect_stdout(capture),contextlib.redirect_stderr(capture),warnings.catch_warnings():
                warnings.simplefilter('always')
                if basis is not None:
                    kwargs['HTL_yields_override']=dict(basis['yields'])
                    kwargs['biocrude_ratios']=dict(basis['component_ratios'])
                    kwargs['central_dry_flowrate']=basis['central_dry_kg_hr']
                    kwargs['pilot_dry_flowrate']=basis['pilot_dry_kg_hr']
                system=module.create_system(**kwargs)
                if basis is None:
                    splitter=system.flowsheet.unit.CrudeSplitter
                    weights={ID:float(f)*float(w) for f,group in zip(splitter.cutoff_fracs,splitter._joint_bins) for ID,w in group.items()}
                    HTL=system.flowsheet.unit.HTL;conditioning=system.flowsheet.unit.FeedstockCond
                    pilot=float(conditioning.feedstock_dry_flowrate)
                    n=float(getattr(HTL,'N_unit',1))
                    keys=('IRR','duration','income_tax','finance_interest','finance_years','finance_fraction',
                          'WC_over_FCI','construction_schedule','startup_months','startup_FOCfrac','startup_VOCfrac','startup_salesfrac')
                    inputs={k:getattr(system.TEA,k) for k in keys if hasattr(system.TEA,k)}
                    inputs.update(COMMON_TEA_OVERRIDES)
                    price=system.flowsheet.stream.scaled_feedstock.price if FEEDSTOCK_PRICE_WET_TONNE is None else FEEDSTOCK_PRICE_WET_TONNE/1000
                    basis=dict(yields=dict(HTL.dw_yields),component_ratios=weights,feed_price_USD_kg=float(price),
                        hours=float(system.operating_hours),tea_inputs=inputs,central_dry_kg_hr=pilot*n,pilot_dry_kg_hr=pilot)
                    (run_dir/'common_basis.json').write_text(json.dumps(basis,indent=2,default=str))
                configure(system,basis)
                settings.update(HTL_yields=str(basis['yields']),Feedstock_price_USD_wet_tonne=basis['feed_price_USD_kg']*1000,
                    Central_dry_kg_hr=basis['central_dry_kg_hr'],Pilot_dry_kg_hr=basis['pilot_dry_kg_hr'],
                    N_HTL=getattr(system.flowsheet.unit.HTL,'N_unit',1))
                system.simulate()
                settings.update(apply_price(system,name,binder._HHV_per_GGE))
        except Exception:
            error=traceback.format_exc();capture.write('\n'+error)
        logpath=run_dir/(name+'_console.log');logpath.write_text(capture.getvalue(),encoding='utf-8')
        if system is not None:
            with contextlib.redirect_stdout(capture),contextlib.redirect_stderr(capture):
                try:
                    snap=collect(system,name,settings,simulation_error=error,console_log=str(logpath))
                    snap['metrics'].update({k:v for k,v in settings.items() if k in ('Price_basis','Priced_stream','Sale_price_USD_kg','Sale_mass_kg_year','Esterification_USD_kg','External_processing_deduction_USD_year','PNNL_USD_GGE','VLSFO_reference_date','Intermediate_route')})
                except Exception:
                    reporting_error=traceback.format_exc();capture.write(reporting_error)
                    snap=failed_snapshot(name,settings,'Reporting failed after system build: '+reporting_error,str(logpath))
            logpath.write_text(capture.getvalue(),encoding='utf-8')
        else:snap=failed_snapshot(name,settings,error,str(logpath))
        snapshots.append(snap)
        print(f"{name}: {snap['metrics'].get('Status')} | IRR {snap['metrics'].get('IRR_pct')} %")
        if basis is None:
            # Other systems cannot meet a common RF/component basis if the joint build failed.
            for other,_,kw in cases[1:]:snapshots.append(failed_snapshot(other,{},'Skipped: common basis could not be constructed',''))
            break
    notes=[('Scenario pricing','External esterification is a price netback only; no onsite equipment or LCA added. Marine historical 2022 benchmark is editable and unnormalized by default. PNNL whole-crude BEV is a sensitivity proxy for lighter cut, not a measured fraction purchase price.'),('Common assay','Four fractions are sensitivity assumptions, not measured assay. Joint uses C18H38 in the 300-343 bin.'),
        ('Matched composition','Joint within-bin weights multiplied by four fractions become common component ratios; coarse binder/SAF bins preserve the same surrogate composition.'),
        ('Routes','SAF: >300 to CHP. Binder: 150-343 biofuel, >343 binder. Joint: 150-300 upgrading; intermediate selected route; >343 binder.'),
        ('Boundary','SAF/joint aggregate gas, char and aqueous output. Binder retains local HTL EC/pass-through and aggregate wastewater scaler. Logistics for aggregate collection are not added.'),
        ('Prices','Fuel prices remain native market assumptions. Binder fixed at configured price. No solve_price/MSP calls; sale intermediate price is explicit.'),
        ('Accounting','Feed stream cost positive=cost, negative=tipping credit. Product value positive=revenue, negative=disposal expense. Gross revenues are separated from net system sales.'),
        ('LCA comparison','Native indicators reported annually and per wet tonne of feed. No invented product allocation, binder displacement or common TRACI factors; missing categories are unavailable.'),
        ('LCA qualification','Inherited LCA definitions/factors are reported, not revalidated. Compare CFs and boundaries before interpreting cross-pathway differences.'),
        ('Capital','Unit/equipment costs are native totals; do not multiply by N_unit again. Land basis aligned to aggregate crude volume.'),
        ('Distillation','FINAL snapshot and logged TRIAL rows share one sheet. VLE diagnostic evaluations do not change chosen LHK/recoveries. Helper may retain only latest invocation.'),
        ('Failures','FAIL snapshots are incomplete/stale; financial totals from failed process runs are diagnostic only. IRR suppressed for failed runs.'),
        ('Validation','IRR root checked against cash-flow NPV within USD1000. Comparison has not been run in the preparation environment.'),
        ('Source templates','User TEA/LCA workbook guided native cashflow, utility, stream-cost and impact-table reporting.')]
    workbook=export(snapshots,run_dir/f'pathway_comparison_{stamp}.xlsx',notes)
    print(f'Workbook: {workbook}')
    return snapshots,workbook

if __name__=='__main__':
    snapshots,workbook=run_comparison()
