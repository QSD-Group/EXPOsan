# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Render saved paired-study observations without constructing or rerunning models."""
import argparse
import csv
from datetime import datetime
import json
from math import isfinite
import os
from pathlib import Path
from zoneinfo import ZoneInfo

DEPLOY = Path(__file__).resolve().parent
RESULTS = DEPLOY/'results'
LABELS = {'htl':'HTL', 'ad':'AD–CHP', 'ad_rng':'AD–RNG',
          'composting':'Composting', 'htg_methane':'HTG methane',
          'htg_hydrogen':'HTG hydrogen', 'htc':'HTC',
          'fast_pyrolysis':'Fast pyrolysis'}
METRICS = ('NPV_usd','IRR_percent','FCI_usd','annual_opex_usd',
           'annual_revenue_usd','GWP_kgCO2e_t_wet')
PRODUCT_METRICS = ('mass_kg_h','hhv_mj_kg','water_wt','ash_dry_wt',
    'dry_matter_wt',
    'carbon_dry_wt','hydrogen_dry_wt','oxygen_dry_wt','nitrogen_dry_wt',
    'sulfur_dry_wt','purity_mol','price_usd_kg','annual_revenue_usd',
    'marketed_kg_year','carbon_nitrogen_ratio','organic_nitrogen_kg_h',
    'ammonium_n_kg_h','nitrate_n_kg_h','normal_volume_nm3_h',
    'gas_CH4_dry_mol','gas_H2_dry_mol','gas_CO2_dry_mol','gas_CO_dry_mol',
    'gas_NH3_dry_mol','gas_H2S_dry_mol')
DRY_FRACTIONS = ('Carbohydrates','Proteins','Ash','Lipids')
FINANCIAL_OUTCOMES = ('NPV_usd','IRR_percent','FCI_usd',
    'reconciled_cash_cost_usd_year','reconciled_receipts_usd_year')


def number(value):
    """Return a finite real observation, retaining missing values as None."""
    return float(value) if isinstance(value,(int,float)) and not isinstance(value,bool) and isfinite(value) else None


def load_json(path, default=None):
    return json.loads(path.read_text(encoding='utf-8-sig')) if path.exists() else default


def block_pairs(metrics, prefix):
    """Read a reporter's two-column cash ledger, preserving duplicate row names."""
    for title, block in metrics.get('reporting_blocks',{}).items():
        if title.startswith(prefix) and block.get('columns') == ['Item','Value']:
            return [(str(row[0]),number(row[1])) for row in block.get('data',[])
                    if len(row)>1 and number(row[1]) is not None]
    return None


def cash_ledger(row):
    """Map retained costs/receipts while preserving native signed HTL accounting.

    HTL negative input costs become feed credits; negative product sales become
    residual-management costs. The bridge preserves sales minus operating cost.
    It does not infer an unavailable equipment or utility burden.
    """
    m = row.get('metrics',{})
    if row.get('status') != 'OK':
        return {},{},'unavailable: native case failed'
    if row['technology'] != 'htl':
        costs, revenues = m.get('annual_costs'),m.get('annual_revenues')
        if not isinstance(costs,dict) or not isinstance(revenues,dict):
            return {},{},'unavailable: annual cash ledgers absent'
        return ({k:number(v) for k,v in costs.items() if number(v) is not None},
                {k:number(v) for k,v in revenues.items() if number(v) is not None},
                'saved alternative annual_costs/annual_revenues')
    inputs,outputs = block_pairs(m,'FEEDS'),block_pairs(m,'PRODUCTS')
    fixed,utilities = number(m.get('tea.FOC')),number(m.get('sys.utility_cost'))
    material,sales,aoc = (number(m.get(k)) for k in ('sys.material_cost','sys.sales','tea.AOC'))
    if inputs is None or outputs is None or any(v is None for v in (fixed,utilities,material,sales,aoc)):
        return {},{},'unavailable: HTL cash bridge operands absent'
    tolerance = max(0.01,abs(aoc)*1e-8,abs(sales)*1e-8)
    if abs(sum(v for _,v in inputs)-material)>tolerance or abs(sum(v for _,v in outputs)-sales)>tolerance:
        return {},{},'unavailable: HTL saved feed/product ledgers do not reconcile'
    costs = {'native_fixed_operations':fixed,'native_net_utilities':utilities}
    revenues = {}
    for direction, entries in (('input',inputs),('output',outputs)):
        for index,(name,value) in enumerate(entries):
            key = f'{direction}:{name}:{index}'
            if direction == 'input':
                (costs if value>=0 else revenues)[key] = abs(value)
            else:
                (revenues if value>=0 else costs)[key] = abs(value)
    residual = sum(revenues.values())-sum(costs.values())-(sales-aoc)
    if abs(residual)>tolerance:
        return {},{},f'unavailable: HTL cash bridge residual {residual:.6g} USD/year'
    return costs,revenues,'reconciled HTL saved cash-ledger bridge; native net utilities retained signed'


def cost_group(name):
    if name in ('maintenance','insurance','labor','land_lease','native_fixed_operations'):
        return 'Fixed operations'
    if name in ('electricity','heat','cooling','diesel','native_net_utilities') or 'natural_gas' in name:
        return 'Energy (net where native)'
    if name in ('residue','liquid') or name.startswith('output:'):
        return 'Residual management'
    if name=='water' or 'process_water' in name.lower():
        return 'Water'
    if name=='other':
        return 'Other operations'
    return 'Materials and handling'


def revenue_group(name):
    if name=='gate_fee' or name.startswith('input:'):
        return 'Gate fees / input credits'
    if 'electricity' in name:
        return 'Electricity sales'
    if name=='heat' or 'heat_export' in name:
        return 'Heat sales'
    return 'Product sales'


def case_observations(row):
    """Extract saved physical and economic observations with named cash bases."""
    m = row.get('metrics',{})
    costs,revenues,basis = cash_ledger(row)
    hours,tonnes = number(m.get('operating_hours')),number(m.get('annual_feed_tonnes'))
    inv = m.get('annual_inventory',{})
    heating = number(inv.get('heat_import_mj'))
    power = number(m.get('net_electricity_kw'))
    result = dict(sample_id=row['sample_id'],technology=row['technology'],
        status=row.get('status','UNOBSERVED'),IRR_status=m.get('IRR_status',''),
        cash_ledger_basis=basis,annual_feed_tonnes=tonnes,
        operating_hours=hours,feedstock_price_usd_t=number(m.get('feedstock_price_usd_t')),
        native_HTL_AOC_usd_year=number(m.get('tea.AOC')),
        native_HTL_sales_usd_year=number(m.get('sys.sales')),
        reconciled_cash_cost_usd_year=sum(costs.values()) if costs else None,
        reconciled_receipts_usd_year=sum(revenues.values()) if revenues else None,
        purchased_heat_mj_t_wet=heating/tonnes if heating is not None and tonnes else None,
        net_electricity_kwh_t_wet=power*hours/tonnes if power is not None and hours and tonnes else None,
        heat_recovered_mj_h=number(m.get('design',{}).get('heat_recovered_mj_h')),
        LCA_status=m.get('LCA_status','unavailable'))
    for key in (*METRICS,'TCI_usd','EBITDA_usd'):
        result[key] = number(m.get(key)) if row.get('status')=='OK' else None
    for key in ('raw_biogas_nm3_h','methane_nm3_h','raw_methane_volume_fraction',
                'retained_ammonia_kg_h'):
        result[key] = number(m.get('design',{}).get(key))
    return result,costs,revenues


def product_observations(row):
    """Retain product mass, quality and pricing without assigning fuel grades."""
    m,tech = row.get('metrics',{}),row['technology']
    if row.get('status')!='OK':
        return []
    quality = m.get('product_quality',m.get('design',{}).get('product_quality',{}))
    markets = m.get('product_markets',m.get('design',{}).get('product_markets',{}))
    records = []
    for name,q in quality.items():
        if not isinstance(q,dict) or number(q.get('mass_kg_h')) is None:
            continue
        market = markets.get(name,{})
        price = market.get('price_basis',{})
        dry = q.get('ultimate_dry_wt',{})
        records.append(dict(product=name,stream_role='main_product',mass_basis='gross whole product',
            mass_kg_h=q.get('mass_kg_h'),hhv_mj_kg=q.get('hhv_mj_kg'),
            water_wt=q.get('water_wt'),ash_dry_wt=q.get('ash_dry_wt'),
            carbon_dry_wt=dry.get('C'),hydrogen_dry_wt=dry.get('H'),oxygen_dry_wt=dry.get('O'),
            nitrogen_dry_wt=dry.get('N'),sulfur_dry_wt=dry.get('S'),
            purity_mol=q.get('purity_mol'),market_status=market.get('status','UNPRICED'),
            market_basis=price.get('basis',''),price_source=price.get('source',''),
            price_basis_json=json.dumps(price,sort_keys=True),
            quality_json=json.dumps(q,sort_keys=True),
            price_usd_kg=m.get('product_prices_usd_kg',{}).get(name),
            annual_revenue_usd=market.get('revenue_usd_year'),
            marketed_kg_year=market.get('marketed_kg_year'),
            normal_volume_nm3_h=q.get('normal_volume_nm3_h'),
            **{'gas_'+k+'_dry_mol':v for k,v in q.get('dry_gas_mol',{}).items()}))
    if tech=='htl':
        for name,mass in m.get('product_kg_h',{}).items():
            price = number(m.get('product_prices_usd_kg',{}).get(name))
            if price is not None and price>0 and number(mass) is not None and mass>0:
                records.append(dict(product=name,stream_role='main_product',mass_basis='positive-price native product',
                    mass_kg_h=mass,price_usd_kg=price,
                    market_basis='HTL configured product pricing/upgrading quality assumption',
                    market_status='NATIVE_PRICE_SCENARIO'))
    elif tech=='composting' and number(m.get('marketable_compost_wet_kg_h')) is not None:
        q = m.get('design',{}).get('product_quality',{})
        records.append(dict(product='compost',stream_role='main_product',mass_basis='marketable wet compost',
            mass_kg_h=m['marketable_compost_wet_kg_h'],
            water_wt=m.get('compost_moisture_wet_mass_fraction'),
            carbon_dry_wt=q.get('carbon_dry_fraction'),nitrogen_dry_wt=q.get('nitrogen_dry_fraction'),
            carbon_nitrogen_ratio=q.get('carbon_nitrogen_ratio'),
            organic_nitrogen_kg_h=q.get('organic_nitrogen_kg_h'),
            ammonium_n_kg_h=q.get('ammonium_n_kg_h'),nitrate_n_kg_h=q.get('nitrate_n_kg_h'),
            market_basis='configured compost bulk-price and marketability scenario',
            market_status='PRICE_SCENARIO'))
    if tech in ('ad','ad_rng'):
        d = m.get('design',{})
        methane = number(d.get('raw_methane_volume_fraction'))
        gas = m.get('process_streams_kg_h',{}).get('raw_biogas',{})
        records.append(dict(product='raw_biogas',stream_role='intermediate',mass_basis='raw intermediate; CH4/CO2 model gas',
            mass_kg_h=sum(gas.values()) if gas else None,component_kg_h_json=json.dumps(gas,sort_keys=True),
            normal_volume_nm3_h=d.get('raw_biogas_nm3_h'),gas_CH4_dry_mol=methane,
            gas_CO2_dry_mol=1-methane if methane is not None else None,
            market_basis='raw gas before CHP combustion or RNG purification',market_status='INTERMEDIATE'))
        for name in ('digestate_cake','centrate'):
            composition = m.get('process_streams_kg_h',{}).get(name,{})
            total = sum(composition.values())
            if total:
                records.append(dict(product=name,stream_role='intermediate' if name=='digestate_cake' else 'residual',
                    mass_basis='gross digestate cake before marketable split' if name=='digestate_cake' else 'liquid digestate residual',
                    mass_kg_h=total,water_wt=composition.get('Water',0.)/total,
                    dry_matter_wt=composition.get('DryMatter',0.)/total,
                    component_kg_h_json=json.dumps(composition,sort_keys=True),
                    composition_basis='aggregate Water/DryMatter balance; no cake/centrate VS, ash or nitrogen-form partition inferred',
                    market_status='INTERMEDIATE' if name=='digestate_cake' else 'RESIDUAL'))
    process_names = set(m.get('process_streams_kg_h',{}))
    process_names.update(('digestate_product','compost_product','residue','combustion_water'))
    for native_name,q in m.get('product_stream_observations',{}).items():
        if number(q.get('mass_kg_h')) is not None and q['mass_kg_h']<=0:
            continue
        candidates = process_names | {r['product'] for r in records}
        name = next((k for k in sorted(candidates,key=len,reverse=True)
                     if native_name==k or native_name.endswith('_'+k)),native_name)
        if name=='compost_product':
            name = 'compost'
        existing = next((r for r in records if r['product']==name),None)
        if existing is None:
            main = ((tech in ('ad','ad_rng') and name in ('digestate_product','rng')) or
                    (tech=='htl' and name in ('biofuel','biobinder','biocrude','recovered_H2','recovered_N','recovered_P','recovered_K')))
            existing = dict(product=name,stream_role='main_product' if main else 'residual',
                mass_basis='native boundary product stream' if main else 'native residual/emission boundary stream',
                market_status='NATIVE_SCENARIO' if main else 'RESIDUAL')
            records.append(existing)
        existing['native_stream_id'] = native_name
        for key in ('mass_kg_h','hhv_mj_kg','ash_dry_wt','purity_mol'):
            if number(q.get(key)) is not None:
                existing[key] = q[key]
        mass,water = number(q.get('mass_kg_h')),number(q.get('water_kg_h'))
        if mass and water is not None:
            existing['water_wt'] = water/mass
        for element,key in zip(('C','H','O','N','S'),
            ('carbon_dry_wt','hydrogen_dry_wt','oxygen_dry_wt','nitrogen_dry_wt','sulfur_dry_wt')):
            if number(q.get('ultimate_dry_wt',{}).get(element)) is not None:
                existing[key] = q['ultimate_dry_wt'][element]
        existing['component_kg_h_json'] = json.dumps(q.get('component_kg_h',{}),sort_keys=True)
        existing['chemical_identities_json'] = json.dumps(q.get('chemical_identities',{}),sort_keys=True)
        existing['composition_basis'] = q.get('composition_basis','native component lumps; identity does not establish experimental speciation')
        if tech in ('ad','ad_rng','composting') and 'DryMatter' in q.get('component_kg_h',{}):
            existing['native_placeholder_hhv_mj_kg'] = q.get('hhv_mj_kg')
            existing['hhv_mj_kg'] = None
            existing['dry_matter_wt'] = q['component_kg_h']['DryMatter']/mass if mass else None
            existing['hhv_basis'] = 'unavailable: formula-free DryMatter adapter placeholder is not a fuel characterization'
    for record in records:
        if record.get('stream_role')!='main_product':
            continue
        if tech=='htl':
            hours,mass,price = number(m.get('operating_hours')),number(record.get('mass_kg_h')),number(record.get('price_usd_kg'))
            if hours is not None and mass is not None and price is not None:
                record['marketed_kg_year'] = mass*hours
                record['annual_revenue_usd'] = mass*hours*price
        elif tech in ('ad','ad_rng','composting'):
            item = 'compost' if tech=='composting' else 'rng' if record['product']=='rng' else 'digestate'
            quantity = number(m.get('design',{}).get('marketable_product_kg_year'))
            if item=='rng':
                mass,hours = number(record.get('mass_kg_h')),number(m.get('operating_hours'))
                quantity = mass*hours if mass is not None and hours is not None else None
            revenue = number(m.get('annual_revenues',{}).get(item))
            record['marketed_kg_year'] = quantity
            record['annual_revenue_usd'] = revenue
            if quantity and revenue is not None:
                record['price_usd_kg'] = revenue/quantity
                record['price_basis_detail'] = 'effective price from saved product receipt divided by saved marketable mass'
            if tech in ('ad','ad_rng'):
                record['market_basis'] = 'configured '+item+' bulk-price and marketability scenario'
    return [dict(sample_id=row['sample_id'],technology=tech,**record) for record in records]


def variable(values):
    """Distinguish resolved numerical variation from constant or rounding-only outputs."""
    if len(values)<2:
        return False
    scale = max(max(abs(v) for v in values),1e-12)
    return max(values)-min(values)>1e-10*scale


def chemistry_analysis(data):
    """Summarize retained product distributions and paired descriptive associations."""
    import numpy as np
    from scipy.stats import spearmanr
    products = [r for r in data['products'] if r['population']=='sampled' and r.get('stream_role')!='residual']
    summaries,components,associations = [],[],[]
    for tech,product in dict.fromkeys((r['technology'],r['product']) for r in products):
        selected = [r for r in products if r['technology']==tech and r['product']==product]
        planned = next(r['planned'] for r in data['coverage'] if r['technology']==tech)
        for metric in PRODUCT_METRICS:
            values = [number(r.get(metric)) for r in selected]
            finite = [v for v in values if v is not None]
            item = dict(technology=tech,product=product,stream_role=selected[0].get('stream_role',''),metric=metric,planned=planned,
                count=len(finite),unavailable_or_failed=planned-len(finite),
                variation_status='varies across sampled cases' if variable(finite) else
                'constant within 1e-10 relative numerical tolerance' if finite else 'unavailable')
            if finite:
                item.update(min=min(finite),max=max(finite),mean=float(np.mean(finite)),
                    p05=float(np.quantile(finite,.05)),median=float(np.median(finite)),
                    p95=float(np.quantile(finite,.95)))
            summaries.append(item)
    for p in data['products']:
        masses = json.loads(p.get('component_kg_h_json','{}'))
        identities = json.loads(p.get('chemical_identities_json','{}'))
        total = number(p.get('mass_kg_h'))
        for name,mass in masses.items():
            identity = identities.get(name,{})
            components.append(dict(population=p['population'],sample_id=p['sample_id'],
                technology=p['technology'],product=p['product'],stream_role=p.get('stream_role',''),component=name,
                mass_kg_h=mass,whole_product_wt=mass/total if total else None,
                formula=identity.get('formula',''),CAS=identity.get('CAS',''),
                common_name=identity.get('common_name',''),basis=p.get('composition_basis','')))
    samples = {}
    path = data['run_dir']/'samples.jsonl'
    if path.exists():
        for line in path.read_text(encoding='utf-8-sig').splitlines():
            if line.strip():
                row = json.loads(line)
                if row['sample_id'] in samples:
                    raise ValueError('Duplicate sample identities in the saved design')
                samples[row['sample_id']] = row
    cases = {(r['technology'],r['sample_id']):r for r in data['observations']
             if r['population']=='sampled' and r['status']=='OK'}
    def record(tech,link,xname,yname,joined):
        finite = [(i,number(x),number(y)) for i,x,y in joined if number(x) is not None and number(y) is not None]
        item = dict(technology=tech,link=link,input=xname,outcome=yname,count=len(finite),
            sample_ids_json=json.dumps([i for i,_,_ in finite]),rho=None,
            status='insufficient paired observations')
        if len(finite)>=3:
            xs,ys = [r[1] for r in finite],[r[2] for r in finite]
            if not variable(xs):
                item['status'] = 'constant input; correlation undefined'
            elif not variable(ys):
                item['status'] = 'constant outcome; correlation undefined'
            else:
                rho = number(float(spearmanr(xs,ys).statistic))
                item.update(rho=rho,status='descriptive paired Spearman association' if rho is not None else 'undefined')
        associations.append(item)
    for tech in data['techs']:
        route = [(sid,c) for (t,sid),c in cases.items() if t==tech]
        for fraction in DRY_FRACTIONS:
            xname = 'feed_dry.'+fraction
            def feed(sid):
                return samples.get(sid,{}).get('dry_composition',{}).get(fraction)
            for metric in (*FINANCIAL_OUTCOMES,'purchased_heat_mj_t_wet','net_electricity_kwh_t_wet',
                           'methane_nm3_h','raw_methane_volume_fraction','retained_ammonia_kg_h'):
                record(tech,'feed_to_process_or_finance',xname,metric,
                    [(sid,feed(sid),c.get(metric)) for sid,c in route])
            for product in dict.fromkeys(p['product'] for p in products if p['technology']==tech):
                selected = [p for p in products if p['technology']==tech and p['product']==product]
                for metric in PRODUCT_METRICS:
                    if any(number(p.get(metric)) is not None for p in selected):
                        record(tech,'feed_to_product',xname,product+'.'+metric,
                            [(p['sample_id'],feed(p['sample_id']),p.get(metric)) for p in selected])
        for product in dict.fromkeys(p['product'] for p in products if p['technology']==tech):
            selected = [p for p in products if p['technology']==tech and p['product']==product]
            for metric in ('mass_kg_h','hhv_mj_kg','nitrogen_dry_wt','ash_dry_wt','water_wt',
                           'purity_mol','price_usd_kg','annual_revenue_usd','gas_CH4_dry_mol',
                           'carbon_nitrogen_ratio'):
                if any(number(p.get(metric)) is not None for p in selected):
                    for outcome in FINANCIAL_OUTCOMES:
                        record(tech,'product_to_finance',product+'.'+metric,outcome,
                            [(p['sample_id'],p.get(metric),cases.get((tech,p['sample_id']),{}).get(outcome)) for p in selected])
    return summaries,components,associations


def cash_summary(data):
    """Summarize the additional reconciled cash presentation absent from summary.csv."""
    import numpy as np
    rows = []
    for tech in data['techs']:
        selected = [r for r in data['observations'] if r['technology']==tech and r['population']=='sampled']
        for metric in ('reconciled_cash_cost_usd_year','reconciled_receipts_usd_year'):
            values = [number(r.get(metric)) for r in selected]
            finite = [v for v in values if v is not None]
            item = dict(technology=tech,metric=metric,count=len(finite),
                unavailable_or_failed=data['config']['n_samples']-len(finite),
                basis='HTL saved-ledger bridge / alternative annual ledgers; native net utilities retain their sign')
            if finite:
                item.update(p05=float(np.quantile(finite,.05)),median=float(np.median(finite)),
                    p95=float(np.quantile(finite,.95)))
            rows.append(item)
    return rows


def prepare(run_dir, allow_incomplete=False):
    """Load saved populations and validate that exported summaries match them."""
    run_dir = Path(run_dir).resolve()
    if not run_dir.is_relative_to(RESULTS.resolve()) or run_dir==RESULTS.resolve():
        raise ValueError('Select one run directory inside deploy/results')
    manifest = load_json(run_dir/'manifest.json')
    if not manifest:
        raise ValueError('Run manifest is missing')
    status = load_json(run_dir/'run_status.json',{})
    if status.get('status')!='COMPLETE' and not allow_incomplete:
        raise ValueError('Run is not COMPLETE; use --allow-incomplete only for an explicitly partial report')
    checkpoint = run_dir/'results.jsonl'
    if not checkpoint.exists():
        raise ValueError('No sampled results.jsonl; a preflight baseline is not an LHS population')
    rows = [json.loads(line) for line in checkpoint.read_text(encoding='utf-8-sig').splitlines() if line.strip()]
    if not rows or any(not isinstance(row.get('sample_id'),int) or row['sample_id']<0 for row in rows):
        raise ValueError('Sample checkpoint must contain positive or zero sample identities, separate from baseline')
    identities = [(r['sample_id'],r['technology']) for r in rows]
    if len(identities)!=len(set(identities)):
        raise ValueError('Duplicate sample/technology identities; resolve checkpoint provenance before analysis')
    config = manifest['config']
    techs = config.get('technologies') or [k for k in LABELS if any(r['technology']==k for r in rows)]
    if set(techs)-set(LABELS) or set(r['technology'] for r in rows)-set(techs):
        raise ValueError('Unsupported or unscheduled technology in saved population')
    baseline = load_json(run_dir/'baseline_results.json',[])
    summary_path = run_dir/'summary.csv'
    if not summary_path.exists():
        raise ValueError('Saved summary.csv is required; export the retained run before analysis')
    with summary_path.open(encoding='utf-8-sig',newline='') as stream:
        summary = list(csv.DictReader(stream))
    for item in summary:
        tech,metric = item['technology'],item['metric']
        selected = [r for r in rows if r['technology']==tech]
        if metric=='status':
            if int(float(item['attempted']))!=len(selected):
                raise ValueError('summary.csv population is stale or includes the baseline')
        else:
            count = sum(r.get('status')=='OK' and number(r.get('metrics',{}).get(metric)) is not None for r in selected)
            if int(float(item.get('count') or 0))!=count:
                raise ValueError(f'summary.csv count disagrees with checkpoint: {tech}/{metric}')
    planned = int(config['n_samples'])
    coverage,metric_coverage,lca_coverage = [],[],[]
    observations,costs,revenues,products,equipment,reasons = [],[],[],[],[],[]
    for row in [*baseline,*rows]:
        observation,cs,rs = case_observations(row)
        population = 'baseline' if row['sample_id']<0 else 'sampled'
        observations.append(dict(population=population,**observation))
        for entries,target,grouping in ((cs,costs,cost_group),(rs,revenues,revenue_group)):
            for name,value in entries.items():
                target.append(dict(population=population,sample_id=row['sample_id'],
                    technology=row['technology'],item=name,group=grouping(name),usd_year=value,
                    basis=observation['cash_ledger_basis']))
        products.extend(dict(population=population,**p) for p in product_observations(row))
        for item in row.get('metrics',{}).get('equipment_audit',[]):
            installed = number(item.get('installed_usd'))
            if installed is not None:
                equipment.append(dict(population=population,sample_id=row['sample_id'],
                    technology=row['technology'],equipment=item.get('name',item.get('unit_id','')),
                    unit_class=item.get('unit_class',''),installed_usd=installed,
                    included_in_system_capital=item.get('included_in_system_capital',True),
                    cost_area=item.get('cost_area',''),source=item.get('source','')))
        if row.get('status')!='OK' or number(row.get('metrics',{}).get('IRR_percent')) is None:
            reasons.append(dict(population=population,sample_id=row['sample_id'],technology=row['technology'],
                case_status=row.get('status'),IRR_status=row.get('metrics',{}).get('IRR_status','unobserved'),
                error=row.get('error',''),log=row.get('log','')))
    for tech in techs:
        selected = [r for r in rows if r['technology']==tech]
        ok = [r for r in selected if r.get('status')=='OK']
        finite_irr = sum(number(r['metrics'].get('IRR_percent')) is not None for r in ok)
        coverage.append(dict(technology=tech,planned=planned,attempted=len(selected),
            succeeded=len(ok),failed_or_timeout=len(selected)-len(ok),
            unobserved=max(planned-len(selected),0),finite_IRR=finite_irr,
            successful_case_IRR_unavailable=len(ok)-finite_irr,
            finite_NPV=sum(number(r['metrics'].get('NPV_usd')) is not None for r in ok),
            baseline_status=next((r.get('status') for r in baseline if r['technology']==tech),'UNOBSERVED')))
        for metric in METRICS:
            finite = sum(number(r['metrics'].get(metric)) is not None for r in ok)
            metric_coverage.append(dict(technology=tech,metric=metric,planned=planned,
                finite=finite,unavailable_in_successful_case=len(ok)-finite,
                native_failed_or_timeout=len(selected)-len(ok),unobserved=max(planned-len(selected),0)))
        for state in dict.fromkeys(r.get('metrics',{}).get('LCA_status','unavailable') for r in selected):
            lca_coverage.append(dict(technology=tech,LCA_status=state,
                count=sum(r.get('metrics',{}).get('LCA_status','unavailable')==state for r in selected)))
    data = dict(run_dir=run_dir,manifest=manifest,config=config,status=status,techs=list(techs),
        rows=rows,baseline=baseline,summary=summary,coverage=coverage,observations=observations,
        costs=costs,revenues=revenues,products=products,equipment=equipment,reasons=reasons,
        metric_coverage=metric_coverage,lca_coverage=lca_coverage,
        property_basis=load_json(run_dir/'property_basis.json',{}))
    data['product_summary'],data['components'],data['associations'] = chemistry_analysis(data)
    data['cash_summary'] = cash_summary(data)
    return data


def write_csv(path, rows, columns=None):
    columns = columns or list(dict.fromkeys(key for row in rows for key in row))
    with path.open('w',encoding='utf-8-sig',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=columns or ['status'])
        writer.writeheader()
        writer.writerows(rows)


def fmt(value, digits=2):
    n = number(value)
    return f'{n:,.{digits}f}' if n is not None else 'unavailable'


def markdown_table(headers, rows):
    def escape(v):
        return str(v).replace('|','/').replace('\n',' ')
    return '\n'.join(['| '+' | '.join(headers)+' |',
                      '| '+' | '.join(['---']*len(headers))+' |',
                      *['| '+' | '.join(escape(v) for v in row)+' |' for row in rows]])


def render_figures(data, output):
    """Render saved summaries and baseline ledgers with explicit missing coverage."""
    os.environ.setdefault('MPLCONFIGDIR',str(DEPLOY/'.cache'/'comparison_analysis_mpl'))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
        'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
    techs = data['techs']
    labels = [LABELS[k] for k in techs]
    baseline = {r['technology']:r for r in data['observations'] if r['population']=='baseline'}
    footer = 'Conditional screening scenarios; missing outcomes retained. Values do not establish technology rankings.'
    figure_paths = []
    fig,axes = plt.subplots(1,3,figsize=(13,5.4),sharey=True)
    for ax,(metric,title,scale,unit) in zip(axes,(
        ('NPV_usd','Net present value',1e6,'million USD'),
        ('IRR_percent','Internal rate of return',1.,'%'),
        ('FCI_usd','Fixed capital investment',1e6,'million USD'))):
        for j,tech in enumerate(techs):
            item = next((s for s in data['summary'] if s['technology']==tech and s['metric']==metric),None)
            if item and all(item.get(k) not in ('',None) for k in ('p05','median','p95')):
                lo,mid,hi = (float(item[k])/scale for k in ('p05','median','p95'))
                ax.plot([lo,hi],[j,j],color='#335f86',lw=2.5)
                ax.plot(mid,j,'o',color='#335f86',ms=5)
                ax.annotate(f"n={int(float(item['count']))}",(hi,j),xytext=(6,0),
                            textcoords='offset points',va='center',fontsize=8)
            else:
                ax.text(.02,j,'No finite sampled value',transform=ax.get_yaxis_transform(),va='center',fontsize=8)
            value = number(baseline.get(tech,{}).get(metric))
            if value is not None:
                ax.plot(value/scale,j,'D',color='#d77d2b',ms=5)
        ax.axvline(0,color='.7',lw=.8)
        ax.set_title(title)
        ax.set_xlabel(unit)
        ax.grid(axis='x',alpha=.18)
        ax.margins(x=.25)
    axes[0].set_yticks(range(len(techs)),labels)
    axes[0].invert_yaxis()
    fig.suptitle('Saved sampled p05–median–p95, with separate baseline diamonds')
    fig.text(.5,.015,footer,ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.045,1,.93))
    filename = '01_financial_spread.png'
    fig.savefig(output/filename)
    plt.close(fig)
    figure_paths.append(filename)
    fig,axes = plt.subplots(1,2,figsize=(12,6.2),sharey=True)
    colors = ['#476f91','#d99150','#75947c','#a187a8','#b1a264','#7e959a']
    for ax,entries,title in zip(axes,(data['costs'],data['revenues']),
                              ('Annual cash costs','Annual receipts')):
        entries = [r for r in entries if r['population']=='baseline']
        groups = list(dict.fromkeys(r['group'] for r in entries))
        pos,neg = np.zeros(len(techs)),np.zeros(len(techs))
        for i,group in enumerate(groups):
            values = np.array([sum(r['usd_year'] for r in entries if r['technology']==tech and r['group']==group)/1e6 for tech in techs])
            left = np.where(values>=0,pos,neg)
            ax.barh(range(len(techs)),values,left=left,color=colors[i%len(colors)],label=group,height=.65)
            pos += np.maximum(values,0.)
            neg += np.minimum(values,0.)
        for j,tech in enumerate(techs):
            key = 'reconciled_cash_cost_usd_year' if title.endswith('costs') else 'reconciled_receipts_usd_year'
            total = number(baseline.get(tech,{}).get(key))
            if total is None:
                ax.text(.02,j,'unavailable',transform=ax.get_yaxis_transform(),va='center',fontsize=8)
        ax.set_title(title+' — baseline')
        ax.set_xlabel('million USD/year')
        ax.grid(axis='x',alpha=.18)
        ax.axvline(0,color='.5',lw=.8)
        if groups:
            ax.legend(loc='upper center',bbox_to_anchor=(.5,-.12),fontsize=8,ncol=2,frameon=False)
    axes[0].set_yticks(range(len(techs)),labels)
    axes[0].invert_yaxis()
    fig.suptitle('Saved cash ledgers: HTL credits and disposal moved to explicit sides')
    fig.text(.5,.015,'Native net utility costs retain their sign. This presentation preserves receipts minus costs.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.12,1,.93))
    filename = '02_baseline_cash_drivers.png'
    fig.savefig(output/filename)
    plt.close(fig)
    figure_paths.append(filename)
    fig,axes = plt.subplots(1,2,figsize=(11,5.4),sharey=True)
    for ax,key,title,unit in (
        (axes[0],'purchased_heat_mj_t_wet','Purchased external heat','MJ / wet metric tonne received'),
        (axes[1],'net_electricity_kwh_t_wet','Net electricity: import (+), export (−)','kWh / wet metric tonne received')):
        for j,tech in enumerate(techs):
            value = number(baseline.get(tech,{}).get(key))
            if value is None:
                ax.text(.02,j,'unavailable',transform=ax.get_yaxis_transform(),va='center',fontsize=8)
            else:
                ax.barh(j,value,color='#476f91',height=.6)
        ax.axvline(0,color='.5',lw=.8)
        ax.set_title(title,fontsize=11)
        ax.set_xlabel(unit)
        ax.grid(axis='x',alpha=.18)
    axes[0].set_yticks(range(len(techs)),labels)
    axes[0].invert_yaxis()
    fig.suptitle('Baseline utility drivers on the common incoming-waste basis')
    fig.text(.5,.015,'HTL purchased heating is unavailable when no physical heat inventory was retained; missing is not zero.',ha='center',fontsize=8)
    fig.tight_layout(rect=(0,.045,1,.93))
    filename = '03_baseline_utilities.png'
    fig.savefig(output/filename)
    plt.close(fig)
    figure_paths.append(filename)
    ps = [r for r in data['products'] if r['population']=='baseline' and
          r.get('stream_role')=='main_product' and number(r.get('mass_kg_h')) is not None]
    if ps:
        fig,axes = plt.subplots(1,2,figsize=(12,max(5.2,len(ps)*.43+1.7)),sharey=True)
        for j,p in enumerate(ps):
            axes[0].barh(j,p['mass_kg_h'],color='#476f91',height=.6)
            hhv = number(p.get('hhv_mj_kg'))
            if hhv is None:
                axes[1].text(.02,j,'not retained',transform=axes[1].get_yaxis_transform(),va='center',fontsize=8)
            else:
                axes[1].barh(j,hhv,color='#75947c',height=.6)
        axes[0].set_yticks(range(len(ps)),[LABELS[p['technology']]+' / '+p['product'] for p in ps])
        axes[0].invert_yaxis()
        axes[0].set_title('Reported product quantity')
        axes[0].set_xlabel('kg/h whole product; see CSV for gross/marketable basis')
        axes[1].set_title('Whole-product HHV where retained')
        axes[1].set_xlabel('MJ/kg; thermochemical values are estimated')
        for ax in axes:
            ax.grid(axis='x',alpha=.18)
        fig.suptitle('Baseline products have different functions and qualification requirements')
        fig.text(.5,.015,'AD–CHP exports are electricity/heat, shown separately. Product mass is not equivalent fuel output.',ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.055,1,.94))
        filename = '04_baseline_products.png'
        fig.savefig(output/filename)
        plt.close(fig)
        figure_paths.append(filename)
    quality_rows = sorted(dict.fromkeys((r['technology'],r['product']) for r in data['product_summary']
        if r['count'] and r.get('stream_role')=='main_product' and r['metric'] in ('hhv_mj_kg','nitrogen_dry_wt','ash_dry_wt')),
        key=lambda item:(techs.index(item[0]),item[1]))
    if quality_rows:
        fig,axes = plt.subplots(1,3,figsize=(13,max(5.2,len(quality_rows)*.45+1.7)),sharey=True)
        for ax,(metric,title,scale) in zip(axes,(
            ('hhv_mj_kg','HHV, MJ/kg whole product',1.),
            ('nitrogen_dry_wt','Nitrogen, wt% dry product',100.),
            ('ash_dry_wt','Ash, wt% dry product',100.))):
            for j,(tech,product) in enumerate(quality_rows):
                item = next((r for r in data['product_summary'] if r['technology']==tech and r['product']==product and r['metric']==metric),None)
                if item and item['count']:
                    lo,mid,hi = (item[k]*scale for k in ('p05','median','p95'))
                    ax.plot([lo,hi],[j,j],color='#335f86',lw=2.5)
                    ax.plot(mid,j,'o',color='#335f86',ms=4)
                else:
                    ax.text(.02,j,'unavailable',transform=ax.get_yaxis_transform(),va='center',fontsize=8)
            ax.set_title(title,fontsize=10)
            ax.grid(axis='x',alpha=.18)
        axes[0].set_yticks(range(len(quality_rows)),[LABELS[t]+' / '+p for t,p in quality_rows])
        axes[0].invert_yaxis()
        fig.suptitle('Saved product-quality p05–median–p95 across sampled incoming compositions')
        fig.text(.5,.015,'Conditional sampled spread; HTL includes unresolved adaptive-cut selection variability. No experimental validation implied.',ha='center',fontsize=8)
        fig.tight_layout(rect=(0,.06,1,.94))
        filename = '05_product_quality_response.png'
        fig.savefig(output/filename)
        plt.close(fig)
        figure_paths.append(filename)
    if data['associations']:
        fig,axes = plt.subplots(1,3,figsize=(12.7,5.9),sharey=True)
        for ax,metric,title in zip(axes,('NPV_usd','IRR_percent','FCI_usd'),('NPV','IRR','FCI')):
            matrix = np.full((len(techs),len(DRY_FRACTIONS)),np.nan)
            for j,tech in enumerate(techs):
                for k,fraction in enumerate(DRY_FRACTIONS):
                    item = next((r for r in data['associations'] if r['technology']==tech and
                        r['link']=='feed_to_process_or_finance' and r['input']=='feed_dry.'+fraction and r['outcome']==metric),None)
                    if item and item['rho'] is not None:
                        matrix[j,k] = item['rho']
            cmap = plt.get_cmap('RdBu_r').copy()
            cmap.set_bad('#eeeeee')
            im = ax.imshow(matrix,cmap=cmap,vmin=-1,vmax=1,aspect='auto')
            for j in range(len(techs)):
                for k in range(len(DRY_FRACTIONS)):
                    value = matrix[j,k]
                    ax.text(k,j,f'{value:.2f}' if np.isfinite(value) else '—',
                        ha='center',va='center',fontsize=8,color='white' if np.isfinite(value) and abs(value)>.65 else '#202020')
            ax.set_xticks(range(len(DRY_FRACTIONS)),['Carbohydrate','Protein','Ash','Lipid'],rotation=45,ha='right')
            ax.set_title(title+' association with feed dry fractions')
        axes[0].set_yticks(range(len(techs)),labels)
        fig.subplots_adjust(left=.13,right=.9,bottom=.22,top=.84,wspace=.12)
        cax = fig.add_axes([.92,.25,.015,.5])
        fig.colorbar(im,cax=cax,label='Spearman ρ')
        fig.suptitle('Paired composition–financial associations; usable counts retained in CSV')
        fig.text(.5,.025,'Compositional closure and deterministic model coupling create associations; these are not controlled causal sensitivities.',ha='center',fontsize=8)
        filename = '06_feed_financial_associations.png'
        fig.savefig(output/filename)
        plt.close(fig)
        figure_paths.append(filename)
    return figure_paths


def report_text(data, figures):
    """Compose a basis-aware explanation of observed outcomes and coverage."""
    config,manifest = data['config'],data['manifest']
    capacity = manifest.get('nominal_capacity',{})
    baseline = {r['technology']:r for r in data['observations'] if r['population']=='baseline'}
    lines = ['# Conditional technology comparison from saved observations','',
        f"Run: `{data['run_dir'].name}`. Analysis generated {datetime.now(ZoneInfo('America/New_York')).isoformat(timespec='seconds')}.",'',
        'These are conditional screening results. Conversion scenarios, property estimates, product outlets and equipment scope differ. The figures describe this configured study; they establish neither a technology ranking nor bankability. No simulations or replacement samples were generated by this analysis.','',
        '## Common basis and population','',
        f"The design contains **{config['n_samples']} paired sampled compositions per route**, plus a separate deterministic baseline. Run state: **{data['status'].get('status','unknown')}**. Baseline rows are excluded from all sampled summary statistics.",'',
        f"Capacity input: {fmt(capacity.get('input_value',config.get('capacity_dry_tpd')))} dry US short tons/day under the legacy uptime convention; hourly conversion divides by uptime {fmt(capacity.get('capacity_uptime_ratio',config.get('capacity_uptime_ratio')))}. Retained nominal feed is {fmt(capacity.get('dry_kg_h'),4)} dry kg/h, {fmt(capacity.get('dry_metric_tpd_operating'),6)} dry metric tonnes/operating day and {fmt(capacity.get('wet_metric_tpd_operating'),6)} wet metric tonnes/operating day. Project operating hours remain the saved shared TEA setting.",'',
        'Wet-mass recipe: '+', '.join(f'{k} {v:.1%}' for k,v in config['recipe_wet_mass'].items())+'. '+
        f"The four dry-composition LHS multipliers use {config.get('variation_mode')} variation ±{config.get('variation')}; moisture and wet blending weights remain fixed. Closure couples the realized fractions. This is a composition study, not a propagated distribution of all technology/cost/price uncertainties.",'',
        markdown_table(['Route','Planned','Attempted','Native OK','Failed/timeout','Unobserved','Finite IRR','IRR unavailable in OK case','Baseline'],
            [[LABELS[r['technology']],r['planned'],r['attempted'],r['succeeded'],r['failed_or_timeout'],r['unobserved'],r['finite_IRR'],r['successful_case_IRR_unavailable'],r['baseline_status']] for r in data['coverage']]),'',
        'A missing IRR in a successful native case is distinct from failed construction/simulation. It is not replaced by zero. Metric-specific usable and unavailable counts are retained in `metric_coverage.csv`; saved diagnostic reasons and case logs are indexed in `status_reasons.csv`. Incomplete/failed outcomes retain their identities; no success quota or resampling is applied.','',
        '## Financial outcomes','',
        'NPV uses the saved common discount rate and financed native cash flows; IRR is the retained native cash-flow solution where finite. FCI is fixed capital, while TCI also includes the model’s working-capital/land convention. Sampled p05/median/p95 values below are reused directly from `summary.csv`; they are empirical spread over the usable cases, not confidence intervals or evidence of stable tails.','']
    stats = []
    for tech in data['techs']:
        for metric,label,scale in (('NPV_usd','NPV, million USD',1e6),('IRR_percent','IRR, %',1.),('FCI_usd','FCI, million USD',1e6)):
            item = next((r for r in data['summary'] if r['technology']==tech and r['metric']==metric),None)
            stats.append([LABELS[tech],label,item.get('count','0') if item else '0',
                *[fmt(float(item[k])/scale) if item and item.get(k) not in ('',None) else 'unavailable' for k in ('p05','median','p95')]])
    lines.extend([markdown_table(['Route','Metric','Finite n','p05','Median','p95'],stats),'',
        'Negative NPV means the saved cash flows do not meet the configured discount-rate target within this scenario. A finite IRR above that target can support the corresponding cash-flow interpretation, but does not validate the conversion yields, prices, equipment quotations or omitted scope.','',
        '**HTL interpretation limit:** repeated launches with identical baseline feed and RF conversion yields selected different adaptive heavy-key cuts, changing fuel output and IRR. The HTL sampled spread therefore includes unresolved solver/cut-selection variability and cannot be attributed purely to incoming composition. This report does not quantify a decomposition of those effects.','',
        '## Baseline cash and utility drivers','',
        'Alternatives save positive operating charges and explicit gate-fee receipts. HTL saves signed input costs in `tea.AOC` and signed output sales in `sys.sales`: negative feedstock prices reduce AOC, while negative disposal-stream prices reduce sales. Comparing these raw totals as if they had identical conventions would be misleading.','',
        'The displayed cash bridge moves negative HTL input costs to receipts and negative output sales to residual-management costs. It is used only when saved feed/product ledgers reconcile to native totals, and preserves `receipts − costs = native sales − native AOC`. Native net utility costs keep their sign. Raw metrics and the bridge’s status remain in `case_metrics.csv`; this does not manufacture a missing native gross-OPEX field.',''])
    rows = []
    for tech in data['techs']:
        b = baseline.get(tech,{})
        rows.append([LABELS[tech],fmt(b.get('feedstock_price_usd_t'),4),
            fmt(b.get('reconciled_cash_cost_usd_year')/1e6 if number(b.get('reconciled_cash_cost_usd_year')) is not None else None),
            fmt(b.get('reconciled_receipts_usd_year')/1e6 if number(b.get('reconciled_receipts_usd_year')) is not None else None),
            fmt(b.get('purchased_heat_mj_t_wet')),fmt(b.get('net_electricity_kwh_t_wet'))])
    lines.extend([markdown_table(['Route','Feed USD/wet t','Cash costs, MUSD/y','Receipts, MUSD/y','Purchased heat, MJ/wet t','Net electricity, kWh/wet t'],rows),'',
        'Feed-price values above are the actual retained shared basis, not a generic tipping-fee default. A negative price denotes the single gate-fee credit; the bridge does not add a second fee. Gross external heating can coexist with internal heat recovery; electricity is net import (+) or export (−). HTL physical purchased heat is reported unavailable if no matching physical inventory was retained.',''])
    tariff = config.get('electricity_price_usd_kwh')
    if number(tariff) is not None:
        lines.extend([f"The explicit common electricity tariff is {fmt(tariff,4)} USD/kWh in the configured constant-dollar basis. The feed/gate fee and tariff are shared scenario inputs, rather than sampled composition parameters.",''])
    cash_rows = [[LABELS[r['technology']],
        'Operating cash costs' if r['metric']=='reconciled_cash_cost_usd_year' else 'Cash receipts',r['count'],
        *[fmt(r.get(k)/1e6) if number(r.get(k)) is not None else 'unavailable' for k in ('p05','median','p95')]]
        for r in data['cash_summary']]
    lines.extend(['### Sampled operating costs and receipts','',
        markdown_table(['Route','Cash measure, million USD/year','Finite n','p05','Median','p95'],cash_rows),'',
        'These additional empirical quantiles use the reconciled saved ledgers, with baseline excluded. Receipts include the single feed/gate-fee credit and configured product/energy sales; operating cash costs exclude financing and capital expenditure. Native net utility charges retain their sign, so this is a harmonized cash presentation rather than a newly inferred all-positive gross engineering OPEX. Missing ledgers remain unavailable.',''])
    for tech in data['techs']:
        cs = [r for r in data['costs'] if r['population']=='baseline' and r['technology']==tech and r['usd_year']>0]
        rs = [r for r in data['revenues'] if r['population']=='baseline' and r['technology']==tech and r['usd_year']>0]
        ctop = sorted(cs,key=lambda r:r['usd_year'],reverse=True)[:3]
        rtop = sorted(rs,key=lambda r:r['usd_year'],reverse=True)[:2]
        if ctop or rtop:
            desc = f"**{LABELS[tech]}:** "
            desc += 'largest retained baseline cash-cost entries are '+', '.join(f"{r['item']} ({r['usd_year']/1e6:.3f} MUSD/year)" for r in ctop)+'. ' if ctop else ''
            desc += 'Leading receipts are '+', '.join(f"{r['item']} ({r['usd_year']/1e6:.3f} MUSD/year)" for r in rtop)+'.' if rtop else 'No positive baseline receipts are retained.'
            lines.extend([desc,''])
    lines.extend(['These are accounting contributions, not causal sensitivity rankings. Capital equipment observations and their source/scope metadata are retained in `equipment_drivers.csv`; installed equipment is distinct from FCI and must not receive a second installation multiplier.','',
        '## Product composition and prices','',
        'Thermochemical product observations retain whole-product water, dry ultimate/ash fractions, HHV and gas purity where available. Oil and char matrices are elemental pseudocomponents; reported gas species do not establish oil molecular composition, hydrochar proximate analysis, contaminant compliance or market qualification. Product mass across these different outlets is not a common fuel-quality functional unit.','',
        'Default oil/gas/biochar receipts are conditional HHV-parity proxies tied to the configured natural-gas price basis. Hydrochar may use a literature **modeled break-even price** converted to an HHV proxy, not an observed market contract. [Lucian and Fiori](https://www.mdpi.com/1996-1073/10/2/211) distinguish a production cost of 157 EUR/t from a break-even price of 200 EUR/t. The saved `price_basis`, source, currency, price year, FX and deflator control the actual scenario. A proxy with `qualified=false` is not buyer acceptance. HTL’s biofuel/binder prices likewise rely on its configured product/upgrading quality assumptions; compost has a separate bulk-price and marketability basis.',''])
    ps = [r for r in data['products'] if r['population']=='baseline' and r.get('stream_role')=='main_product']
    lines.extend([markdown_table(['Route / product','Quantity kg/h','Whole HHV MJ/kg','Price USD/kg','Market basis/status'],
        [[LABELS[r['technology']]+' / '+r['product'],fmt(r.get('mass_kg_h')),fmt(r.get('hhv_mj_kg')),fmt(r.get('price_usd_kg'),4),
          r.get('market_basis','')+'; '+r.get('market_status','')] for r in ps]),'',
        'AD–CHP exports electricity and heat; its digestate market value is a separate input. No fuel-grade, soil-amendment or renewable-policy premium is inferred here. Product details and sampled observations remain in `product_observations.csv`.','',
        '## What the modeled products contain','',
        '**HTL:** the retained native biofuel cut contains phenolic and nitrogen-heterocycle surrogates as well as hydrocarbons; binder includes fatty-acid/amide and nitrogen-heterocycle surrogates. These chemical lumps define its stream representation. Chemical IDs, formulas and component quantities below are captured native model identities, not a measured chromatogram. A boiling cut, higher modeled HHV or an HHV-linked USD/kg price does not establish diesel/SAF qualification, upgrading success, storage stability or buyer acceptance.','',
        '**AD–CHP:** the Buswell-based raw intermediate gas is represented by methane and carbon dioxide. Its methane fraction and gas quantity respond to the incoming biochemical mixture and the COD-weighted biodegradability correction. CHP consumes that gas to deliver electricity/heat. Captured digestate cake and centrate provide Water/DryMatter masses, not a resolved volatile-solids, ash or nitrogen-species allocation between these streams. Whole-digestate retained ammonia/carbon cannot be assigned to cake without an additional partition model. Formula-free DryMatter adapter HHV placeholders are excluded from product-quality statistics.','',
        '**Composting:** PaCTEA tracks residual organic pools and microbial states, carbon, organic nitrogen, ammonium-N and nitrate-N. The reported dry C/N fractions and C:N ratio describe this state-based representation, not a molecular GC–MS analysis. Fixed carbon/volatile matter, maturity, contaminant certification and nutrient availability are not inferred from those states.','',
        '**HTG methane and hydrogen:** gas species are explicit CH4, H2, CO2 and CO, with modeled NH3/H2S partition and water. Each route has its own specified methanation/purification scenario; actual final bulk mole fractions can change with feed composition even when recovery settings stay fixed. Feed sulfur is uncharacterized: zero S/H2S follows from the selected sulfur-free feed-surrogate formulas and is not evidence of sulfur-free waste or compliant product gas. The models do not claim additional light-hydrocarbon species, electrolyte speciation or gas-grid/fuel-cell qualification. Unconverted organics, nitrogenous water and minerals remain residual streams.','',
        '**HTC:** the product is a transformed CHONS hydrochar matrix mixed with mineral ash and finishing moisture. It is neither pure carbon nor unchanged lipid/protein/carbohydrate feed. Its estimated dry elemental composition, ash dilution and whole-product HHV follow the incoming feed and the specified elemental retentions. Proximate fixed carbon and volatile matter require separate measurements.','',
        '**Fast pyrolysis:** raw oil is an oxygen-/nitrogen-containing elemental organic pseudocomponent plus retained water and entrained material; biochar is a distinct transformed matrix with ash. These bulk balances do not identify individual acids, sugars, phenols, amides or hydrocarbons, and do not justify substituting a wood-oil GC–MS profile for mixed wet waste.','',
        '### Baseline native component identities','',
        'The largest five captured components of each main product are shown on a whole-product mass basis. Blank or surrogate formulas remain unresolved; formula strings on matrix pseudocomponents describe the model representation rather than molecular speciation. Full residual and intermediate observations remain in the CSV tables.',''])
    identities = []
    for p in ps:
        selected = [r for r in data['components'] if r['population']=='baseline' and
            r['technology']==p['technology'] and r['product']==p['product']]
        for r in sorted(selected,key=lambda r:r['mass_kg_h'],reverse=True)[:5]:
            identities.append([LABELS[r['technology']]+' / '+r['product'],r['component'],
                r['common_name'] or 'unresolved / model lump',r['formula'] or 'unresolved',
                fmt(r['mass_kg_h']),fmt(r['whole_product_wt']*100) if r['whole_product_wt'] is not None else 'unavailable'])
    lines.extend([markdown_table(['Route / product','Native ID','Native identity','Formula / surrogate','kg/h','Whole-product wt%'],identities),'',
        '## Did incoming composition change the products?','',
        'The table below uses saved sampled product quantities, bulk quality and effective prices. A narrow or constant series is retained explicitly. The sampled inputs vary incoming composition; specified conversion retentions, recovery fractions, operating settings and reference energy-price coefficients remain fixed unless the saved manifest explicitly varies them. Estimated HHV and an HHV-based effective USD/kg price can therefore change while the underlying USD/MJ price stays fixed. HTL additionally has the unresolved adaptive-cut variability described above, so its observed spread is not a pure feed response.',''])
    variability = []
    displayed = {'mass_kg_h':('Quantity, kg/h',1.),'hhv_mj_kg':('Whole HHV, MJ/kg',1.),
        'nitrogen_dry_wt':('Nitrogen, dry wt%',100.),'ash_dry_wt':('Ash, dry wt%',100.),
        'price_usd_kg':('Effective price, USD/kg',1.),'gas_CH4_dry_mol':('CH4, dry mol%',100.),
        'carbon_nitrogen_ratio':('C:N mass ratio',1.)}
    for r in data['product_summary']:
        if r['count'] and r['metric'] in displayed:
            label,scale = displayed[r['metric']]
            variability.append([LABELS[r['technology']]+' / '+r['product'],label,r['count'],
                *[fmt(r[k]*scale,4) for k in ('p05','median','p95')],r['variation_status']])
    lines.extend([markdown_table(['Route / product','Observation','Finite n','p05','Median','p95','Sampled spread'],variability),'',
        'Full min/max/mean/quantiles, unavailable counts, CHONS fractions, gas mole fractions, nitrogen forms and prices are in `product_summary.csv`. Observed component identities and whole-product mass fractions are in `product_components.csv`; absent molecular detail remains absent.','',
        '## Feed → product → finance associations','',
        'Associations join the saved input design, product observations and financial results by the actual same `sample_id`. SciPy `spearmanr` calculates descriptive rank associations for finite paired observations. Undefined constant inputs/outcomes are labeled, not assigned zero correlation. Each coefficient retains its paired count and exact sample identities in `paired_associations.csv`.','',
        'Dry fractions sum to one: raising one necessarily changes others after closure. Product properties, energy use, sales and cash flows are also coupled by the equations; HTL adds unresolved adaptive-cut selection variability. These associations are consequently descriptive paired observations, not isolated causal effects, controlled sensitivity results, global parameter-importance rankings or validated probabilities. Many associations are inspected; no hypothesis-testing significance claim is made.',''])
    for tech in data['techs']:
        available = [r for r in data['associations'] if r['technology']==tech and r['rho'] is not None]
        fp = sorted((r for r in available if r['link']=='feed_to_product'),key=lambda r:abs(r['rho']),reverse=True)
        pf = sorted((r for r in available if r['link']=='product_to_finance' and r['outcome']=='IRR_percent'),key=lambda r:abs(r['rho']),reverse=True)
        sentences = []
        for group in (fp,pf):
            if group:
                r = group[0]
                sentences.append(f"{r['input']} ↔ {r['outcome']}: ρ={r['rho']:.3f}, paired n={r['count']}")
        if sentences:
            lines.extend([f"**{LABELS[tech]}:** selected strongest observed links by absolute association are "+'; '.join(sentences)+'. These linked outputs share the same closed input design.',''])
    lines.extend([
        '## Environmental and scientific coverage','',
        'Existing LCA exports contain explicit missing-factor/exclusion statuses and different native boundaries. Utility-only subtotals are not complete life-cycle GWP. `lca_coverage.csv` counts their retained status labels. This report does not display a climate ranking or combine incompatible avoided-burden factors. See the original `lca_factor_audit.csv`, `lca_factors_required.json` and `comparison_gaps.csv` before interpreting any retained environmental subtotal.',''])
    properties = data.get('property_basis',{})
    if properties:
        reference = properties.get('htl_reference_package_inspection',{})
        thermal = properties.get('thermal_baseline_package',{}).get('feed_surrogates',{})
        property_rows = []
        for component in ('Lipids','Proteins','Carbohydrates'):
            h,t = reference.get(component,{}),thermal.get(component,{})
            property_rows.append([component,h.get('common_name') or 'unavailable',h.get('CAS') or 'unavailable',
                h.get('formula') or 'unavailable',t.get('resolved_name') or 'unavailable',
                t.get('CAS') or 'unavailable',t.get('formula') or 'unavailable'])
        lines.extend(['### Feed identities and property-package limits','',
            markdown_table(['Feed pool','Inspected HTL identity','HTL CAS','HTL formula','Thermal identity','Thermal CAS','Thermal formula'],property_rows),'',
            'The retained `property_basis.json` records fresh matching-runtime package inspection. Formula-only HTL Proteins lookup resolves urethane (CAS 51-79-6), while the new thermal routes explicitly select L-alanine (CAS 56-41-7). Both have C3H7NO2 and the same molecular weight, so matching elemental feed definitions does not imply matching heat capacity, formation enthalpy, combustion or phase-property models. Palmitic acid and glucose identities match in this inspection.','',
            'Use of the same ThermoSTEAM/chemicals framework does not make these property packages identical. The hydrothermal routes also apply a CoolProp real-water correction, and oil/char matrices use their stated empirical property approximations. AD/composting use aggregate COD/state balances. The native HTL parent source and its property lookup were retained for this screening batch. These differences preclude a strict claim of thermodynamic equivalence across routes.',''])
    vessel_rows = []
    for row in data['baseline']:
        if row['technology'] not in ('htg_methane','htg_hydrogen'):
            continue
        for item in row.get('metrics',{}).get('equipment_audit',[]):
            if item.get('unit_class')!='Conversion' or 'Vertical pressure vessel' not in item.get('name',''):
                continue
            design = item.get('design',{})
            limits = item.get('native_correlation',{}).get('native_bounds',{}).get('Vertical vessel weight',[])
            weight,upper = number(design.get('Weight')),number(limits[1]) if len(limits)>1 else None
            route_config = config.get(row['technology'],{}).get('config',{})
            pressure,temperature = number(route_config.get('pressure_pa')),number(route_config.get('temperature_k'))
            vessel_rows.append([LABELS[row['technology']],fmt(weight),fmt(upper),
                'outside upper bound' if weight is not None and upper is not None and weight>upper else 'within recorded weight bound',
                fmt(design.get('Wall thickness')),fmt(pressure/1e6) if pressure is not None else 'unavailable',
                fmt(temperature-273.15) if temperature is not None else 'unavailable',
                fmt(item.get('installed_usd')/1e6) if number(item.get('installed_usd')) is not None else 'unavailable'])
    if vessel_rows:
        lines.extend(['### Baseline reactor cost-correlation limits','',
            markdown_table(['Route','Single vessel weight, lb','Native maximum, lb','Weight check','Wall, in','Pressure, MPa','Temperature, C','Vessel installed total, MUSD'],vessel_rows),'',
            'The baseline hydrogen reactor exceeds the recorded native vertical-vessel weight correlation maximum; the methane reactor is close to that maximum. Costs outside a correlation bound are extrapolations. The native design thickness and cost correlations do not establish pressure-code compliance, metallurgy suitability, creep/corrosion life or a vendor quotation for these thick-wall high-temperature services. Installed vessel totals include the configured parallel vessels and remain distinct from total plant FCI.',''])
    lines.extend([
        'Conversion/property scenarios, catalyst and salt-removal design, exchanger assumptions, high-pressure vessel metallurgy, product quality, residual treatment and site-specific equipment remain consequential engineering gaps. Shared financial assumptions alone do not remove these gaps. A completed numerical batch is not independent scientific validation.','',
        '## Figures and retained analysis tables',''])
    for filename in figures:
        lines.extend([f'![{filename.replace("_"," ")}]({filename})',''])
    lines.extend(['Tables: `coverage.csv`, `metric_coverage.csv`, `lca_coverage.csv`, `metric_summary.csv`, `case_metrics.csv`, `baseline_metrics.csv`, '
                  '`annual_cost_drivers.csv`, `annual_revenue_drivers.csv`, `cash_summary.csv`, `equipment_drivers.csv`, '
                  '`product_observations.csv`, `product_summary.csv`, `product_components.csv`, '
                  '`paired_associations.csv`, and `status_reasons.csv`.', '',
                  'All tables derive from the named run’s saved checkpoint and reporter observations. No probability distribution, missing process cost or unobserved product characteristic is inferred.',''])
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True,help='One saved run inside deploy/results')
    parser.add_argument('--allow-incomplete',action='store_true',help='Explicitly render an incomplete sampled checkpoint')
    args = parser.parse_args(argv)
    data = prepare(args.run_dir,args.allow_incomplete)
    output = data['run_dir']/'analysis'
    output.mkdir(exist_ok=True)
    write_csv(output/'coverage.csv',data['coverage'])
    write_csv(output/'metric_coverage.csv',data['metric_coverage'])
    write_csv(output/'lca_coverage.csv',data['lca_coverage'])
    write_csv(output/'metric_summary.csv',data['summary'])
    write_csv(output/'case_metrics.csv',data['observations'])
    write_csv(output/'baseline_metrics.csv',[r for r in data['observations'] if r['population']=='baseline'])
    write_csv(output/'annual_cost_drivers.csv',data['costs'])
    write_csv(output/'annual_revenue_drivers.csv',data['revenues'])
    write_csv(output/'cash_summary.csv',data['cash_summary'])
    write_csv(output/'equipment_drivers.csv',data['equipment'])
    write_csv(output/'product_observations.csv',data['products'])
    write_csv(output/'product_summary.csv',data['product_summary'])
    write_csv(output/'product_components.csv',data['components'])
    write_csv(output/'paired_associations.csv',data['associations'])
    write_csv(output/'status_reasons.csv',data['reasons'])
    figures = render_figures(data,output)
    report = output/'report.md'
    report.write_text(report_text(data,figures),encoding='utf-8')
    print(report)
    print(f"Saved {len(data['rows'])} sampled cases across {len(data['techs'])} routes; baseline remains separate.")
    return output


if __name__=='__main__':
    main()
