# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Revalue saved HTC cash flows at alternate prices without process simulations."""
import argparse
import csv
import json
from pathlib import Path

import compare_technologies_lhs as study


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-dir',type=Path,required=True)
    args = parser.parse_args(argv)
    run = args.run_dir.resolve()
    if not run.is_relative_to(study.RESULTS.resolve()):
        raise ValueError('Select a saved run inside deploy/results')
    status = json.loads((run/'run_status.json').read_text())
    if status['status']!='COMPLETE':
        raise ValueError('The sampled study must be complete')
    manifest = json.loads((run/'manifest.json').read_text())
    config = manifest['config']
    study.bootstrap(config.get('python_paths',()))
    tools = study.recovery_tools()
    from scipy.optimize import brentq
    import numpy as np
    from importlib import import_module
    finance_module = import_module(tools.__name__+'._finance')
    reference = config['fast_pyrolysis']['markets']['bio_oil']
    if reference['unit']!='USD/MJ-HHV':
        raise ValueError('Energy parity requires a saved USD/MJ-HHV reference')
    # https://www.eia.gov/dnav/ng/hist/rngwhhda.htm : 2020 Henry Hub
    # 2.03 USD/MMBtu = 0.001924068754236 USD/MJ; retain the saved coefficient.
    parity_price = reference['price']*reference['usd_per_currency']*reference['year_factor']
    baseline = json.loads((run/'baseline_results.json').read_text())
    rows = baseline+[json.loads(s) for s in (run/'results.jsonl').read_text().splitlines() if s.strip()]
    observations = []
    for row in rows:
        if row['technology']!='htc' or row['status']!='OK':
            continue
        m = row['metrics']
        f = m['financial_basis']
        keys = ('discount_rate','income_tax','operating_days','depreciation',
            'construction_schedule','startup_months','startup_FOCfrac','startup_VOCfrac',
            'startup_salesfrac','working_capital_fraction','finance_interest','finance_years',
            'finance_fraction')
        e = tools.Economics.from_htl(**{k:f[k] for k in keys},
            start_year=f['duration'][0],lifetime=f['duration'][1]-f['duration'][0])
        fixed = m['FCI_usd']
        land = m['TCI_usd']-fixed*(1+e.working_capital_fraction)
        costs,revenues = m['annual_costs'],m['annual_revenues']
        foc = sum(costs[k] for k in ('maintenance','insurance','labor','land_lease'))
        voc = sum(costs.values())-foc
        sales = sum(revenues.values())
        def revalue(sales_value):
            return finance_module.cashflow_analysis(e,fixed,land,foc,voc,sales_value,return_tea=True)
        reproduced,tea = revalue(sales)
        if not np.allclose(reproduced['cashflows'],m['cashflow_usd'],rtol=1e-9,atol=.01):
            raise ValueError(f'Cash-flow reproduction failed for HTC sample {row["sample_id"]}')
        if not np.isclose(reproduced['NPV'],m['NPV_usd'],rtol=1e-9,atol=.01):
            raise ValueError('Saved NPV did not reconcile')
        quality = m['product_quality']['hydrochar']
        market = m['product_markets']['hydrochar']
        annual_mass = market['marketed_kg_year']
        annual_energy = annual_mass*quality['hhv_mj_kg']
        parity_revenue = annual_energy*parity_price
        base_product = revenues['hydrochar']
        fee = revenues['gate_fee']
        scenarios = [('saved_prices_and_gate_fee',sales,base_product,fee),
            ('energy_parity_and_gate_fee',sales-base_product+parity_revenue,parity_revenue,fee),
            ('saved_prices_without_gate_fee',sales-fee,base_product,0.),
            ('energy_parity_without_gate_fee',sales-base_product+parity_revenue-fee,parity_revenue,0.)]
        other_receipts = sales-base_product
        def npv_at_product_price(price):
            return revalue(other_receipts+annual_mass*price)[0]['NPV']
        if npv_at_product_price(0)>=0:
            threshold = 0.
        else:
            upper = max(base_product/annual_mass,1.)
            while npv_at_product_price(upper)<0:
                upper *= 2
                if upper>1e5:
                    raise ValueError('No bounded product price threshold')
            threshold = brentq(npv_at_product_price,0.,upper,xtol=1e-10)
        for name,scenario_sales,product_receipts,scenario_fee in scenarios:
            analysis,tea = revalue(scenario_sales)
            irr,reason = study.native_irr(tea)
            observations.append(dict(population='baseline' if row['sample_id']<0 else 'sampled',
                sample_id=row['sample_id'],technology='htc',scenario=name,
                FCI_usd=fixed,annual_cash_costs_usd=sum(costs.values()),
                annual_receipts_usd=scenario_sales,annual_product_receipts_usd=product_receipts,
                annual_gate_fee_usd=scenario_fee,
                effective_product_price_usd_kg=product_receipts/annual_mass,
                effective_product_price_usd_gj=product_receipts/annual_energy*1000,
                NPV_usd=analysis['NPV'],IRR_percent=irr,IRR_status=reason,
                product_price_for_zero_NPV_usd_t=threshold*1000,
                product_price_for_zero_NPV_usd_gj=threshold/quality['hhv_mj_kg']*1000))
    output = run/'analysis'
    output.mkdir(exist_ok=True)
    with (output/'htc_price_scenarios.csv').open('w',encoding='utf-8-sig',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(observations[0]))
        writer.writeheader();writer.writerows(observations)
    lines = ['# HTC price and gate-fee scenarios','',
        'These are cash-flow counterfactuals on saved HTC plant snapshots. Feed, conversion, utility quantities, operating costs, capital, tax, debt and startup policy stay fixed. Each original annual cash-flow vector and NPV was reproduced before revaluation. No additional process runs or LHS samples were generated.','',
        'Energy parity applies the same saved 2020 natural-gas USD/MJ-HHV coefficient used for fast-pyrolysis products and HTG gases. The original hydrochar proxy remains the configured literature modeled break-even proxy, with its saved EUR/USD conversion and price-year adjustment. Neither scenario establishes a qualified outlet.','',
        '| Scenario | Baseline price USD/t | Baseline NPV MUSD | Baseline IRR % | Sampled IRR p05 / median / p95 | Finite IRR n |',
        '| --- | --- | --- | --- | --- | --- |']
    for scenario in dict.fromkeys(r['scenario'] for r in observations):
        b = next(r for r in observations if r['population']=='baseline' and r['scenario']==scenario)
        selected = [r['IRR_percent'] for r in observations if r['population']=='sampled' and r['scenario']==scenario and r['IRR_percent'] is not None]
        quantiles = ' / '.join(f'{v:.2f}' for v in np.quantile(selected,[.05,.5,.95])) if selected else 'unavailable'
        baseline_irr = f'{b["IRR_percent"]:.2f}' if b['IRR_percent'] is not None else 'unavailable'
        lines.append(f'| {scenario} | {b["effective_product_price_usd_kg"]*1000:.2f} | {b["NPV_usd"]/1e6:.2f} | {baseline_irr} | {quantiles} | {len(selected)} |')
    b = next(r for r in observations if r['population']=='baseline')
    lines.extend(['',f'At the saved gate fee and all fixed baseline plant costs, the HTC hydrochar price required for NPV=0 at the configured discount rate is **{b["product_price_for_zero_NPV_usd_t"]:.2f} USD/t** ({b["product_price_for_zero_NPV_usd_gj"]:.2f} USD/GJ-HHV). This is a model financial threshold, not a market-price estimate.','',
        'Removing gate fees is an explicit alternate contract scenario; no fee is added to negative feed prices. The common fee is identical across routes in the primary study, so it contributes materially to viability without explaining route-to-route differences on its own. A change in modeled sale price here does not change reactor size, drying demand or product composition.',''])
    (output/'htc_price_scenarios.md').write_text('\n'.join(lines),encoding='utf-8')
    print(output/'htc_price_scenarios.md')


if __name__=='__main__':
    main()
