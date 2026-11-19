# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Run composition-based recovery cases without initializing HTL/ML models."""
import argparse
import json
from pathlib import Path
import sys
from dataclasses import replace

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from technologies import (create_ad_system, create_composting_system, ADConfig,
                              CompostConfig, Feedstock, Economics, FEEDSTOCKS,
                              create_fast_pyrolysis_system, create_htc_system, create_htg_system,
                              FastPyrolysisConfig, HTCConfig, HTGConfig)
    from technologies._composition import feedstock_module
else:
    from . import (create_ad_system, create_composting_system, ADConfig,
                   CompostConfig, Feedstock, Economics, FEEDSTOCKS,
                   create_fast_pyrolysis_system, create_htc_system, create_htg_system,
                   FastPyrolysisConfig, HTCConfig, HTGConfig)
    from ._composition import feedstock_module


def main():
    parser = argparse.ArgumentParser(description='Run composition-based waste-treatment technology cases.')
    parser.add_argument('technology', choices=('ad','composting','fast_pyrolysis','htc','htg'))
    capacity = parser.add_mutually_exclusive_group(required=True)
    capacity.add_argument('--capacity-tpd', type=float)
    capacity.add_argument('--biogas-nm3-h', type=float)
    parser.add_argument('--feedstock', choices=feedstock_module().SUPPORTED_FEEDSTOCKS, default='food')
    parser.add_argument('--route', choices=('chp','rng','methane','hydrogen'))
    parser.add_argument('--cepci', type=float, default=596.2)
    parser.add_argument('--factors', type=Path, help='JSON mapping activity keys to kg CO2e/unit')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--config', type=Path, help='Thermochemical configuration JSON')
    parser.add_argument('--markets', type=Path, help='ProductMarket contract or price-proxy mappings keyed by product')
    parser.add_argument('--economics', type=Path, help='Economics overrides; finance defaults follow HTL')
    parser.add_argument('--cooling-price-usd-mj', type=float, default=0.)
    parser.add_argument('--characterization', type=Path, help='JSON overrides for the canonical feedstock characterization')
    parser.add_argument('--bulking-composition', type=Path, help='JSON wet mass fractions for the bulking material')
    parser.add_argument('--ad-biodegradability', type=float, help='Independently determined degradable COD/total COD')
    parser.add_argument('--reactor-conversion', type=float, help='Achieved fraction of corrected BMP at operating conditions')
    parser.add_argument('--product-price-usd-t', type=float, default=0.)
    parser.add_argument('--marketable-fraction', type=float, default=1.)
    parser.add_argument('--empirical', action='store_true', help='Explicitly select the legacy screening conversion')
    args = parser.parse_args()
    economic_overrides = json.loads(args.economics.read_text(encoding='utf-8-sig')) if args.economics else {}
    economic_overrides.setdefault('target_cepci',args.cepci)
    e = Economics.from_htl(**economic_overrides)
    characterization = json.loads(args.characterization.read_text(encoding='utf-8-sig')) if args.characterization else None
    if args.technology in ('fast_pyrolysis','htc','htg'):
        if args.capacity_tpd is None:
            parser.error('Thermochemical technologies require --capacity-tpd')
        if args.empirical or args.product_price_usd_t or args.marketable_fraction != 1:
            parser.error('Thermochemical product pricing requires --markets with an explicit price and quality basis')
        if args.route and args.route not in ('methane','hydrogen'):
            parser.error('Thermochemical --route must be methane or hydrogen for HTG')
        if args.route and args.technology != 'htg':
            parser.error('--route is only applicable to AD or HTG')
        configs = dict(fast_pyrolysis=FastPyrolysisConfig,htc=HTCConfig,htg=HTGConfig)
        factories = dict(fast_pyrolysis=create_fast_pyrolysis_system,htc=create_htc_system,htg=create_htg_system)
        settings = json.loads(args.config.read_text(encoding='utf-8-sig')) if args.config else {}
        if args.technology == 'htg' and args.route:
            if 'route' in settings and settings['route'] != args.route:
                parser.error('--route and configuration route disagree')
            settings['route'] = args.route
        markets = json.loads(args.markets.read_text(encoding='utf-8-sig')) if args.markets else None
        try:
            result = factories[args.technology](args.capacity_tpd,feedstock=args.feedstock,
                config=configs[args.technology](**settings),economics=e,characterization=characterization,
                markets=markets,cooling_price_usd_mj=args.cooling_price_usd_mj)
        except (ValueError,TypeError) as error:
            parser.error(str(error))
        output = json.loads(result.to_json())
        if args.factors:
            output['lca'] = result.lca(json.loads(args.factors.read_text(encoding='utf-8-sig')))
        content = json.dumps(output,indent=2,allow_nan=False)
        if args.output:
            args.output.write_text(content+'\n',encoding='utf-8')
        else:
            print(content)
        return
    if args.route and args.route not in ('chp','rng'):
        parser.error('AD --route must be chp or rng')
    if args.config or args.markets:
        parser.error('--config and --markets apply to thermochemical technologies')
    try:
        f = (FEEDSTOCKS[args.feedstock] if args.empirical else
             Feedstock.from_name(args.feedstock,characterization=characterization))
        f = replace(f,product_price_usd_t=args.product_price_usd_t,
                    marketable_fraction=args.marketable_fraction)
    except (ValueError, KeyError) as error:
        parser.error(str(error))
    if args.technology == 'ad':
        try:
            result = create_ad_system(args.capacity_tpd, biogas_nm3_h=args.biogas_nm3_h,
                feedstock=f, config=ADConfig(route=args.route or 'chp',
                    conversion_model='empirical' if args.empirical else 'cornell',
                    biodegradability=args.ad_biodegradability,reactor_conversion=args.reactor_conversion),economics=e)
        except ValueError as error:
            parser.error(str(error))
    else:
        if args.biogas_nm3_h is not None:
            parser.error('Composting requires --capacity-tpd')
        bulking = json.loads(args.bulking_composition.read_text(encoding='utf-8-sig')) if args.bulking_composition else None
        try:
            result = create_composting_system(args.capacity_tpd,feedstock=f,economics=e,
                config=CompostConfig(conversion_model='empirical' if args.empirical else 'pactea'),
                bulking_composition=bulking)
        except ValueError as error:
            parser.error(str(error))
    output = json.loads(result.to_json())
    if args.factors:
        output['lca'] = result.lca(json.loads(args.factors.read_text(encoding='utf-8')))
    content = json.dumps(output, indent=2, allow_nan=False)
    if args.output:
        args.output.write_text(content+'\n', encoding='utf-8')
    else:
        print(content)


if __name__ == '__main__':
    main()
