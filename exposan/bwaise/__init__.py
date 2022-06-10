#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Yalin Li <mailto.yalin.li@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

# %%

# =============================================================================
# Path management
# =============================================================================

import os
bwaise_path = os.path.dirname(__file__)
data_path = os.path.join(bwaise_path, 'data')
results_path = os.path.join(bwaise_path, 'results')
figures_path = os.path.join(bwaise_path, 'figures')
# To save simulation data and generated figures
if not os.path.isdir(results_path): os.mkdir(results_path)
if not os.path.isdir(figures_path): os.mkdir(figures_path)
del os

# %%

# =============================================================================
# Load components and systems
# =============================================================================

import qsdsan as qs
from . import _cmps, _lca_data, systems, models
from ._process_settings import *
from ._cmps import *
from ._lca_data import *
from .systems import *
from .models import *
from exposan.bwaise._cmps import create_components

currency = qs.currency

# %%

# =============================================================================
# Util functions
# =============================================================================
from qsdsan.utils import time_printer

def update_lca_data(kind):
    '''
    Load impact indicator and impact item data.

    Parameters
    ----------
    kind : str
        "original" loads the data from Trimmer et al.
        (TRACI, ecoinvent v3.2),
        "new" loads the data for ReCiPe and TRACI
        (ecoinvent 3.7.1, at the point of substitution).
    '''
    global lca_data_kind

    if lca_data_kind != kind:
        load_lca_data(kind)
        batch_create_stream_items(kind)

        for lca in (lcaA, lcaB, lcaC):
            for i in lca.lca_streams:
                # To refresh the impact items
                source_ID = i.stream_impact_item.source.ID
                i.stream_impact_item.source = ImpactItem.get_item(source_ID)

        Biogas_CFs = ImpactItem.get_item('Biogas_item').CFs
        for k, v in Biogas_CFs.items():
            Biogas_CFs[k] = v * get_biogas_factor()

        for i in sysA, sysB, sysC:
            i.simulate()

        lca_data_kind = kind


def get_total_inputs(unit, multiplier=1):
    if len(unit.ins) == 0: # Excretion units do not have ins
        ins = unit.outs
    else:
        ins = unit.ins
    inputs = {}
    inputs['COD'] = sum(i.COD*i.F_vol/1e3 for i in ins)
    inputs['N'] = sum(i.TN*i.F_vol/1e3 for i in ins)
    inputs['NH3'] = sum(i.imass['NH3'] for i in ins)
    inputs['P'] = sum(i.TP*i.F_vol/1e3 for i in ins)
    inputs['K'] = sum(i.TK*i.F_vol/1e3 for i in ins)
    for i, j in inputs.items():
        inputs[i] = j * multiplier
    return inputs


def get_recovery(ins, outs, multiplier=1):
    try: iter(outs)
    except: outs = (outs,)

    non_g = tuple(i for i in outs if (i.phase != 'g' and 'los' not in i.ID))
    recovery = {}
    recovery['COD'] = sum(i.COD*i.F_vol/1e3 for i in non_g)
    recovery['N'] = sum(i.TN*i.F_vol/1e3 for i in non_g)
    recovery['NH3'] = sum(i.imass['NH3'] for i in non_g)
    recovery['P'] = sum(i.TP*i.F_vol/1e3 for i in non_g)
    recovery['K'] = sum(i.TK*i.F_vol/1e3 for i in non_g)

    for i, j in recovery.items():
        inputs = get_total_inputs(ins, multiplier)
        recovery[i] /= inputs[i]

    return recovery


sys_dct = {
    'ppl': dict(sysA=get_ppl('exist'), sysB=get_ppl('alt'), sysC=get_ppl('exist')),
    'input_unit': dict(sysA=A1, sysB=B1, sysC=C1),
    'liq_unit': dict(sysA=A13, sysB=B13, sysC=C14),
    'sol_unit': dict(sysA=A12, sysB=B12, sysC=C13),
    'gas_unit': dict(sysA=None, sysB=B15, sysC=None),
    'stream_dct': dict(sysA=streamsA, sysB=streamsB, sysC=streamsC),
    'TEA': dict(sysA=teaA, sysB=teaB, sysC=teaC),
    'LCA': dict(sysA=lcaA, sysB=lcaB, sysC=lcaC),
    'cache': dict(sysA={}, sysB={}, sysC={}),
    }


def cache_recoveries(sys):
    sys_dct['ppl'][sys.ID] = ppl = get_ppl('alt') if sys.ID=='sysB' else get_ppl('exist')
    total_COD = get_total_inputs(sys_dct['input_unit'][sys.ID], ppl)['COD']

    if sys_dct['gas_unit'][sys.ID]:
        gas_mol = sys_dct['gas_unit'][sys.ID].outs[0].imol['CH4']
        gas_COD = gas_mol*1e3*biogas_energy/14e3/total_COD
    else:
        gas_COD = 0

    sys_dct['cache'][sys.ID] = cache = {
        'liq': get_recovery(ins=sys_dct['input_unit'][sys.ID],
                            outs=sys_dct['liq_unit'][sys.ID].ins,
                            multiplier=ppl),
        'sol': get_recovery(ins=sys_dct['input_unit'][sys.ID],
                            outs=sys_dct['sol_unit'][sys.ID].ins,
                            multiplier=ppl),
        'gas': dict(COD=gas_COD, N=0, P=0, K=0)
        }
    return cache


sysA._set_facilities([*sysA.facilities, lambda: cache_recoveries(sysA)])
sysB._set_facilities([*sysB.facilities, lambda: cache_recoveries(sysB)])
sysC._set_facilities([*sysC.facilities, lambda: cache_recoveries(sysC)])


def get_summarizing_functions(system):
    func_dct = {}
    func_dct['get_annual_net_cost'] = lambda tea, ppl: (tea.annualized_equipment_cost-tea.net_earnings)/ppl
    func_dct['get_annual_CAPEX'] = lambda tea, ppl: tea.annualized_equipment_cost/ppl
    func_dct['get_annual_OPEX'] = lambda tea, ppl: tea.AOC/ppl
    func_dct['get_annual_sales'] = lambda tea, ppl: tea.sales/ppl

    for i in ('COD', 'N', 'P', 'K'):
        func_dct[f'get_liq_{i}_recovery'] = \
            lambda sys, i: sys_dct['cache'][sys.ID]['liq'][i]
        func_dct[f'get_sol_{i}_recovery'] = \
            lambda sys, i: sys_dct['cache'][sys.ID]['sol'][i]
        func_dct[f'get_gas_{i}_recovery'] = \
            lambda sys, i: sys_dct['cache'][sys.ID]['gas'][i]
        func_dct[f'get_tot_{i}_recovery'] = \
            lambda sys, i: \
                sys_dct['cache'][sys.ID]['liq'][i] + \
                sys_dct['cache'][sys.ID]['sol'][i] + \
                sys_dct['cache'][sys.ID]['gas'][i]

    return func_dct


def print_summaries(systems):
    if not isinstance(systems, Iterable):
        systems = (systems, )

    for sys in systems:
        func = get_summarizing_functions(sys)
        sys.simulate()
        ppl = sys_dct['ppl'][sys.ID]
        print(f'\n---------- Summary for {sys} ----------\n')
        for i in ('COD', 'N', 'P', 'K'):
            print(f'Total {i} recovery is {func[f"get_tot_{i}_recovery"](sys, i):.1%}, '
                  f'{func[f"get_liq_{i}_recovery"](sys, i):.1%} in liquid, '
                  f'{func[f"get_sol_{i}_recovery"](sys, i):.1%} in solid, '
                  f'{func[f"get_gas_{i}_recovery"](sys, i):.1%} in gas.')
        print('\n')

        tea = sys_dct['TEA'][sys.ID]
        tea.show()

        unit = f'{currency}/cap/yr'
        print(f'\nNet cost: {func["get_annual_net_cost"](tea, ppl):.1f} {unit}.')
        print(f'Capital: {func["get_annual_CAPEX"](tea, ppl):.1f} {unit}.')
        print(f'Operating: {func["get_annual_OPEX"](tea, ppl):.1f} {unit}.')
        print(f'Sales: {func["get_annual_sales"](tea, ppl):.1f} {unit}.')
        print('\n')

        lca = sys_dct['LCA'][sys.ID]
        lca.show()
        print('\n')

        for ind in lca.indicators:
            unit = f'{ind.unit}/cap/yr'
            print(f'\nImpact indicator {ind.ID}:')

            val = lca.total_impacts[ind.ID]/lca.lifetime/ppl
            print(f'\nNet emission: {val:.1f} {unit}.')

            val = lca.total_construction_impacts[ind.ID]/lca.lifetime/ppl
            print(f'Construction: {val:.1f} {unit}.')

            val = lca.total_transportation_impacts[ind.ID]/lca.lifetime/ppl
            print(f'Transportation: {val:.1f} {unit}.')

            val = lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='direct_emission')[ind.ID] \
                /lca.lifetime/ppl
            print(f'Direct emission: {val:.1f} {unit}.')

            val = lca.get_stream_impacts(stream_items=lca.stream_inventory, kind='offset')[ind.ID] \
                /lca.lifetime/ppl
            print(f'Offset: {val:.1f} {unit}.')

            val = lca.total_other_impacts[ind.ID]/lca.lifetime/ppl
            print(f'Other: {val:.1} {unit}.\n')


def save_all_reports():
    import os
    if not os.path.isdir(results_path):
        os.path.mkdir(results_path)
    for i in (sysA, sysB, sysC, lcaA, lcaB, lcaC):
        if isinstance(i, System):
            i.simulate()
            i.save_report(os.path.join(results_path, f'{i.ID}_report.xlsx'))
        else:
            i.save_report(os.path.join(results_path, f'{i.system.ID}_lca.xlsx'))

__all__ = ('sysA', 'sysB', 'sysC', 'teaA', 'teaB', 'teaC', 'lcaA', 'lcaB', 'lcaC',
           'print_summaries', 'save_all_reports', 'update_lca_data',
           *(i.ID for i in sysA.units),
           *(i.ID for i in sysB.units),
           *(i.ID for i in sysC.units),
           )


@time_printer
def evaluate(model, samples=None):
    if samples is not None:
        model.load_samples(samples)
    model.evaluate()

def get_key_metrics(model, alt_names={}):
    key_metrics = [i for i in model.metrics if 'total' in i.name.lower()]
    key_metrics += [i for i in model.metrics if 'net' in i.name.lower()]
    for old, new in alt_names.items():
        for i in key_metrics:
            i.name = i.name.replace(old, new)
    return key_metrics






__all__ = (
    'evaluate',
    'get_key_metrics',
    'bwaise_path',
    'data_path',
    'results_path',
    'figures_path',
	*_cmps.__all__,
    *_lca_data.__all__,
	*systems.__all__,
    *models.__all__,
	)