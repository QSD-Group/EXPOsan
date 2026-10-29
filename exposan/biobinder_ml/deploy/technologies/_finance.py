# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Use the HTL financial policy and BioSTEAM cash-flow engine for alternatives."""
import ast
from pathlib import Path
from types import SimpleNamespace

FINANCIAL_KEYS = ('IRR', 'duration', 'income_tax', 'operating_days', 'depreciation',
                  'construction_schedule', 'startup_months', 'startup_FOCfrac',
                  'startup_VOCfrac', 'startup_salesfrac', 'WC_over_FCI',
                  'finance_interest', 'finance_years', 'finance_fraction',
                  'property_insurance','labor_burden','maintenance')


def _literal_settings(path, variable):
    """Read literal policy values without importing process simulation modules."""
    tree = ast.parse(path.read_text(encoding='utf-8-sig'))
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign) or not any(
                isinstance(t, ast.Name) and t.id == variable for t in node.targets):
            continue
        value = node.value
        pairs = (zip(value.keys, value.values) if isinstance(value, ast.Dict) else
                 ((ast.Constant(k.arg), k.value) for k in value.keywords)
                 if isinstance(value, ast.Call) else ())
        values = {}
        for key, val in pairs:
            if ast.literal_eval(key) not in FINANCIAL_KEYS:
                continue
            try:
                values[ast.literal_eval(key)] = ast.literal_eval(val)
            except (ValueError, TypeError):
                raise ValueError(f'Financial setting in {path} is no longer literal; pass the actual HTL TEA')
        return values
    raise ValueError(f'Cannot locate {variable} in {path}')


def htl_financial_assumptions(tea=None):
    if tea is None:
        root = Path(__file__).resolve().parents[3]
        values = _literal_settings(root/'htl/_tea.py', 'default_kwargs')
        values.update(_literal_settings(root/'saf/_process_settings.py', 'tea_kwargs'))
    else:
        values = {key:getattr(tea,key) for key in FINANCIAL_KEYS}
    duration = values.pop('duration')
    values['start_year'], values['lifetime'] = int(duration[0]), int(duration[1]-duration[0])
    values['discount_rate'] = values.pop('IRR')
    values['working_capital_fraction'] = values.pop('WC_over_FCI')
    values['insurance_fraction'] = values.pop('property_insurance')
    values['maintenance_fraction'] = values.pop('maintenance')
    values['construction_schedule'] = tuple(values['construction_schedule'])
    return values


def cashflow_analysis(e, fixed, land, foc, voc, sales, *, native_system=None, return_tea=False):
    """Run BioSTEAM's construction, startup, depreciation, tax and debt engine.

    Bundled alternative-plant capital is already installed: HTL equipment
    markups are not reapplied. Land is a nondepreciable cash purchase at the
    end of construction and recovered at project end, without debt financing.
    """
    import numpy as np
    import biosteam as bst

    class SnapshotTEA(bst.TEA):
        @property
        def sales(self):
            return sales

        @property
        def VOC(self):
            return voc

        @property
        def TCI(self):
            return fixed*(1+self.WC_over_FCI)+land

        def _DPI(self, installed):
            return fixed

        def _TDC(self, dpi):
            return dpi

        def _FCI(self, tdc):
            return tdc

        def _FOC(self, fci):
            return foc

        def get_cashflow_table(self):
            table = super().get_cashflow_table()
            start = len(e.construction_schedule)-1
            table.iloc[start,table.columns.get_loc('Working capital [MM$]')] += land/1e6
            table.iloc[-1,table.columns.get_loc('Working capital [MM$]')] -= land/1e6
            table.iloc[start,table.columns.get_loc('Cash flow [MM$]')] -= land/1e6
            table.iloc[-1,table.columns.get_loc('Cash flow [MM$]')] += land/1e6
            discounted = table['Cash flow [MM$]'].to_numpy()/(1+self.IRR)**self._get_duration_array()
            table['Net present value (NPV) [MM$]'] = discounted
            table['Cumulative NPV [MM$]'] = np.cumsum(discounted)
            return table

        @property
        def cashflow_array(self):
            return self.get_cashflow_table()['Cash flow [MM$]'].to_numpy()*1e6

        @property
        def NPV(self):
            return float(np.sum(self.cashflow_array/(1+self.IRR)**self._get_duration_array()))

    snapshot = SimpleNamespace(cost_units=(), feeds=(), products=(),
                               sales=sales, material_cost=voc, utility_cost=0.,
                               installed_equipment_cost=fixed)
    tea = SnapshotTEA(native_system if native_system is not None else snapshot, IRR=e.discount_rate,
                      duration=(e.start_year,e.start_year+e.lifetime),
                      depreciation=e.depreciation, income_tax=e.income_tax,
                      operating_days=e.operating_days, lang_factor=None,
                      construction_schedule=e.construction_schedule,
                      startup_months=e.startup_months, startup_FOCfrac=e.startup_FOCfrac,
                      startup_VOCfrac=e.startup_VOCfrac, startup_salesfrac=e.startup_salesfrac,
                      WC_over_FCI=e.working_capital_fraction,
                      finance_interest=e.finance_interest, finance_years=e.finance_years,
                      finance_fraction=e.finance_fraction)
    table = tea.get_cashflow_table()
    flows = table['Cash flow [MM$]'].to_numpy()*1e6
    periods = tea._get_duration_array()
    npv = float(np.sum(flows/(1+e.discount_rate)**periods))
    analysis = dict(cashflows=flows.tolist(), cashflow_years=table.index.tolist(),
                discount_periods=periods.tolist(), NPV=npv,
                annual_tax=float(table['Tax [MM$]'].iloc[-1]*1e6),
                depreciation_usd=(table['Depreciation [MM$]']*1e6).tolist(),
                loan_draw_usd=(table['Loan [MM$]']*1e6).tolist(),
                loan_payment_usd=(table['Loan payment [MM$]']*1e6).tolist(),
                tax_usd=(table['Tax [MM$]']*1e6).tolist(),
                financial_assumptions=dict(duration=list(tea.duration),
                    depreciation=e.depreciation, construction_schedule=list(e.construction_schedule),
                    startup_months=e.startup_months, startup_FOCfrac=e.startup_FOCfrac,
                    startup_VOCfrac=e.startup_VOCfrac, startup_salesfrac=e.startup_salesfrac,
                    finance_interest=e.finance_interest, finance_years=e.finance_years,
                    finance_fraction=e.finance_fraction, income_tax=e.income_tax,
                    discount_rate=e.discount_rate, operating_days=e.operating_days,
                    working_capital_fraction=e.working_capital_fraction,
                    cashflow_engine='BioSTEAM TEA'))
    return (analysis,tea) if return_tea else analysis
