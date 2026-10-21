# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Paired California food/sludge incentive comparisons.

Uses the eight process cases and all native reporting from compare_six_systems.
Each process is simulated once, then its TEA is evaluated without/with support.
Edit the scenario assumptions below before running. Output workbooks and logs
use the existing comparison_results/<timestamp> directory convention.
The default is 16 sludge scenarios at the CDF median. Use --both-feedstocks
for 32 scenarios: food at the CA food-waste gate fee and sludge at the CA landfill fee.
"""

import ast
from contextlib import contextmanager
from pathlib import Path
import traceback

import numpy as np
import pandas as pd
from biosteam._tea import taxable_earnings_with_fowarded_losses
import compare_six_systems as comparison
from comparison_reporting import collect, export, failed_snapshot
from exposan.biobinder_ml.deploy._incentives import (
    annual_quantity, capital_grant, clean_fuel_support, section45z_support,
)


START_YEAR = 2027
LIFETIME = 30
STATE = 'CA'
SLUDGE_CDF_QUANTILE = 0.5
OUTPUT_DIR = Path(__file__).resolve().parent / 'comparison_results'

# https://www.irs.gov/instructions/i7218
# https://www.irs.gov/pub/irs-drop/n-26-53.pdf
# 45Z: $1/physical US gal alternative amount before annual indexing, with
# prevailing-wage/apprenticeship compliance assumed; eligible sales end 2029.
# Dollars are converted to the model's constant 2020-USD basis below.
CREDIT_45Z_USD_GAL = 1.0
CREDIT_45Z_END_YEAR = 2029
CREDIT_45Z_REALIZATION = 0.95
MJ_PER_MMBTU = 1055.05585262
CI_45Z_G_CO2E_MJ = 40.0
CI_LCFS_G_CO2E_MJ = 40.0
# 40 g CO2e/MJ is an illustrative qualified-SAF pathway assumption, not an
# estimated or certified CI. 45Z requires the permitted 45ZCF-GREET/SAF
# methodology or a provisional emissions rate; LCFS requires CA-GREET and
# CARB pathway certification. These separate inputs may differ. Native QSDsan
# results require boundary/CF revision before either regulatory application:
# replace generic landfill displacement with the permitted sludge counterfactual,
# use the required feedstock transport, H2, CA electricity, coproduct allocation,
# fuel distribution and use-stage accounting. A 'GREET' indicator label alone
# does not implement either regulatory method. Native LCA reporting is preserved.
# Diagnostic CI below is annual total native GWP (kg CO2e)*1000 divided by
# annual LHV of exported combustible liquid products (MJ); it includes native
# avoided-disposal credits and is NOT used as either policy CI.

# https://ww2.arb.ca.gov/sites/default/files/2025-07/2025_lcfs_fro_oal-approved_unofficial_07012025.pdf
# Table 3: jet benchmarks, g CO2e/MJ; published schedule without additional
# Automatic Acceleration Mechanism triggers. 2045 value continues thereafter.
JET_BENCHMARK = dict(zip(range(2027, 2046), (
    78.63, 77.10, 75.57, 74.03, 69.27, 64.51, 59.75, 54.99, 50.23,
    45.47, 40.71, 35.95, 31.19, 26.44, 23.26, 20.09, 16.92, 13.74, 10.57)))
# https://ww2.arb.ca.gov/our-work/programs/low-carbon-fuel-standard
# Illustrative constant credit price, not a market quotation. All eligible jet
# is assumed sold for California aviation use; project retains the full credit.
LCFS_USD_TONNE = 75.0
LCFS_RETAINED_SHARE = 1.0
QUALIFIED_STREAMS = ('jet',)
# Finished SAF only: no credit on binder, crude, marine intermediate, refinery
# sales, or an untreated 'biofuel' cut. No automatic D-code/RIN is assumed.
# https://www.epa.gov/renewable-fuel-standard-program/approved-pathways-renewable-fuel

# https://www.energy.gov/eere/bioenergy/funding-opportunities
# Approved demonstration capital-award scenario: 20% of eligible upstream
# installed equipment, capped at $5M/project in 2020 USD. These are scenario
# assumptions, not a universal DOE award rate/cap. No per-pilot multiplication.
# Taxable grant treatment: gross asset/depreciation basis remains unchanged;
# receipt enters taxable cash flow in each construction expenditure year.
CAPITAL_GRANT_FRACTION = 0.20
CAPITAL_GRANT_CAP_USD = 5_000_000.0
GRANT_ELIGIBLE_UNITS = ('FeedstockCond', 'HTL')

# https://www.rd.usda.gov/programs-services/energy-programs/biorefinery-renewable-chemical-and-biobased-product-manufacturing-assistance-program
# Approved eligible facility/renewable-material financing scenario: negotiated
# 3% interest, 10-year amortization; not a USDA entitlement or guaranteed rate.
# Native finance fraction is retained on capital net of the grant. Guarantee
# value is never booked as a grant or sale. No CWSRF loan is assigned to the
# entire private upgrading refinery; utility financing requires an asset split.
SUPPORTED_LOAN_INTEREST = 0.03
SUPPORTED_LOAN_YEARS = 10


def sludge_cdf_baseline():
    """Read the existing empirical CDF without importing the UA simulation."""
    source = Path(__file__).resolve().parents[2] / 'analyses' / 'SAF_UA.py'
    tree = ast.parse(source.read_text(encoding='utf-8-sig'))
    arrays = {}
    for node in tree.body:
        if isinstance(node, ast.Assign) and isinstance(node.value, ast.Call):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in (
                        'SLUDGE_CDF_P', 'SLUDGE_CDF_X'):
                    arrays[target.id] = np.asarray(
                        ast.literal_eval(node.value.args[0]), dtype=float)
    p, x = arrays['SLUDGE_CDF_P'], arrays['SLUDGE_CDF_X']
    if len(p) != len(x) or not np.all(np.diff(p) > 0):
        raise ValueError('Invalid sludge CDF probability grid')
    # CDF values are 2020 USD/wet metric tonne. Preserve the module's empirical
    # interpolation (including its local nonmonotonic price segment), without
    # applying Badgett_factor a second time. p=0.5 gives -36.733333 USD/t wet.
    return float(np.interp(SLUDGE_CDF_QUANTILE, p, x))


def policy_dollar_factor():
    """Convert 2026 nominal policy dollars to native constant 2020 dollars."""
    from qsdsan.utils import indices
    pce = indices.tea_indices['PCEPI']
    # https://www.bea.gov/data/personal-consumption-expenditures-price-index
    # If a 2026 annual PCEPI is not present, use the latest index <=2026 as
    # an explicit constant-dollar proxy; record the actual year in the workbook.
    available = [int(y) for y in pce if str(y).isdigit() and int(y) <= 2026]
    year = max(available)
    return float(pce[2020] / pce[year]), year


def make_policy_tea(base):
    class CaliforniaIncentiveTEA(base):
        policy_enabled = False
        policy_rows = None

        def _fill_tax_and_incentives(self, incentives, taxable_cashflow,
                                     nontaxable_cashflow, tax, depreciation):
            super()._fill_tax_and_incentives(
                incentives, taxable_cashflow, nontaxable_cashflow, tax,
                depreciation)
            if not self.policy_enabled:
                return
            taxable_support = np.zeros_like(incentives)
            credit = np.zeros_like(incentives)
            first = self.duration[0] - len(self.construction_schedule)
            for row in self.policy_rows:
                index = row['Year'] - first
                if 0 <= index < len(incentives):
                    taxable_support[index] += row['Taxable_receipts_USD']
                    credit[index] += row['Tax_credit_receipts_USD']
            # LCFS sales and taxable grant receipts enter earnings BEFORE loss
            # carryforward; 45Z transfer proceeds are realized tax-credit cash.
            raw, _, _ = self._taxable_nontaxable_depreciation_cashflows()
            original_tax = self.income_tax * taxable_earnings_with_fowarded_losses(raw)
            supported_tax = self.income_tax * taxable_earnings_with_fowarded_losses(
                raw + taxable_support)
            tax[:] += supported_tax - original_tax
            incentives[:] += taxable_support + credit

        def get_cashflow_table(self):
            table = super().get_cashflow_table()
            if not self.policy_enabled:
                return table
            lcfs = np.zeros(len(table))
            grant = lcfs.copy()
            first = int(table.index[0])
            for row in self.policy_rows:
                index = row['Year'] - first
                if 0 <= index < len(table):
                    if row['Program'] == 'CA LCFS':
                        lcfs[index] += row['Taxable_receipts_USD']
                    elif row['Program'] == 'Capital award':
                        grant[index] += row['Taxable_receipts_USD']
            table['Sales [MM$]'] += lcfs / 1e6
            table['Incentives [MM$]'] -= (lcfs + grant) / 1e6
            table['Taxable capital grant [MM$]'] = grant / 1e6
            table['LCFS sales included above [MM$]'] = lcfs / 1e6
            raw, _, _ = self._taxable_nontaxable_depreciation_cashflows()
            raw = raw + lcfs + grant
            earnings = taxable_earnings_with_fowarded_losses(raw)
            table['Taxed earnings [MM$]'] = earnings / 1e6
            losses = np.zeros(len(table))
            losses[1:] = (raw-earnings).cumsum()[:-1]
            table['Forwarded losses [MM$]'] = losses / 1e6
            return table

    return CaliforniaIncentiveTEA


def build_policy(system, dollar_factor):
    """Return separate dated funding/revenue entries and project inputs."""
    tea = system.TEA
    first, last = tea.duration
    eligible_cost = sum(float(getattr(system.flowsheet.unit, name).installed_cost)
                        for name in GRANT_ELIGIBLE_UNITS)
    grant = capital_grant(eligible_cost, CAPITAL_GRANT_FRACTION,
                          eligible=True, cap=CAPITAL_GRANT_CAP_USD)
    rows = []

    def add(year, program, kind, quantity, rate, taxable=0., credit=0., stream=''):
        rows.append(dict(Year=int(year), Program=program, Treatment=kind,
                         Recipient='HTL project', Stream=stream,
                         Eligible_quantity=quantity,
                         Quantity_unit=('USD eligible capital' if program == 'Capital award'
                                        else 'US gal' if program == '45Z' else 't CO2e'),
                         Rate=rate,
                         Rate_unit=('USD/USD' if program == 'Capital award'
                                    else 'USD/US gal' if program == '45Z' else 'USD/t CO2e'),
                         Taxable_receipts_USD=float(taxable),
                         Tax_credit_receipts_USD=float(credit)))

    for offset, share in enumerate(tea.construction_schedule):
        add(first-len(tea.construction_schedule)+offset, 'Capital award',
            'Taxable capital grant', eligible_cost * share,
            CAPITAL_GRANT_FRACTION, taxable=grant * share)

    for ID in QUALIFIED_STREAMS:
        stream = getattr(system.flowsheet.stream, ID, None)
        if stream is None or stream not in system.products or stream.F_mass <= 0:
            continue
        gallons = annual_quantity(stream, tea, basis='US_gal')
        energy = float(stream.LHV) / 1000 * tea.operating_hours
        if not np.isfinite(energy) or energy <= 0 or gallons <= 0:
            raise ValueError(f'{ID} requires positive physical volume and LHV')
        # https://www.irs.gov/instructions/i7218
        # kg CO2e/mmBTU = g CO2e/MJ * 1055.05585262 / 1000.
        # Round emissions factor to nearest 0.1 (half-up), not downward.
        rate_kg_mmbtu = max(0., CI_45Z_G_CO2E_MJ) * MJ_PER_MMBTU / 1000
        for year in range(first, last):
            # Match native TEA first-year sales proration; subsequent years
            # represent the same annual steady-state operating-hour basis.
            startup = (1.-tea._startup_time*(1.-tea.startup_salesfrac)
                       if year == first else 1.)
            if 2026 <= year <= CREDIT_45Z_END_YEAR:
                rate = section45z_support(
                    1., rate_kg_mmbtu, CREDIT_45Z_USD_GAL * dollar_factor,
                    baseline=50., eligible=True, realization=CREDIT_45Z_REALIZATION)
                add(year, '45Z', 'Realized tax credit', gallons*startup,
                    rate, credit=gallons*startup*rate, stream=ID)
            if year >= min(JET_BENCHMARK):
                benchmark = JET_BENCHMARK[min(year, max(JET_BENCHMARK))]
                tonnes = clean_fuel_support(energy*startup, benchmark,
                                           CI_LCFS_G_CO2E_MJ, 1., eligible=True)
                rate = LCFS_USD_TONNE * LCFS_RETAINED_SHARE * dollar_factor
                add(year, 'CA LCFS', 'Taxable tradable-credit revenue',
                    tonnes, rate, taxable=tonnes*rate, stream=ID)
    return rows, dict(Eligible_capital_USD=eligible_cost,
                      Capital_grant_USD=grant,
                      Qualified_streams=','.join(sorted({r['Stream'] for r in rows if r['Stream']})))


def diagnostic_ci(system):
    liquids = [s for s in system.products
               if s.ID in ('gasoline', 'jet', 'diesel', 'biofuel', 'intermediate', 'biobinder')]
    energy = sum(max(0., float(s.LHV)) for s in liquids) / 1000 * system.operating_hours
    if energy <= 0:
        return None
    return float(system.LCA.get_total_impacts(annual=True)['GWP']) * 1000 / energy


@contextmanager
def comparison_configuration(paired, ledger, factor, index_year,
                             feedstock_price_wet_tonne=None,
                             price_assumptions=None, run_label=None, feedstock='sludge'):
    """Temporarily extend the native runner; restore every binding afterwards."""
    saved = {key: getattr(comparison, key) for key in (
        'FEEDSTOCK', 'FEEDSTOCK_PRICE_WET_TONNE', 'COMMON_TEA_OVERRIDES',
        'INCLUDE_SECOND_REFINERY_ROUTE', 'OUTPUT_DIR', 'collect', 'export')}
    modules = (comparison.joint, comparison.saf, comparison.binder)
    classes = [(module, module.BiobinderTEA) for module in modules]
    from exposan.biobinder_ml import feedprice_dct
    price_assumptions = dict(price_assumptions or {})
    if feedstock == 'food' and feedstock_price_wet_tonne is None:
        feedstock_price_wet_tonne, reference_assumptions = ca_food_fee_basis()
        price_assumptions = dict(reference_assumptions, **price_assumptions)
    feed_price = ((sludge_cdf_baseline() if feedstock == 'sludge'
                   else float(feedprice_dct[feedstock])*1000)
                  if feedstock_price_wet_tonne is None else float(feedstock_price_wet_tonne))
    uses_cdf = feedstock == 'sludge' and feedstock_price_wet_tonne is None
    if not np.isfinite(feed_price):
        raise ValueError('Feedstock price must be finite USD/wet metric tonne')

    def collect_pair(system, case, settings, simulation_error='', console_log=''):
        tea = system.TEA
        common = dict(settings, State=STATE,
                      Policy_CI_45Z_g_MJ=CI_45Z_G_CO2E_MJ,
                      Policy_CI_LCFS_g_MJ=CI_LCFS_G_CO2E_MJ,
                      Sludge_CDF_quantile=(SLUDGE_CDF_QUANTILE
                                          if uses_cdf else None),
                      Distillation_policy='min_purchase_cost',
                      Distillation_helper='exposan.biobinder_ml._distill_utils_flex',
                      Policy_USD_to_2020_USD=factor, Policy_PCEPI_proxy_year=index_year)
        common.update(price_assumptions)
        snapshots = []
        previous_finance = (tea.finance_interest, tea.finance_years, tea.finance_fraction)
        try:
            for enabled in (False, True):
                name = case + ('_with' if enabled else '_without')
                local = dict(common, Incentives=enabled)
                error = simulation_error
                try:
                    tea.policy_enabled = enabled and not error
                    tea.policy_rows = []
                    if tea.policy_enabled:
                        rows, details = build_policy(system, factor)
                        tea.policy_rows = rows
                        local.update(details)
                        tea.finance_interest = SUPPORTED_LOAN_INTEREST
                        tea.finance_years = SUPPORTED_LOAN_YEARS
                        # Reduce debt/equity need by the actual capital award;
                        # keep gross FCI and taxable-grant depreciation basis.
                        tea.finance_fraction = previous_finance[2] * (
                            tea.FCI-details['Capital_grant_USD']) / tea.FCI
                        for row in rows:
                            ledger.append(dict(Case=name, **row))
                    local.update(Finance_interest=tea.finance_interest,
                                 Finance_years=tea.finance_years,
                                 Finance_fraction=tea.finance_fraction)
                    snap = collect(system, name, local, simulation_error=error,
                                   console_log=console_log)
                    snap['metrics'].update(local)
                    if not error:
                        snap['metrics']['Native_diagnostic_CI_g_MJ'] = diagnostic_ci(system)
                        receipts = sum(r['Taxable_receipts_USD'] + r['Tax_credit_receipts_USD']
                                       for r in tea.policy_rows)
                        snap['metrics']['Lifetime_support_receipts_USD'] = receipts
                        for program, label in (('45Z', 'Tax_credit'), ('CA LCFS', 'LCFS_revenue'),
                                               ('Capital award', 'Capital_grant')):
                            snap['metrics']['Lifetime_'+label+'_USD'] = sum(
                                r['Taxable_receipts_USD'] + r['Tax_credit_receipts_USD']
                                for r in tea.policy_rows if r['Program'] == program)
                        years = np.arange(tea.duration[0]-len(tea.construction_schedule), tea.duration[1])
                        cash = np.asarray(tea.cashflow_array).copy()
                        discounted = cash / (1.+tea.IRR)**tea._get_duration_array()
                        if not np.isclose(discounted.sum(), tea.NPV, atol=1e-4):
                            raise RuntimeError('Native cashflow/NPV reconciliation failed')
                        # Native breakdowns remain intact; this block preserves
                        # the exact cash array used by NPV/IRR, including land.
                        snap['blocks'].append(('NPV CASHFLOW RECONCILIATION', pd.DataFrame(
                            dict(Year=years, Cash_flow_USD=cash, Discounted_USD=discounted))))
                        snap['blocks'].append(('INCENTIVE RECEIPTS', pd.DataFrame(tea.policy_rows)))
                except Exception:
                    snap = failed_snapshot(name, local, traceback.format_exc(), console_log)
                snapshots.append(snap)
                paired.append(snap)
        finally:
            tea.policy_enabled = False
            tea.policy_rows = None
            tea.finance_interest, tea.finance_years, tea.finance_fraction = previous_finance
        return snapshots[0]

    def export_pairs(original, path, notes=None):
        # Include both variants even if the native runner skipped a build.
        have = {s['case'] for s in paired}
        for snap in original:
            if snap['case'].endswith('_without'):
                continue
            for suffix in ('_without', '_with'):
                name = snap['case'] + suffix
                if name not in have:
                    paired.append(failed_snapshot(name, {}, snap['metrics'].get(
                        'Error', 'Process build failed'), ''))
        path = Path(path).with_name(Path(path).name.replace(
            'pathway_comparison_', f'pathway_incentives_CA_{feedstock}_'))
        if run_label:
            path = path.with_name(path.stem + '_' + run_label + path.suffix)
        export(paired, path, notes)
        with pd.ExcelWriter(path, engine='openpyxl', mode='a') as writer:
            pd.DataFrame(ledger).to_excel(writer, sheet_name='Incentive_ledger', index=False)
            pd.DataFrame([
                ('State', STATE), ('Start_year', START_YEAR), ('Lifetime_years', LIFETIME),
                ('Feedstock', feedstock),
                ('Sludge_CDF_p', SLUDGE_CDF_QUANTILE if uses_cdf else None),
                ('Sludge_price_2020_USD_wet_tonne', feed_price if feedstock == 'sludge' else None),
                ('Feedstock_price_2020_USD_wet_tonne', feed_price),
                ('CI_45Z_g_CO2e_MJ', CI_45Z_G_CO2E_MJ), ('CI_LCFS_g_CO2e_MJ', CI_LCFS_G_CO2E_MJ),
                ('45Z_USD_gal_before_factor', CREDIT_45Z_USD_GAL),
                ('45Z_realization', CREDIT_45Z_REALIZATION), ('45Z_last_year', CREDIT_45Z_END_YEAR),
                ('LCFS_USD_tonne', LCFS_USD_TONNE), ('LCFS_retained_share', LCFS_RETAINED_SHARE),
                ('Capital_grant_fraction', CAPITAL_GRANT_FRACTION), ('Capital_grant_cap_2020_USD', CAPITAL_GRANT_CAP_USD),
                ('Loan_interest', SUPPORTED_LOAN_INTEREST), ('Loan_years', SUPPORTED_LOAN_YEARS),
                ('USD_to_2020_USD', factor), ('PCEPI_proxy_year', index_year),
            ] + list(price_assumptions.items()), columns=['Assumption', 'Value']).to_excel(
                writer, sheet_name='Policy_inputs', index=False)
            from copy import copy
            for title in ('Incentive_ledger', 'Policy_inputs'):
                sheet = writer.book[title]
                sheet.freeze_panes = 'A2'
                sheet.sheet_view.showGridLines = False
                sheet.auto_filter.ref = sheet.dimensions
                for cell in sheet[1]:
                    cell.font = copy(writer.book['Comparison']['A1'].font)
                    cell.fill = copy(writer.book['Comparison']['A1'].fill)
                for column in sheet.columns:
                    sheet.column_dimensions[column[0].column_letter].width = min(
                        48, max(16, max(len(str(c.value or '')) for c in column[:50])+2))
                for row in sheet.iter_rows(min_row=2):
                    for cell in row:
                        if isinstance(cell.value, float):
                            cell.number_format = '#,##0.00;[Red](#,##0.00)'
        return path

    try:
        comparison.FEEDSTOCK = feedstock
        comparison.FEEDSTOCK_PRICE_WET_TONNE = feed_price
        comparison.COMMON_TEA_OVERRIDES = dict(saved['COMMON_TEA_OVERRIDES'],
                                             duration=(START_YEAR, START_YEAR+LIFETIME))
        comparison.INCLUDE_SECOND_REFINERY_ROUTE = True
        comparison.OUTPUT_DIR = OUTPUT_DIR
        comparison.collect = collect_pair
        comparison.export = export_pairs
        for module, base in classes:
            module.BiobinderTEA = make_policy_tea(base)
        yield
    finally:
        for key, value in saved.items():
            setattr(comparison, key, value)
        for module, base in classes:
            module.BiobinderTEA = base


def run_comparison(*, feedstock_price_wet_tonne=None, price_assumptions=None,
                   run_label=None, feedstock='sludge'):
    """Run eight paired cases; optional feed price is 2020 USD/wet metric tonne."""
    paired, ledger = [], []
    factor, index_year = policy_dollar_factor()
    with comparison_configuration(paired, ledger, factor, index_year,
                                  feedstock_price_wet_tonne,
                                  price_assumptions, run_label, feedstock):
        _, workbook = comparison.run_comparison()
    if len(paired) != 16:
        raise RuntimeError(f'Expected 16 scenarios, got {len(paired)}')
    return paired, workbook


def ca_food_fee_basis():
    """California regional food-waste gate fee in constant 2020 USD/t wet."""
    # https://www.lacsd.org/services/solid-waste/tipping-fees-for-solid-waste-and-recyclables
    # Puente Hills MRF, Whittier: Food Waste $103.70/ton, effective 2026-07-15.
    # Use the food-only rate, not the $130.89/ton food/green-waste program.
    # Interpret the US operator's ton as 2,000 lb of as-received material:
    # 1 short ton = 0.90718474 metric tonne; $103.70/ton = $114.30968/t wet.
    # This is a regional gate-fee benchmark, not a statewide average or an
    # HTL contract. Assume full gate-fee capture, without adding hauling fees,
    # franchise-fee credits or a second avoided-disposal revenue stream.
    # The 2026 annual PCEPI is unavailable in the local index series; use its
    # latest available year <=2026, recorded explicitly as a deflator proxy.
    reported_fee = 103.70
    factor, index_year = policy_dollar_factor()
    metric_fee = reported_fee / 0.90718474
    price = -metric_fee * factor
    assumptions = dict(
        Fee_scenario='CA Puente Hills MRF food-only gate fee',
        Fee_reference_URL='https://www.lacsd.org/services/solid-waste/tipping-fees-for-solid-waste-and-recyclables',
        Fee_effective_date='2026-07-15', Fee_nominal_year=2026,
        Fee_reported_USD_wet_short_ton=reported_fee,
        Fee_ton_interpretation='US short ton, as received',
        Fee_reported_USD_wet_metric_tonne=metric_fee,
        Fee_includes_hauling_and_tipping=False, Fee_captured_share=1.0,
        Fee_USD_to_2020_USD=factor, Fee_PCEPI_proxy_year=index_year,
        Tipping_receipt_2020_USD_wet_metric_tonne=-price)
    return price, assumptions


def run_ca_food_fee_test():
    """Run all eight paired food cases using the California gate-fee benchmark."""
    price, assumptions = ca_food_fee_basis()
    print(f'CA food gate fee: $103.70/wet short ton; '
          f'model feed price {price:.6f} 2020 USD/wet metric tonne')
    return run_comparison(feedstock='food', feedstock_price_wet_tonne=price,
                          price_assumptions=assumptions, run_label='CA_food_gate_fee')


def run_ca_landfill_fee_test():
    """Run the paired comparison using a California avoided-disposal fee."""
    from qsdsan.utils.indices import tea_indices
    # Observed landfill management cost: $130/wet US short ton, including
    # hauling and tipping; 2018-2023 survey average. The $95/WT beneficial-use
    # average and projected $150-260/WT landfill costs are separate scenarios.
    # A short ton is 907.18474 kg; $130/WT = $143.3004707/wet metric tonne.
    # Treat 2023 as a dollar-year proxy for the mixed survey period; convert to
    # constant 2020 USD using PCEPI. Assume the full avoided cost is captured
    # as an HTL service fee, not that the report establishes a gate-only fee.
    # Replace the existing negative feed price; do not add another disposal
    # credit. Native feedstock/crude transportation accounting is retained.
    reported_fee = 130.0
    source_year = 2023
    factor = tea_indices['PCEPI'][2020] / tea_indices['PCEPI'][source_year]
    metric_fee = reported_fee / 0.90718474
    price = -metric_fee * factor
    assumptions = dict(
        Fee_scenario='CA landfill avoided management cost',
        Fee_reported_USD_wet_short_ton=reported_fee,
        Fee_report_period='2018-2023', Fee_dollar_year_proxy=source_year,
        Fee_includes_hauling_and_tipping=True, Fee_captured_share=1.0,
        Fee_reported_USD_wet_metric_tonne=metric_fee,
        Fee_USD_to_2020_USD=factor,
        Tipping_receipt_2020_USD_wet_metric_tonne=-price)
    print(f'CA fee test: ${reported_fee:.2f}/wet short ton = '
          f'${metric_fee:.2f}/wet metric tonne; '
          f'model feed price {price:.6f} 2020 USD/wet metric tonne')
    return run_comparison(feedstock_price_wet_tonne=price,
                          price_assumptions=assumptions,
                          run_label='CA_landfill_fee_test')


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ca-landfill-fee-test', action='store_true',
                        help='Use $130/wet short ton avoided landfill-management fee')
    parser.add_argument('--ca-food-fee-test', action='store_true',
                        help='Run food at the Puente Hills MRF $103.70/short ton gate fee')
    parser.add_argument('--both-feedstocks', action='store_true',
                        help='Run food at the CA food gate fee and sludge at the CA landfill fee')
    options = parser.parse_args()
    if options.both_feedstocks:
        food_snapshots, food_workbook = run_ca_food_fee_test()
        snapshots, workbook = run_ca_landfill_fee_test()
        print(f'Food workbook: {food_workbook}')
        print(f'Sludge workbook: {workbook}')
    else:
        snapshots, workbook = (run_ca_food_fee_test() if options.ca_food_fee_test else
                               run_ca_landfill_fee_test() if options.ca_landfill_fee_test
                               else run_comparison())
