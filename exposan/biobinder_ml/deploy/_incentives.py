# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

"""Opt-in HTL incentive scenarios for BioSTEAM/QSDsan.

All money is USD, volumes are physical US gallons (never GGE), energy is
MJ LHV, and CI is g CO2e/MJ. Approved pathway CI is a separate input from
a system-level LCA result.

Use ``cls=with_incentives(BiobinderTEA)`` when calling ``create_tea``;
then set ``tea.incentive_schedule = lambda: {2027: annual_support, ...}``.
The callable is evaluated on every TEA calculation, including solve_price,
so it can read current simulated flows. Values are realized AFTER-TAX cash
receipts; do not also add them to stream prices or taxable sales. Grants
requiring basis/depreciation changes need a project-specific TEA subclass.
No financial helper changes LCA characterization factors or inventories.
"""

from dataclasses import dataclass
from math import floor, isfinite

__all__ = (
    'annual_quantity', 'production_support', 'rin_support',
    'clean_fuel_support', 'section45z_support', 'capital_grant',
    'Support', 'support_schedule', 'with_incentives',
    'set_financing', 'lca_carbon_intensity',
)


def _nonnegative(value, name):
    value = float(value)
    if not isfinite(value) or value < 0:
        raise ValueError(f'{name} must be finite and nonnegative')
    return value


def _fraction(value, name):
    value = _nonnegative(value, name)
    if value > 1:
        raise ValueError(f'{name} must be between zero and one')
    return value


def annual_quantity(stream, tea, *, basis='kg', eligible_fraction=1.,
                    lhv_MJ_per_kg=None):
    """Annual output from current kg/h or m3/h stream flow and TEA hours.

    ``basis`` is 'kg', 'US_gal', or 'MJ_LHV'. LHV must be supplied explicitly
    for energy output; do not substitute a stream's HHV. Eligibility and
    destination allocation are represented by eligible_fraction.
    """
    hours = _nonnegative(tea.operating_hours, 'operating_hours')
    fraction = _fraction(eligible_fraction, 'eligible_fraction')
    if basis == 'US_gal':
        hourly = _nonnegative(stream.F_vol, 'F_vol') / 0.003785411784
    elif basis in ('kg', 'MJ_LHV'):
        hourly = _nonnegative(stream.F_mass, 'F_mass')
        if basis == 'MJ_LHV':
            if lhv_MJ_per_kg is None:
                raise ValueError('MJ_LHV requires lhv_MJ_per_kg')
            hourly *= _nonnegative(lhv_MJ_per_kg, 'lhv_MJ_per_kg')
    else:
        raise ValueError("basis must be 'kg', 'US_gal', or 'MJ_LHV'")
    return hourly * hours * fraction


# USDA advanced biofuel production payments:
# https://www.rd.usda.gov/programs-services/energy-programs/advanced-biofuel-payment-program
# Examples of state product credits (check product, recipient and award caps):
# https://www.revenue.state.mn.us/sustainable-aviation-fuel-credit
# https://tax.illinois.gov/forms/sales/aviation/st-70-instructions.html
# https://opportunityiowa.gov/business/financial-assistance/sustainable-aviation-fuel-production-tax-credit-program
# https://legislature.maine.gov/statutes/36/title36sec5219-XX.html
def production_support(quantity, rate, *, eligible=False, retained_share=1.,
                       realization=1., cap=None):
    """Annual receipts; quantity and $/quantity rate must use matching units.

    Cap applies to gross program award before ownership/transfer discounts.
    For taxable revenue, account for taxes before passing receipts to TEA.
    """
    gross = _nonnegative(quantity, 'quantity') * _nonnegative(rate, 'rate')
    if cap is not None:
        gross = min(gross, _nonnegative(cap, 'cap'))
    share = _fraction(retained_share, 'retained_share')
    realized = _fraction(realization, 'realization')
    return gross * share * realized if eligible else 0.


# EPA approved feedstock/process/fuel pathways:
# https://www.epa.gov/renewable-fuel-standard-program/approved-pathways-renewable-fuel
def rin_support(gallons, price_per_RIN, equivalence_value, **kwargs):
    """Approved eligible gallons times RIN equivalence times $/RIN.

    No D-code or pathway approval is inferred. kwargs match production_support.
    """
    return production_support(
        gallons, _nonnegative(price_per_RIN, 'price_per_RIN')
        * _nonnegative(equivalence_value, 'equivalence_value'), **kwargs)


# Destination-specific clean fuel programs; use each program's actual rules:
# https://ww2.arb.ca.gov/our-work/programs/low-carbon-fuel-standard/about
# https://www.oregon.gov/deq/ghgp/cfp/pages/cfp-overview.aspx
# https://ecology.wa.gov/air-climate/reducing-greenhouse-gas-emissions/clean-fuel-standard
# https://www.env.nm.gov/clean-transportation-fuel-program/
def clean_fuel_support(energy_MJ, benchmark_ci, pathway_ci,
                       price_per_tonne, **kwargs):
    """Illustrative positive CI-difference revenue, not a compliance engine.

    Caller must apply jurisdiction-specific EER, pathway, destination and
    reporting rules before using this result. Excess CI earns zero here;
    compliance deficits are not modeled by this function.
    """
    if not all(isfinite(float(i)) for i in (benchmark_ci, pathway_ci)):
        raise ValueError('CI values must be finite')
    rate = max(0., benchmark_ci - pathway_ci) / 1e6
    rate *= _nonnegative(price_per_tonne, 'price_per_tonne')
    return production_support(energy_MJ, rate, **kwargs)


# IRS 45Z overview and emissions guidance:
# https://www.irs.gov/credits-deductions/clean-fuel-production-credit
# https://www.irs.gov/pub/irs-drop/n-26-53.pdf
# Proposed implementation regulations (distinguish proposal from enacted law):
# https://www.federalregister.gov/documents/2026/02/04/2026-02246/section-45z-clean-fuel-production-credit
def section45z_support(gallons, emissions_rate, applicable_amount, *,
                       baseline, negative_rate_allowed=False, **kwargs):
    """Illustrative 45Z calculation with nearest-0.1 factor rounding (half-up).

    emissions_rate and baseline must share the approved statutory units
    (not automatically g/MJ). Caller supplies the applicable indexed amount,
    certification and negative-rate permission. Calendar windows belong in
    Support; no current-law dates/rates are hard-coded.
    """
    baseline = _nonnegative(baseline, 'baseline')
    if baseline == 0 or not isfinite(float(emissions_rate)):
        raise ValueError('positive baseline and finite emissions_rate required')
    if emissions_rate < 0 and not negative_rate_allowed:
        raise ValueError('negative emissions rate requires explicit permission')
    # https://www.irs.gov/instructions/i7218
    factor = floor(max(0., (baseline - emissions_rate) / baseline) * 10 + .5) / 10
    return production_support(
        gallons, _nonnegative(applicable_amount, 'applicable_amount') * factor,
        **kwargs)


# Capital assistance / principal-forgiveness screening sources:
# https://www.epa.gov/cwsrf/cwsrf-eligibilities-factsheets
# https://www.epa.gov/cwsrf/cwsrf-emerging-contaminants-frequent-questions-and-answers
# https://www.rd.usda.gov/programs-services/energy-programs/rural-energy-america-program-renewable-energy-systems-energy-efficiency-improvement-guaranteed-loans
def capital_grant(eligible_capital, award_fraction, *, eligible=False, cap=None):
    """One-time gross award; does not reduce FCI, debt, or tax basis.

    Apply once in the actual receipt year. Account separately for grant
    taxation, basis adjustments, eligible asset allocation and debt sizing.
    A loan guarantee is not a capital grant.
    """
    return production_support(eligible_capital,
                              _fraction(award_fraction, 'award_fraction'),
                              eligible=eligible, cap=cap)


@dataclass(frozen=True)
class Support:
    """Ownership ledger entry with inclusive calendar-year receipt window.

    annual_cash is an after-tax USD amount or callable(year) -> amount.
    Distinct IDs identify distinct claims, including recipient/destination;
    duplicate IDs are rejected. Caller verifies stacking across facilities.
    lifetime_cap limits realized cash, not a statutory gross-award cap.
    Model statutory caps before transfer discounts in production_support.
    """
    claim_id: str
    recipient: str
    start_year: int
    end_year: int
    annual_cash: object
    eligible: bool = False
    lifetime_cap: float = None


def support_schedule(supports):
    """Combine enabled ledger claims as {calendar_year: after_tax_USD}.

    Recompute inside a TEA schedule callable for flow-dependent scenarios.
    Expiry is applied to actual calendar years, never years since startup.
    """
    schedule, seen = {}, set()
    for claim in supports:
        if not claim.claim_id or not claim.recipient or claim.claim_id in seen:
            raise ValueError('claims require unique IDs and explicit recipients')
        seen.add(claim.claim_id)
        if (not isinstance(claim.start_year, int)
                or not isinstance(claim.end_year, int)
                or claim.end_year < claim.start_year):
            raise ValueError('claim years must be ordered integers')
        remaining = (float('inf') if claim.lifetime_cap is None else
                     _nonnegative(claim.lifetime_cap, 'lifetime_cap'))
        if not claim.eligible:
            continue
        for year in range(claim.start_year, claim.end_year + 1):
            cash = (claim.annual_cash(year) if callable(claim.annual_cash)
                    else claim.annual_cash)
            cash = min(_nonnegative(cash, 'annual_cash'), remaining)
            schedule[year] = schedule.get(year, 0.) + cash
            remaining -= cash
    return schedule


# BioSTEAM TEA implementation, including the native tax/incentive hook:
# https://github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/biosteam/_tea.py
def with_incentives(base_tea):
    """Return TEA subclass integrating realized receipts into native cash flow.

    Works with the local BioSTEAM five-argument tax/incentive hook and derived
    QSDsan TEAs. Use as create_tea(..., cls=with_incentives(BiobinderTEA)).
    Set incentive_schedule AFTER construction to a mapping or zero-argument
    callable returning a mapping. Default is no support. Native NPV/IRR and
    solve_price use the same receipts. Existing base-class incentives persist.
    The first array year is duration[0] - len(construction_schedule).
    No startup proration is inferred; supply actual first-year eligible output.
    """
    class IncentivizedTEA(base_tea):
        incentive_schedule = None

        def _fill_tax_and_incentives(self, incentives, taxable_cashflow,
                                     nontaxable_cashflow, tax, depreciation):
            super()._fill_tax_and_incentives(
                incentives, taxable_cashflow, nontaxable_cashflow, tax,
                depreciation)
            schedule = self.incentive_schedule
            if callable(schedule):
                schedule = schedule()
            if schedule is None:
                return
            first_year = self.duration[0] - len(self.construction_schedule)
            for year, cash in schedule.items():
                if not isinstance(year, int):
                    raise ValueError('schedule keys must be integer calendar years')
                cash = _nonnegative(cash, 'scheduled cash')
                index = year - first_year
                if 0 <= index < len(incentives):
                    incentives[index] += cash

    IncentivizedTEA.__name__ = f'Incentivized{base_tea.__name__}'
    return IncentivizedTEA


# Financing programs (loans/guarantees are distinct from capital grants):
# https://www.rd.usda.gov/programs-services/energy-programs/biorefinery-renewable-chemical-and-biobased-product-manufacturing-assistance-program
# https://www.epa.gov/cwsrf/about-cwsrf
# https://www.epa.gov/water-infrastructure/water-infrastructure-wifia-program
def set_financing(tea, *, interest, years, fraction):
    """Set concessional debt terms in-place; return previous terms for reset.

    A guarantee changes debt availability/terms, never product revenue.
    fraction is BioSTEAM's project finance fraction, not the lender's federal
    guarantee percentage. Combine multiple tranches outside this helper.
    """
    interest = _nonnegative(interest, 'interest')
    fraction = _fraction(fraction, 'fraction')
    if not isinstance(years, int) or years <= 0:
        raise ValueError('years must be a positive integer')
    previous = dict(interest=tea.finance_interest, years=tea.finance_years,
                    fraction=tea.finance_fraction)
    tea.finance_interest = interest
    tea.finance_years = years
    tea.finance_fraction = fraction
    return previous


# QSDsan LCA implementation: annual totals and product allocation:
# https://github.com/QSD-Group/QSDsan/blob/main/qsdsan/_lca.py
def lca_carbon_intensity(lca, annual_energy_MJ, *, indicator='GWP',
                         allocated_streams=None):
    """Diagnostic g CO2e/MJ using a QSDsan kg-CO2e indicator.

    Without allocated_streams this uses total system annual impacts. Supply
    product streams for QSDsan's energy allocation. Neither result establishes
    a regulatory pathway CI. Indicator units must be kg CO2e; verify them.
    Financial support does not change this physical environmental result.
    """
    energy = _nonnegative(annual_energy_MJ, 'annual_energy_MJ')
    if energy == 0:
        raise ValueError('annual_energy_MJ must be positive')
    if allocated_streams is None:
        impacts = lca.get_total_impacts(annual=True)
    else:
        impacts = lca.get_allocated_impacts(
            streams=allocated_streams, allocate_by='energy', annual=True)
    impact = float(impacts[indicator])
    if not isfinite(impact):
        raise ValueError('LCA impact must be finite')
    return impact * 1000 / energy
