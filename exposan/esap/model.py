"""Discounted-cash-flow minimum selling prices for ESAP products."""

import math
import qsdsan as qs

__all__ = ('calculate_msp',)


def calculate_msp(system, *, discount_rate, lifetime, operating_hours,
                  tea_kwargs=None):
    """Simulate an ESAP system and solve its final-product price at zero NPV.

    Parameters
    ----------
    system : qsdsan.System
        System returned by ``exposan.esap._systems.create_system``. Existing
        feed and coproduct prices are retained. Only ProductConverter's outlet
        is priced as the target product; mother liquor is not saleable product.
    discount_rate : float
        Target annual return as a fraction (e.g., 0.10).
    lifetime : int
        Operating lifetime, in years.
    operating_hours : float
        Annual operating hours, greater than zero and at most 8760.
    tea_kwargs : dict, optional
        Additional ``qsdsan.TEA`` assumptions, such as ``income_tax``,
        ``start_year``, ``CAPEX``, ``annual_labor``, ``annual_maintenance``,
        ``system_add_OPEX``, and ``construction_schedule``. Unspecified values
        use QSDsan defaults: equipment-based capital, zero labor/maintenance,
        no tax, and straight-line depreciation. Cash-flow options default to
        no debt, no working capital, and no startup period. These are preliminary
        assumptions, not a complete project cost estimate. ``uptime_ratio``,
        ``simulate_system``, ``discount_rate``, ``IRR``, ``lifetime``, and
        ``operating_days`` must not be supplied here.

    Returns
    -------
    dict
        ``MSP_USD_per_kg_product`` includes water in the final sold stream
        (including 60 wt% lactate solution). ``MSP_USD_per_kg_N`` divides the
        same product revenue by elemental nitrogen in that stream, excluding
        nitrogen lost to mother liquor. Annual production, residual NPV, and
        the TEA object are also returned for inspection.

    Notes
    -----
    Assigns the solved selling price to the product stream and attaches a new
    TEA to the system. Uses the active utility prices. Upstream recovery costs,
    waste disposal, packaging, and freight are included only if explicitly
    represented by stream prices or additional OPEX supplied by the caller.
    """
    if not math.isfinite(operating_hours) or not 0 < operating_hours <= 8760:
        raise ValueError('operating_hours must be greater than 0 and at most 8760.')
    if not math.isfinite(discount_rate) or discount_rate < 0:
        raise ValueError('discount_rate must be finite and nonnegative.')
    if isinstance(lifetime, bool) or int(lifetime) != lifetime or lifetime <= 0:
        raise ValueError('lifetime must be a positive integer.')
    options = dict(tea_kwargs or {})
    reserved = {'uptime_ratio', 'simulate_system', 'discount_rate', 'IRR',
                'lifetime', 'operating_days'}
    conflicts = reserved.intersection(options)
    if conflicts:
        raise ValueError(f'Use explicit arguments instead of: {sorted(conflicts)}')
    # Set these explicitly to avoid inheriting cash-flow options from prior TEAs.
    for key, value in dict(
            startup_months=0, startup_FOCfrac=1, startup_VOCfrac=1,
            startup_salesfrac=1, WC_over_FCI=0, finance_interest=0,
            finance_years=0, finance_fraction=0).items():
        options.setdefault(key, value)
    # ESAP includes native BioSTEAM equipment. Declare zero additional
    # unit OPEX for that equipment so QSDsan TEA can account for this flowsheet;
    # utility/material costs remain calculated by BioSTEAM as usual.
    for unit in system.units:
        if not hasattr(unit, 'add_OPEX'):
            unit.add_OPEX = {}
        if not hasattr(unit, 'uptime_ratio'):
            unit.uptime_ratio = operating_hours / 8760
    product = system.flowsheet.unit.ProductConverter.outs[0]
    tea = qs.TEA(system, discount_rate=discount_rate, lifetime=int(lifetime),
                 simulate_system=True, **options)
    # Assign after construction so the annualization basis is explicit.
    tea.uptime_ratio = operating_hours / 8760
    product_mass = float(product.F_mass)
    nitrogen_mass = float(product.get_atomic_flow('N') * 14.0067)
    if not math.isfinite(product_mass) or product_mass <= 0:
        raise ValueError('The final product must have positive finite mass flow.')
    if not math.isfinite(nitrogen_mass) or nitrogen_mass <= 0:
        raise ValueError('The final product must contain nitrogen.')
    price = float(tea.solve_price(product))
    if not math.isfinite(price):
        raise RuntimeError('The MSP solver returned a nonfinite price.')
    product.price = price
    return {
        'MSP_USD_per_kg_product': price,
        'MSP_USD_per_kg_N': price * product_mass / nitrogen_mass,
        'annual_product_kg': product_mass * tea.operating_hours,
        'annual_N_kg': nitrogen_mass * tea.operating_hours,
        'NPV_USD': float(tea.NPV),
        'tea': tea,
    }
