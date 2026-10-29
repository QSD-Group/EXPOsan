# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Shared characterization, HTL-aligned economics and life-cycle inventory.

Mass is kg, time is operating hour, gas volume is Nm3 (0 C, 1 atm),
capacity is wet metric tonnes/day during operation, and money is constant USD.
"""
from dataclasses import dataclass, field, asdict
from math import isfinite
import json
import warnings


def positive(name, value, allow_zero=False):
    if not isfinite(value) or value < 0 or (not allow_zero and value == 0):
        raise ValueError(f'{name} must be finite and {"nonnegative" if allow_zero else "positive"}')


def fraction(name, value, upper=1.):
    if not isfinite(value) or not 0 <= value <= upper:
        raise ValueError(f'{name} must be between 0 and {upper}')


@dataclass(frozen=True)
class Feedstock:
    """Empirical characterization; BMP is achieved Nm3 CH4/kg incoming VS.

    Prices and marketable fraction encode contamination, maturity, nutrient
    quality and access to a permitted outlet. They do not certify product quality.
    """
    name: str = 'food'
    total_solids: float = 0.30
    vs_fraction: float = 0.8845
    methane_yield: float = 0.35
    carbon_in_vs: float = 0.50
    product_price_usd_t: float = 0.
    marketable_fraction: float = 1.
    composition: dict = field(default_factory=dict)
    ad_biodegradability: float | None = None
    ad_reactor_conversion: float | None = None
    characterization_source: str = ''
    carbohydrate_fractions: dict = field(default_factory=dict)
    ad_biodegradability_metadata: dict = field(default_factory=dict)

    def __post_init__(self):
        for key in ('total_solids', 'vs_fraction', 'carbon_in_vs', 'marketable_fraction'):
            fraction(key, getattr(self, key))
        positive('methane_yield', self.methane_yield, True)
        if not isfinite(self.product_price_usd_t):
            raise ValueError('Product price must be finite')
        for key in ('ad_biodegradability','ad_reactor_conversion'):
            if getattr(self,key) is not None:
                fraction(key,getattr(self,key))
        if self.composition:
            if (any(not isfinite(v) or v < 0 for v in self.composition.values()) or
                abs(sum(self.composition.values())-1) > 1e-9):
                raise ValueError('composition must contain nonnegative wet fractions summing to one')
            if abs(self.composition['Water']-(1-self.total_solids)) > 1e-9:
                raise ValueError('Composition water must match total_solids')
            if abs(self.composition['Ash']-self.total_solids*(1-self.vs_fraction)) > 1e-9:
                raise ValueError('Composition ash must match vs_fraction')

    @classmethod
    def from_name(cls, name, *, characterization=None, **kwargs):
        """Read the current canonical record; never cache composition by name."""
        from ._composition import characterize
        record, composition = characterize(name, characterization)
        from ._composition import feedstock_module
        resolved = feedstock_module().get_recovery_composition(record,resolve_fibers=True)
        fibers = {key:resolved[key] for key in ('Cellulose','Hemicellulose','Lignin') if key in resolved}
        ts = 1-composition['Water']
        if ts <= 0:
            raise ValueError('Feedstock must contain dry solids')
        return cls(name=name, total_solids=ts, vs_fraction=1-composition['Ash']/ts,
                   composition=composition,
                   ad_biodegradability=record.get('ad_biodegradability'),
                   ad_reactor_conversion=record.get('ad_reactor_conversion'),
                   characterization_source=record.get('characterization_source', ''),
                   carbohydrate_fractions=fibers,
                   ad_biodegradability_metadata={key:record.get('ad_biodegradability_'+key,'')
                       for key in ('source','substrate','basis','status')},
                   **kwargs)


# https://www.nrel.gov/docs/fy13osti/57082.pdf section 6.1.4:
# food TS=30%, VS/TS=88.45%. Other presets and achieved BMPs are scenarios.
FEEDSTOCKS = {
    'food': Feedstock(),
    'manure': Feedstock('manure', .08, .80, .23),
    'sludge': Feedstock('sludge', .20, .70, .22, marketable_fraction=0.),
    'green': Feedstock('green', .50, .85, .12),
}


@dataclass(frozen=True)
class Economics:
    operating_days: float = 330.
    lifetime: int = 30
    discount_rate: float = .10
    income_tax: float = .21
    target_cepci: float = 596.2  # Constant 2020 USD; set explicitly for another basis.
    # Scenario adders on installed cost; zero avoids duplicating bundled study costs.
    engineering_fraction: float = 0.
    contingency_fraction: float = 0.
    working_capital_fraction: float = .05
    maintenance_fraction: float = .03
    insurance_fraction: float = .007
    labor_burden: float = .90
    labor_usd_year: float = 150_000.
    other_opex_usd_t: float = 5.
    feedstock_price_usd_t: float = -40.  # Negative is the single gate-fee income.
    electricity_buy_usd_kwh: float = .10
    electricity_sell_usd_kwh: float = .07
    heat_buy_usd_mj: float = .01
    heat_sell_usd_mj: float = 0.
    diesel_usd_l: float = 1.
    water_usd_t: float = 1.
    residue_disposal_usd_t: float = 40.
    liquid_disposal_usd_t: float = 5.
    land_usd_acre: float = 20_000.
    land_lease_usd_acre_year: float = 0.
    purchase_land: bool = True
    site_development_usd_acre: float = 0.
    start_year: int = 2020
    depreciation: str = 'MACRS7'
    construction_schedule: tuple = (.08, .60, .32)
    startup_months: float = 6.
    startup_FOCfrac: float = 1.
    startup_VOCfrac: float = .75
    startup_salesfrac: float = .5
    finance_interest: float = .08
    finance_years: int = 10
    finance_fraction: float = .60

    def __post_init__(self):
        positive('operating_days', self.operating_days)
        if self.operating_days > 365 or not isinstance(self.lifetime, int) or self.lifetime < 1:
            raise ValueError('Use <=365 operating days and a positive integer lifetime')
        for key in ('income_tax', 'engineering_fraction', 'contingency_fraction',
                    'working_capital_fraction', 'maintenance_fraction', 'insurance_fraction','labor_burden'):
            fraction(key, getattr(self, key))
        for key, value in asdict(self).items():
            if key in ('feedstock_price_usd_t', 'purchase_land', 'depreciation',
                       'construction_schedule'):
                continue
            positive(key, value, key != 'target_cepci')
        if not isfinite(self.feedstock_price_usd_t):
            raise ValueError('Feedstock price must be finite')
        if self.purchase_land and self.land_lease_usd_acre_year:
            raise ValueError('Choose land purchase or lease, not both')
        if (not self.construction_schedule or
            any(not isfinite(i) or i < 0 for i in self.construction_schedule) or
            abs(sum(self.construction_schedule)-1) > 1e-9):
            raise ValueError('construction_schedule must be nonnegative and sum to one')
        if not 0 <= self.startup_months <= 12:
            raise ValueError('startup_months must be between zero and twelve')
        for key in ('startup_FOCfrac','startup_VOCfrac','startup_salesfrac','finance_fraction'):
            fraction(key, getattr(self,key))
        if not isinstance(self.finance_years, int) or not 1 <= self.finance_years <= self.lifetime:
            raise ValueError('finance_years must fit within the project lifetime')

    @classmethod
    def from_htl(cls, tea=None, **overrides):
        """Copy financial policy from an HTL TEA, or its current source defaults.

        Plant-specific labor, equipment factors and land are not transferred.
        Pass the actual HTL TEA to match any run-specific financial overrides.
        """
        from ._finance import htl_financial_assumptions
        values = htl_financial_assumptions(tea)
        values.update(overrides)
        return cls(**values)


@dataclass(frozen=True)
class Equipment:
    name: str
    reference_size: float
    size_unit: str
    installed_usd: float
    base_cepci: float
    exponent: float = .60
    source: str = ''

    def scale(self, size, target_cepci):
        positive('equipment size', size)
        positive('reference_size', self.reference_size)
        positive('base_cepci', self.base_cepci)
        positive('target_cepci', target_cepci)
        positive('installed_usd', self.installed_usd, True)
        fraction('exponent', self.exponent)
        ratio = size / self.reference_size
        if ratio < .1 or ratio > 10:
            warnings.warn(f'{self.name}: scaling {ratio:.2g} times reference; obtain a new cost basis',
                          stacklevel=2)
        return dict(name=self.name, size=size, size_unit=self.size_unit,
                    reference_size=self.reference_size, exponent=self.exponent,
                    reference_installed_usd=self.installed_usd,
                    reference_cepci=self.base_cepci, target_cepci=target_cepci,
                    installed_usd=self.installed_usd * ratio**self.exponent * target_cepci/self.base_cepci,
                    source=self.source)


@dataclass
class Result:
    technology: str
    feedstock: Feedstock
    economics: Economics
    capacity_tpd: float
    equipment: list
    streams_kg_h: dict
    inventory_annual: dict
    design: dict
    land_acres: float
    revenues: dict = field(default_factory=dict)
    extra_costs: dict = field(default_factory=dict)
    notes: list = field(default_factory=list)

    @property
    def annual_tonnes(self):
        return self.capacity_tpd * self.economics.operating_days

    @property
    def hours(self):
        return self.economics.operating_days * 24

    def tea(self, *, native_system=None, return_tea=False):
        """Constant-dollar TEA using the same cash-flow engine as HTL.

        Construction, startup, depreciation, tax-loss carryforward and debt
        follow Economics. No replacement capex, incentives or inflation.
        """
        e, q, inv = self.economics, self.annual_tonnes, self.inventory_annual
        installed = sum(i['installed_usd'] for i in self.equipment)
        if native_system is not None:
            native_installed = native_system.installed_equipment_cost
            if abs(native_installed-installed)>max(1e-6,installed*1e-9):
                raise ValueError('Native installed capital differs from the screening correlation')
            installed = native_installed
        fixed = installed*(1+e.engineering_fraction+e.contingency_fraction)
        fixed += self.land_acres*e.site_development_usd_acre
        land = self.land_acres*e.land_usd_acre if e.purchase_land else 0.
        wc = fixed*e.working_capital_fraction
        costs = dict(maintenance=installed*e.maintenance_fraction,
                     insurance=fixed*e.insurance_fraction, labor=e.labor_usd_year*(1+e.labor_burden),
                     other=q*e.other_opex_usd_t,
                     feedstock=max(e.feedstock_price_usd_t, 0)*q,
                     electricity=inv.get('electricity_import_kwh', 0)*e.electricity_buy_usd_kwh,
                     heat=inv.get('heat_import_mj', 0)*e.heat_buy_usd_mj,
                     diesel=inv.get('diesel_l', 0)*e.diesel_usd_l,
                     water=inv.get('water_kg', 0)/1000*e.water_usd_t,
                     residue=inv.get('residue_kg', 0)/1000*e.residue_disposal_usd_t,
                     liquid=inv.get('liquid_kg', 0)/1000*e.liquid_disposal_usd_t,
                     land_lease=0. if e.purchase_land else self.land_acres*e.land_lease_usd_acre_year,
                     **self.extra_costs)
        revenue = dict(gate_fee=max(-e.feedstock_price_usd_t, 0)*q,
                       electricity=-inv.get('electricity_export_kwh', 0)*e.electricity_sell_usd_kwh,
                       heat=-inv.get('heat_export_mj', 0)*e.heat_sell_usd_mj,
                       **self.revenues)
        opex, sales = sum(costs.values()), sum(revenue.values())
        ebitda = sales-opex
        from ._finance import cashflow_analysis
        fixed_opex = sum(costs[k] for k in ('maintenance','insurance','labor','land_lease'))
        finance, native_tea = cashflow_analysis(e, fixed, land, fixed_opex, opex-fixed_opex, sales,
                                               native_system=native_system,return_tea=True)
        finance['financial_assumptions'].update(target_cepci=e.target_cepci,
            electricity_buy_usd_kwh=e.electricity_buy_usd_kwh,
            electricity_sell_usd_kwh=e.electricity_sell_usd_kwh,
            installed_capital_basis='bundled installed correlations; no additional equipment BM factor')
        r, n = e.discount_rate, e.lifetime
        crf = r/(1-(1+r)**-n) if r else 1/n
        annualized_capital = (fixed+land+wc-(land+wc)/(1+r)**n)*crf
        analysis = dict(installed_cost=installed, fixed_capital=fixed, land=land,
                    working_capital=wc, total_investment=fixed+land+wc,
                    annual_costs=costs, annual_revenues=revenue, annual_opex=opex,
                    annual_revenue=sales, EBITDA=ebitda,
                    **finance,
                    annualized_capital=annualized_capital,
                    levelized_net_cost_usd_t=(annualized_capital+opex-sales)/q,
                    net_operating_cost_usd_t=(opex-sales)/q)
        return (analysis,native_tea) if return_tea else analysis

    def lca(self, factors, *, avoided_disposal_kg_co2e_t=0.):
        """Characterize inventory using explicit kg CO2e/activity-unit factors.

        Every nonzero activity requires a factor, including explicit zeros.
        Export quantities are negative credits. Biogenic CO2 is separately
        inventoried. Supply consistent electricity displacement and GWP factors.
        """
        missing = [k for k,v in self.inventory_annual.items() if v and k not in factors]
        if missing:
            raise ValueError(f'Missing characterization factors: {missing}')
        if any(not isfinite(v) for v in factors.values()) or not isfinite(avoided_disposal_kg_co2e_t):
            raise ValueError('Characterization factors must be finite')
        impacts = {k:v*factors.get(k, 0.) for k,v in self.inventory_annual.items()}
        impacts['avoided_disposal'] = -self.annual_tonnes*avoided_disposal_kg_co2e_t
        total = sum(impacts.values())
        energy = self.design.get('exported_energy_mj_year', 0.)
        product = self.design.get('marketable_product_kg_year', 0.)
        return dict(annual_kg_co2e=total, breakdown=impacts,
                    kg_co2e_t_wet=total/self.annual_tonnes,
                    g_co2e_mj=total*1000/energy if energy else None,
                    kg_co2e_kg_product=total/product if product else None)

    def to_json(self):
        return json.dumps(dict(**asdict(self), tea=self.tea()), indent=2, allow_nan=False)

    def to_qsdsan(self, factors, **kwargs):
        """Build native streams, utilities and LCA; requires BioSTEAM/QSDsan."""
        from ._native import to_qsdsan
        return to_qsdsan(self, factors, **kwargs)
