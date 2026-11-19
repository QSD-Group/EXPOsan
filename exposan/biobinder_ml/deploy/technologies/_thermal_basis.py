# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Explicit cost, market and operating bases for thermochemical systems."""
from dataclasses import dataclass, field, asdict
from math import isfinite
from ._common import positive, fraction


@dataclass(frozen=True)
class CostBasis:
    """Equipment quotation/correlation; FX is USD per source-currency unit.

    Apply FX in the quotation year, then CEPCI to the target cost year.
    ``installation_factor`` applies only to purchase-cost quotations.
    """
    value: float
    reference_size: float
    size_unit: str
    cost_year: int
    reference_cepci: float
    source: str
    locator: str
    exponent: float
    kind: str = 'purchase'
    installation_factor: float = 1.
    currency: str = 'USD'
    usd_per_currency: float = 1.
    fx_source: str = ''
    fx_year: int | None = None
    scaling_basis: str = 'published exponent'
    scope: str = ''

    def __post_init__(self):
        for key in ('value', 'reference_size', 'reference_cepci', 'usd_per_currency'):
            positive(key, getattr(self, key))
        fraction('exponent', self.exponent)
        if self.installation_factor < 1 or not isfinite(self.installation_factor):
            raise ValueError('installation_factor must be finite and >=1')
        if self.kind not in ('purchase', 'installed'):
            raise ValueError('Cost kind must be purchase or installed')
        if self.kind == 'installed' and self.installation_factor != 1:
            raise ValueError('Installed cost cannot receive another installation factor')
        if not self.source.startswith('https://') or not self.locator:
            raise ValueError('Equipment cost needs an active source URL and exact locator')
        if self.currency != 'USD' and (not self.fx_source.startswith('https://') or self.fx_year != self.cost_year):
            raise ValueError('Foreign costs require source-year USD/FX rate and source URL')
        if self.currency == 'USD' and self.usd_per_currency != 1:
            raise ValueError('USD cost must have a unit currency conversion')

    def estimate(self, size, target_cepci, target_year):
        positive('equipment size', size)
        positive('target_cepci', target_cepci)
        ratio = size/self.reference_size
        amount = self.value*self.usd_per_currency*ratio**self.exponent*target_cepci/self.reference_cepci
        purchase = amount if self.kind == 'purchase' else None
        return dict(**asdict(self), size=size, target_year=target_year,
                    target_cepci=target_cepci, scale_ratio=ratio,
                    purchase_usd=purchase,
                    installed_usd=amount*self.installation_factor,
                    applicability='scaled reference; vendor confirmation required' if ratio < .1 or ratio > 10 else 'reference scaling range not independently validated')


@dataclass(frozen=True)
class ProductMarket:
    """An explicit outlet and price basis, with quality requirements.

    Price is per kg whole product or per MJ HHV. Price-year adjustment uses an
    explicit commodity/deflator factor; equipment CEPCI is never used for sales.
    Documentary qualification covers unmodeled contaminants and buyer acceptance.
    """
    price: float
    unit: str
    currency: str
    price_year: int
    target_year: int
    source: str
    outlet: str
    basis: str = 'qualified_contract'
    usd_per_currency: float = 1.
    fx_source: str = ''
    fx_year: int | None = None
    year_factor: float = 1.
    year_factor_source: str = ''
    qualified: bool = False
    qualification_source: str = ''
    saleable_fraction: float = 1.
    min_purity_mol: float | None = None
    max_water_wt: float | None = None
    max_ash_dry_wt: float | None = None
    max_nitrogen_dry_wt: float | None = None
    max_sulfur_dry_wt: float | None = None
    max_co_mol: float | None = None
    min_hhv_mj_kg: float | None = None

    def __post_init__(self):
        for key in ('price', 'usd_per_currency', 'year_factor'):
            positive(key, getattr(self, key), key == 'price')
        fraction('saleable_fraction', self.saleable_fraction)
        if self.unit not in ('USD/kg', 'USD/MJ-HHV', 'currency/kg', 'currency/MJ-HHV'):
            raise ValueError('Price unit must identify kg product or MJ-HHV')
        if self.basis not in ('qualified_contract','energy_proxy','literature_proxy'):
            raise ValueError('Market basis must distinguish contracts and price proxies')
        if not self.source.startswith('https://') or not self.outlet:
            raise ValueError('Market needs a source URL and plant-gate outlet/boundary')
        if self.currency != 'USD' and (not self.fx_source.startswith('https://') or self.fx_year != self.price_year):
            raise ValueError('Foreign sale price requires price-year FX source and rate')
        if self.currency == 'USD' and self.usd_per_currency != 1:
            raise ValueError('USD sale price must use unit currency conversion')
        if self.currency != 'USD' and self.unit.startswith('USD/'):
            raise ValueError('Foreign price units must start with currency/, before FX conversion')
        if self.price_year != self.target_year and not self.year_factor_source.startswith('https://'):
            raise ValueError('Different price years need an explicit sourced price-year factor')
        if self.price_year == self.target_year and self.year_factor != 1:
            raise ValueError('Same-year price must use year_factor=1')
        if self.qualified and not self.qualification_source.startswith('https://'):
            raise ValueError('Qualification needs a buyer/specification evidence URL')
        for key in ('min_purity_mol', 'max_water_wt', 'max_ash_dry_wt',
                    'max_nitrogen_dry_wt', 'max_sulfur_dry_wt', 'max_co_mol'):
            value = getattr(self, key)
            if value is not None:
                fraction(key, value)
        if self.min_hhv_mj_kg is not None:
            positive('min_hhv_mj_kg', self.min_hhv_mj_kg, True)


def default_product_markets(product_names, target_year=2020, saleable_fraction=1.):
    """Return explicitly conditional energy/literature revenue proxies.

    Passing ``markets={}`` to a factory selects the zero-product-credit scenario.
    Conversion of price years uses CPI, never equipment CEPCI.
    """
    if target_year != 2020:
        raise ValueError('Default price proxies are in 2020 USD; supply explicit markets for another cost year')
    markets = {}
    for name in product_names:
        if name == 'hydrochar':
            # https://www.mdpi.com/1996-1073/10/2/211 Section3.6.3:
            # 200 EUR/t break-even proxy, economic basis2015, HHV25.64 GJ/t.
            # Use 200/25.64=7.800312 EUR/GJ; the paper also prints8.3 EUR/GJ,
            # which does not reconcile with its quoted mass price and HHV.
            # This is a modeled break-even proxy, not an observed sales contract.
            # https://www.bundesbank.de/resource/blob/966466/83893e360b65c144492a0bf32a7fb6b0/472B63F073F071307366337C94F8C870/2025-09-18-09-05-16-lange-zeitreihen-data.pdf
            # Annual2015 EUR/USD =1.1095 USD/EUR.
            # https://www.bls.gov/regions/mid-atlantic/data/ConsumerPriceIndexAnnualandSemiAnnual_Table.htm
            # CPI-U2015=237.017,2020=258.811; constant-USD proxy adjustment.
            markets[name] = ProductMarket(price=200/1000/25.64,unit='currency/MJ-HHV',
                currency='EUR',price_year=2015,target_year=2020,
                source='https://www.mdpi.com/1996-1073/10/2/211',
                outlet='Conditional solid-fuel hydrochar proxy: 200 EUR/t at25.64 MJ/kg HHV; modeled2015 break-even, adjusted to actual product HHV including water/ash dilution; no pellet or contaminant certification',
                basis='literature_proxy',usd_per_currency=1.1095,fx_year=2015,
                fx_source='https://www.bundesbank.de/resource/blob/966466/83893e360b65c144492a0bf32a7fb6b0/472B63F073F071307366337C94F8C870/2025-09-18-09-05-16-lange-zeitreihen-data.pdf',
                year_factor=258.811/237.017,
                year_factor_source='https://www.bls.gov/regions/mid-atlantic/data/ConsumerPriceIndexAnnualandSemiAnnual_Table.htm',
                saleable_fraction=saleable_fraction)
        else:
            # https://www.eia.gov/dnav/ng/hist/rngwhhda.htm : Henry Hub2020
            # annual2.03 USD/MMBtu.1 MMBtu=1055.05585262 MJ (IT Btu).
            # https://www.eia.gov/energyexplained/units-and-calculators/british-thermal-units.php
            # Energy-parity scenario for oil/gases/biochar; no petroleum, SAF,
            # hydrogen-purity, soil-amendment or renewable-credit premium.
            markets[name] = ProductMarket(price=2.03/1055.05585262,unit='USD/MJ-HHV',
                currency='USD',price_year=2020,target_year=2020,
                source='https://www.eia.gov/dnav/ng/hist/rngwhhda.htm',
                outlet='Conditional plant-gate fuel-energy parity with2020 Henry Hub natural gas; modeled whole-product HHV; no delivery, connection, qualification or renewable-credit premium',
                basis='energy_proxy',saleable_fraction=saleable_fraction)
    return markets


@dataclass(frozen=True)
class ThermalConfig:
    """Process scenario; retentions are fractions of each feed element.

    Defaults are deliberately uncalibrated elemental-allocation scenarios.
    They support construction and balance checks, not validated waste conversion.
    """
    # https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-23053.pdf
    # Fast-pyrolysis design temperature 500 C; atmospheric pressure is a scenario.
    temperature_k: float = 773.15
    pressure_pa: float = 101325.
    residence_time_h: float = .5  # h; generic scenario, specialized by each technology.
    # https://www.mdpi.com/1996-1073/10/2/211 : wet-slurry/centrifuge/dryer topology.
    # The fractions below are engineering scenarios, NOT reported waste measurements.
    slurry_solids: float = .15  # kg dry matter/kg reactor slurry.
    dryer_moisture: float = .10  # kg water/kg pyrolysis feed after drying.
    cake_moisture: float = .50  # kg water/kg separated cake.
    product_moisture: float = .10  # finishing target; no fuel-grade certification.
    solids_capture: float = .99  # specified recovery; no particle-size model.
    # https://www.psi.ch/en/cpe/hydropilot : reforming/methanation process basis.
    # Conversion, extent, CO allocation and separation below are uncalibrated
    # scenario specifications, not PSI performance measurements or equilibrium.
    gasification_conversion: float = .90
    methanation_completion: float = .97
    co_carbon_fraction: float = .02
    gas_recovery: float = .90
    impurity_recovery: float = .001
    product_pressure_pa: float = 3e6  # Pa; scenario delivery pressure.
    separation_pressure_pa: float = 2e6  # Pa; scenario separator/PSA pressure.
    salt_capture: float = .85  # scenario fraction of aggregate ash to brine.
    brine_water_per_kg_ash: float = 4.  # kg/kg; scenario brine dilution.
    # https://www.mdpi.com/1996-1073/10/2/211 Section3.6: 75% reactor filling.
    working_volume_fraction: float = .75
    max_vessel_m3: float = 100.  # m3; engineering modularization scenario.
    land_acres: float = 3.  # acres; common scenario, not a site plot-plan estimate.
    cost_year: int = 2020
    conversion_source: str = ''
    conversion_basis: str = 'uncalibrated elemental-retention scenario'
    # https://www.ieabioenergy.com/wp-content/uploads/2015/02/XII1-Oasma-Bio_oil-is-not-Bio_oil.pdf
    # Phenolic oil family motivates a guaiacol latent-heat surrogate only. This
    # does not assert guaiacol composition or measured waste-oil phase behavior.
    oil_latent_reference: str = 'Guaiacol'
    ultimate_daf: dict = field(default_factory=dict)
    # https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-23053.pdf
    # Table5 is low-ash wood only. The following elemental retentions are explicitly
    # uncalibrated scenarios rather than transferring its wood yields to wet wastes.
    oil_retention: dict = field(default_factory=lambda: dict(C=.60, H=.60, O=.45, N=.60, S=.50))
    char_retention: dict = field(default_factory=lambda: dict(C=.20, H=.05, O=.08, N=.25, S=.35))
    aqueous_retention: dict = field(default_factory=dict)
    cost_overrides: dict = field(default_factory=dict)

    def __post_init__(self):
        for key in ('temperature_k', 'pressure_pa', 'residence_time_h', 'slurry_solids',
                    'product_pressure_pa', 'separation_pressure_pa', 'max_vessel_m3'):
            positive(key, getattr(self, key))
        for key in ('dryer_moisture', 'cake_moisture', 'product_moisture', 'solids_capture',
                    'gasification_conversion', 'methanation_completion', 'co_carbon_fraction',
                    'gas_recovery', 'impurity_recovery', 'salt_capture',
                    'working_volume_fraction', 'slurry_solids'):
            fraction(key, getattr(self, key))
        if self.slurry_solids == 0 or self.working_volume_fraction == 0:
            raise ValueError('Working volume and slurry solids must be positive')
        if max(self.dryer_moisture, self.cake_moisture, self.product_moisture) >= 1:
            raise ValueError('Moisture targets must be below one')
        if self.product_moisture > self.cake_moisture:
            raise ValueError('Product drying target must not exceed cake moisture')
        for key in ('land_acres', 'brine_water_per_kg_ash'):
            positive(key, getattr(self, key), True)
        for mapping in (self.oil_retention,self.char_retention,self.aqueous_retention):
            if set(mapping)-set(('C','H','O','N','S')):
                raise ValueError('Element retentions support CHONS only')
        for element in ('C', 'H', 'O', 'N', 'S'):
            retentions = [mapping.get(element, 0.) for mapping in
                          (self.oil_retention, self.char_retention, self.aqueous_retention)]
            for retention in retentions:
                fraction('element retention', retention)
            if sum(retentions) > 1+1e-12:
                raise ValueError(f'{element} product retentions exceed feed availability')
        if self.conversion_source and not self.conversion_source.startswith('https://'):
            raise ValueError('Conversion source must be an active URL')
        if self.ultimate_daf:
            if set(self.ultimate_daf) - set(('C', 'H', 'O', 'N', 'S')):
                raise ValueError('Ultimate DAF supports CHONS only; ash is separate')
            if any(not isfinite(v) or v < 0 for v in self.ultimate_daf.values()) or abs(sum(self.ultimate_daf.values())-1) > 1e-8:
                raise ValueError('Ultimate DAF must contain CHONS mass fractions summing to one')
            positive('Ultimate carbon', self.ultimate_daf.get('C', 0.))
