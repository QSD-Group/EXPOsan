"""User-supplied acid and utility prices for preliminary ESAP economics.

Prices are USD/kg purchased product, using the arithmetic mean of the
supplier rows in TEA.xlsx. Only the resulting mean prices are stored here.
Each quote receives equal weight regardless of package size. No freight or
cost-year adjustment is applied. Common supplier-product purities are 99%
maleic, 88% liquid lactic, 96% sulfuric, and 100% citric acid. The 60% lactic
powder with calcium lactate carrier is excluded from the model price inputs.
Workbook edits are not loaded automatically: update these mean prices when the workbook changes.
"""

from qsdsan import CEPCI_by_year

__all__ = ('equipment_cost_year', 'equipment_cost_index', 'acid_prices',
           'acid_purities', 'acid_densities',
           'electricity_price', 'electricity_price_source',
           'gas_escalation_factor', 'electricity_escalation_factor',
           'steam_regeneration_price', 'cooling_water_regeneration_price',
           'chilled_water_heat_transfer_price',
           'natural_gas_price', 'natural_gas_price_per_kmol',
           'natural_gas_price_per_Mcf', 'natural_gas_price_source')

# Temporary equipment-cost basis: latest year in QSDsan's bundled CEPCI table
# (currently 2023: 797.9), NOT a 2026 index. QSDsan's table has no source citation.
# TODO: Replace with a verified Chemical Engineering CEPCI value for the study
# period, documenting its year/month and preliminary/final status.
# https://www.chemengonline.com/pci-home/
equipment_cost_year = max(CEPCI_by_year)
equipment_cost_index = CEPCI_by_year[equipment_cost_year]

# User-selected study value: 14.29 cents/kWh, converted to USD/kWh.
# The EIA table updates over time; no specific region/sector/year was supplied.
electricity_price = 0.1429
electricity_price_source = (
    'https://www.eia.gov/electricity/monthly/epm_table_grapher.php?t=table_5_06_b'
)

# U.S. industrial, July 2026; nominal USD per thousand standard cubic feet.
natural_gas_price_per_Mcf = 4.58
natural_gas_price_source = 'https://www.eia.gov/dnav/ng/ng_pri_sum_dcu_nus_m.htm'
# EIA volume basis: 60 deg F, 14.73 psia (absolute).
# https://www.eia.gov/naturalgas/monthly/pdf/front_matter.pdf
# ESAP's dryer models fuel as pure CH4. Use ideal-gas methane at this
# reference state (Z=1), not volume at the dryer inlet temperature/pressure.
_gas_reference_T = (60 - 32) * 5 / 9 + 273.15  # K
_gas_reference_P = 14.73 * 6894.757293168       # Pa
_methane_MW = 16.04246                        # kg/kmol
_ft3_to_m3 = 0.3048**3                        # m3/ft3
_ideal_gas_constant = 8314.46261815324        # Pa*m3/(kmol*K)
_kmol_per_Mcf = _gas_reference_P * (1000 * _ft3_to_m3) / (
    _ideal_gas_constant * _gas_reference_T
)
natural_gas_price_per_kmol = natural_gas_price_per_Mcf / _kmol_per_Mcf
natural_gas_price = natural_gas_price_per_kmol / _methane_MW  # USD/kg CH4

# Preliminary 2026 service-price escalation; 2016 is an ASSUMED price basis.
# BioSTEAM _heat_utility.py cites Seider (2016), Table 17.1, using Busche (1995)
# data; the textbook publication year does not establish the actual price year.
# https://biosteam.readthedocs.io/en/latest/_modules/biosteam/_heat_utility.html
# EIA 2016 U.S. industrial annual prices: gas 3.51 USD/Mcf; power 0.0676 USD/kWh.
# https://www.eia.gov/dnav/ng/hist/n3035us3a.htm
# https://www.eia.gov/electricity/annual/table.php?t=epa_02_07.html
# Ratios use the selected study prices, not a complete 2026 annual average.
gas_escalation_factor = natural_gas_price_per_Mcf / 3.51
electricity_escalation_factor = electricity_price / 0.0676
# Explicit BioSTEAM baseline constants prevent compounding across system builds.
# Steam follows gas; all cooling-water/chilled-water costs follow electricity
# as a simplifying proxy, including their non-electric water/chemical costs.
steam_regeneration_price = 0.2378 * gas_escalation_factor  # USD/kmol LP steam
cooling_water_regeneration_price = 4.8785e-4 * electricity_escalation_factor  # USD/kmol
chilled_water_heat_transfer_price = 5e-6 * electricity_escalation_factor  # USD/kJ

# TODO: Replace proxy escalation with updated price data. Some resources:
# https://www.energy.gov/sites/prod/files/2014/05/f15/tech_brief_true_cost.pdf
# cooling water: tower makeup/blowdown, chemicals, pump/fan electricity;
# https://www.energy.gov/cmei/femp/best-management-practice-10-cooling-tower-management
# chilled water: chiller COP at operating temperatures, pumps and heat rejection;
# https://betterbuildingssolutioncenter.energy.gov/better-plants/process-cooling-and-hvac
# See BioSTEAM._heat_utility for default price library and utility cost calculation.

# Approved purchased-solution densities, used ONLY for volume package pricing.
# Sulfuric: 1.84 kg/L at 25 C for Lab Alley's nominal 96% (95-98%) product;
# applied uniformly to the sulfuric quotes despite their different assays.
# https://www.laballey.com/products/sulfuric-acid-acs-grade
# Lactic: 1.20 kg/L is an approximation from Sigma's approximately 85% solution
# applied to 88% products without concentration/temperature correction. It is
# within Lab Alley's broad 1.0-1.3 kg/L range for its 88% product.
# https://www.sigmaaldrich.com/deepweb/assets/sigmaaldrich/product/documents/208/625/l1250pis.pdf
# https://www.laballey.com/products/lactic-acid-88-fcc-grade
# These purchased-product densities do not set modeled makeup-stream density.
acid_densities = {
    'MaleicAcid': None, 'LacticAcid': 1.20,
    'H2SO4': 1.84, 'CitricAcid': None,
}

# Acid purities from supplier websites, assuming weight percentage.
acid_purities = {
    'MaleicAcid': 0.99, 'LacticAcid': 0.88,
    'H2SO4': 0.96, 'CitricAcid': 1.0,
}

# USD/kg purchased product; equal-weight means of TEA.xlsx unit prices.
acid_prices = {
    'MaleicAcid': 93.2465,
    'LacticAcid': 17.13191358140667,
    'H2SO4': 3.9377111829669595,
    'CitricAcid': 6.833934581882406,
}
