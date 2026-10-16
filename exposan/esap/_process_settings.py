"""User-supplied acid and utility prices for preliminary ESAP economics.

Quotes are USD/kg purchased material, with no freight adjustment or cost-year
escalation. Supplier listing dates are unspecified except for the older sulfuric
acid listing (supplier last online April 2024). These are study inputs, not
independently verified current market prices. Citric acid is assumed 100% pure
because the supplied quote specifies anhydrous material without an assay.
"""

__all__ = ('acid_price_quotes', 'acid_prices', 'acid_price_ranges',
           'electricity_price', 'electricity_price_source',
           'natural_gas_price', 'natural_gas_price_per_kmol',
           'natural_gas_price_per_Mcf', 'natural_gas_price_source')

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
_kmol_per_Mcf = _gas_reference_P * (1000 * 0.3048**3) / (
    8314.46261815324 * _gas_reference_T
)
natural_gas_price_per_kmol = natural_gas_price_per_Mcf / _kmol_per_Mcf
natural_gas_price = natural_gas_price_per_kmol / _methane_MW  # USD/kg CH4

acid_price_quotes = {
    'MaleicAcid': dict(
        price_range=(2.20, 2.50), purity=0.99, minimum_order_kg=1000,
        grade='99% powder', supplier='Ningbo Dongbo New Energy',
        source='https://www.made-in-china.com/products-search/hot-china-products/Maleic_Acid_Price.html',
    ),
    'LacticAcid': dict(
        price_range=(1.20, 1.76), purity=0.80, minimum_order_kg=1000,
        grade='80% solution; industrial/cosmetic/food grades advertised',
        supplier='Boqi Chemicals (Shanghai)',
        source='https://www.made-in-china.com/products-search/hot-china-products/Lactic_Acid_80%25.html',
    ),
    'H2SO4': dict(
        price_range=(0.21, 0.21), purity=0.98, minimum_order_kg=20000,
        grade='98%, industrial grade', supplier='Shijiazhuang Xinlongwei Chemical',
        source='https://www.diytrade.com/china/pd/12832305/Industrial_Grade_Sulfuric_Acid_98.html',
    ),
    'CitricAcid': dict(
        price_range=(0.625, 0.660), purity=1.0, minimum_order_kg=25000,
        grade='Anhydrous, food grade', supplier='Dalian Huihai Shengtong Biotechnology',
        source='https://www.made-in-china.com/new-products/citric-acid-anhydrous_7.html',
    ),
}

# Convert purchased-material quotes to USD/kg pure acid for create_system.
acid_price_ranges = {
    ID: tuple(price / quote['purity'] for price in quote['price_range'])
    for ID, quote in acid_price_quotes.items()
}
acid_prices = {
    ID: (low + high) / 2
    for ID, (low, high) in acid_price_ranges.items()
}
