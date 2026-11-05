
# ESAP minimum selling prices

The MEE `Overall` evaporation specification is water removed divided by feed
water, calculated from the target solute mass fraction and capped at 0.95.

`model.calculate_msp` solves the final product price at zero project NPV and
reports USD/kg final product and USD/kg elemental N in the sold product.
Final product mass includes residual moisture or solution water; nitrogen lost
to mother liquor is excluded from both saleable production denominators.

```python
from exposan.esap._systems import create_system
from exposan.esap.model import calculate_msp

# Illustrative assumptions only; replace with the study's agreed inputs.
tea_assumptions = dict(
    start_year=2026,
    income_tax=0.0,
    annual_labor=0.0,          # USD/year; fill in staffing costs
    annual_maintenance=0.0,    # fraction of fixed capital/year
    system_add_OPEX=0.0,       # USD/year; add disposal, overhead, etc.
    depreciation='SL',
    construction_schedule=(0, 1),
    WC_over_FCI=0.0,
    finance_fraction=0.0,      # all equity
)

for product in (
    'ammonium_sulfate', 'ammonium_citrate',
    'ammonium_lactate', 'ammonium_maleate',
):
    system = create_system(product, simulate=False)
    result = calculate_msp(
        system,
        discount_rate=0.10,   # illustrative target return
        lifetime=20,         # operating years
        operating_hours=8000,
        tea_kwargs=tea_assumptions,
    )
    print(f"{product}: ${result['MSP_USD_per_kg_product']:.3f}/kg product; "
          f"${result['MSP_USD_per_kg_N']:.3f}/kg N; "
          f"NPV residual = ${result['NPV_USD']:.2f}")
```

The function retains feed and coproduct prices, applies the active utility
prices, simulates, and assigns the solved price to the final product. Acid
prices, electricity, and natural gas are set by `create_system`. Natural gas
uses the July 2026 U.S. industrial EIA quote of $4.58/thousand cubic feet,
converted to $/kg assuming ideal methane at 60°F and 14.73 psia. The source
and conversion are retained in `_process_settings.py`. Low-pressure steam is
scaled from BioSTEAM's baseline by the study gas price / 2016 U.S. industrial
gas price (4.58/3.51). Cooling water and chilled water use the electricity ratio
(0.1429/0.0676). Treating the defaults as 2016 prices is an approved preliminary
assumption, not a verified original price year. Temperatures and efficiencies
retain their active settings. Sources, baseline rates, and a TODO for DOE-based
site-specific service costs are in `_process_settings.py`. Other utility prices
remain at their active values. QSDsan's default capital basis is installed equipment
cost; use `CAPEX` in `tea_kwargs` for a separately estimated total investment.
Upstream recovery, labor, maintenance, waste treatment, packaging, freight,
working capital, taxes, and financing need explicit study assumptions before
these MSPs can be treated as project estimates.

Equipment purchase costs use QSDsan's latest bundled annual CEPCI, currently
797.9 for 2023, as a temporary basis rather than a verified 2026 index.
`create_system` sets the shared `qs.CEPCI`/`bst.CE` before constructing units,
including auxiliaries and the custom crystal separator. The bundled table
does not cite its source; replace this basis when a verified study-period
CEPCI is available. Existing systems must be re-simulated to refresh costs.

Acid purchase prices are USD/kg purchased product, with equal-weight arithmetic
means stored in `_process_settings.py`. Supplier package quotes and unit-price
calculations remain in `TEA.xlsx`. Liquid lactic quotes are included; the
60% powder with calcium lactate carrier is excluded. Common purchased-product
purities are maleic 99%, lactic 88%, sulfuric 96%, and citric 100%.
Volume quotes use the documented purchased-solution densities (lactic 1.20 kg/L,
sulfuric 1.84 kg/L); mass quotes need no density. These densities do not set the
modeled makeup-stream density. Spreadsheet edits require updating the code's
mean prices; the workbook is not read at runtime.

`create_system(acid_prices=..., acid_purities=...)` charges required acid mass
divided by supplier purity, multiplied by the purchased-product unit price.
The resulting stream price is this cost divided by total modeled stream mass.
Custom `acid_prices` now use the purchased-product basis; to supply pure-acid
prices, explicitly set the corresponding `acid_purities` entries to 1.

# ESAP Guides

Choose the system you are working with:

- [Trickling-filter validation](tricklingfilter_README.md)
- [SVCW treatment system](SVCW_README.md)
- [Sand-filter validation](sandfilter_README.md)

Each guide follows the same order: simulation example; layout and ownership;
functions and classes; environment and inputs; key assumptions and conceptual
models; reproducibility, outputs, and limits. Documentation is scoped to
ESAP-created code, project modifications, and calibration decisions.

## Tests

- [test_esap_sand_filter.py](../../tests/test_esap_sand_filter.py): feed definitions,
  component splitting and dynamic backwash/recycle behavior.
- [test_esap_svcw.py](../../tests/test_esap_svcw.py): plant construction, calibration,
  AD/DW consistency and measured-concentration plotting.
- [test_esap_trickling_filter.py](../../tests/test_esap_trickling_filter.py): TF
  physics, D1/F1 input compatibility, validation profiles, solver safeguards,
  PAO startup and sweep exports.
- [test_esap_provenance.py](../../tests/test_esap_provenance.py): source/input
  archives, hashes, metadata and overwrite protection.

Run from the EXPOsan repository root:

```bash
MPLBACKEND=Agg python -m pytest -q \
  tests/test_esap_sand_filter.py tests/test_esap_svcw.py \
  tests/test_esap_trickling_filter.py tests/test_esap_provenance.py
```

[Review record](REVIEW_20261015.md) and
[historical source-snapshot index](results/source_snapshots/SOURCE_INDEX.json).

