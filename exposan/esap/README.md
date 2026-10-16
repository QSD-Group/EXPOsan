# ESAP minimum selling prices

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
and conversion are retained in `_process_settings.py`. Other utility prices remain
at their active values. QSDsan's default capital basis is installed equipment
cost; use `CAPEX` in `tea_kwargs` for a separately estimated total investment.
Upstream recovery, labor, maintenance, waste treatment, packaging, freight,
working capital, taxes, and financing need explicit study assumptions before
these MSPs can be treated as project estimates.
