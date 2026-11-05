# Sand-Filter Validation

**Purpose:** Use this workflow to test how solids-capture and backwash settings
affect effluent quality and backwash loads before coupling the filter to SVCW.
It provides a controlled system for checking candidate filtration settings and
mass partitioning. For the filter in the full plant, use the [SVCW guide](SVCW_README.md).

**Contents:** [1. Example](#1-simulation-example) |
[2. Layout](#2-layout-and-ownership) |
[3. Functions](#3-function-and-class) |
[4. Inputs](#4-environment-and-inputs) |
[5. Assumptions](#5-key-assumptions-and-conceptual-models) |
[6. Reproducibility](#6-reproducibility-outputs-and-limits)

## 1. Simulation Example

### Constant-Splitting Comparison

```python
from exposan import esap

system = esap.create_sand_filter_validation_system(
    wastewater_strength='mid', model='constant_splitting',
    sf_kwargs={
        'isdynamic': False,
        'backwash_solids_mass_fraction': 0.60,
        'backwash_flow_fraction': 0.02,
    },
)
system.simulate()
SF = system.flowsheet.unit.SF
for label, stream in [('Feed', SF.ins[0]),
                      ('Effluent', SF.outs[0]), ('Backwash', SF.outs[1])]:
    print(label, {
        'COD': stream.COD, 'BOD': stream.BOD, 'TSS': stream.get_TSS(),
        'TN': stream.TN, 'PO4-P': float(stream.iconc['S_PO4']),
        'TP': stream.TP,
    })  # Concentrations in mg/L of the indicated basis.
```

The two fractions deliberately select the SVCW splitting recipe for this example;
the validation factory does not introduce an independent fitted SF parameter set.
`isdynamic=False` makes this a static comparison without a time-history claim.

### Empirical-Fouling Topology Check

```python
from exposan import esap

dynamic_system = esap.create_sand_filter_validation_system(
    wastewater_strength='low', model='empirical_fouling',
)
dynamic_system.simulate(method='BDF', t_span=(0, 0.1), rtol=1e-5, atol=1e-8)
sol = dynamic_system.scope.sol
assert sol.success, sol.message
SF_dynamic = dynamic_system.flowsheet.unit.SF
print('Captured-solids inventory:', SF_dynamic.captured_solids_inventory)
```

The short second example tests construction/integration only, not long-term
fouling calibration. Dynamic tracking is attached by the factory for this mode.

## 2. Layout and Ownership

The IDs below are repeated in Section 3, linking source responsibilities to the
functions that implement them.

| ID | ESAP Source | Purpose: When to Use It | Function Group |
| --- | --- | --- | --- |
| SF1 | [SandFilter_validation.py](systems/SandFilter_validation.py) | Used to test how filter settings affect effluent quality, nutrient partitioning and backwash loads. | [3.1](#31-sf1-validation-factory) |
| SF2 | [TF feed helpers](systems/TricklingFilter_validation.py) | Used to keep influent assumptions identical when comparing SF and TF validation results. | [3.2](#32-sf2-shared-project-inputs) |
| SF3 | [provenance.py](provenance.py) | Used to preserve a filtration test's code, inputs and settings for later reproduction. | [3.3](#33-sf3-shared-source-capture) |

The factory uses the filtration unit without duplicating its equations in ESAP.
This guide describes the project's mode selection and connections, not unrelated
unit-model internals. There is no SF-specific heatmap or plotting script in the
current ESAP package. SVCW's constant-split selection is configured separately in
[SVCW.py](systems/SVCW.py).

## 3. Function and Class

The descriptions distinguish what a function is **used for** from how it is
implemented. The public factory is the normal starting point; shared private
helpers support it rather than forming separate calibration scripts.

### 3.1. SF1: Validation Factory

Source: SF1 in the [layout map](#2-layout-and-ownership).

**Script purpose:** `SandFilter_validation.py` is used to evaluate candidate
filter settings independently of upstream biology and plant recycle. This helps
you decide whether a change in solids capture or backwash explains an effluent
COD, TSS, TN or TP result before making that change in the full SVCW system.

**`create_sand_filter_validation_system` is used to create the controlled test
case for that comparison.** It builds a fresh flowsheet, validates the
requested wastewater profile, creates project-owned components/feed, and wires
the selected filtration mode. It copies `sf_kwargs` and permits a `model` entry
there to override the top-level mode, so configuration dictionaries can be reused
without in-place mutation or conflicting hidden settings.

For constant splitting, the project connects one raw-water inlet to two named
outlets, `effluent` and `backwash`. For empirical fouling, it additionally creates
`filtered_water_recycle` and connects that outlet to the backwash inlet. This
explicit recycle is necessary to close the selected mode's water path; treating
the third outlet as discharge would construct a different test system.

The factory attaches stream/unit tracking only when the unit is dynamic. Static
comparisons therefore do not pretend to provide a simulated time history. One
factory with explicit topology branches was retained rather than separate copied
factories whose wastewater definitions could drift during later calibration.
No new ESAP filtration class or automatic parameter-fitting routine is defined
here. Run the factory from a script/notebook, vary the selected settings, and
compare the resulting streams against the intended performance or available data.

### 3.2. SF2: Shared Project Inputs

Source: SF2 in the [layout map](#2-layout-and-ownership). The shared definitions
are described in more detail in [TF Section 3.1](tricklingfilter_README.md#31-tf1-system-construction-and-runtime).

**Purpose of this code group:** Used to hold the influent and nutrient-content
assumptions fixed while comparing different filtration settings or TF/SF results.
Otherwise a feed change could be mistaken for an improvement in filter performance.

| Helper or Alias | Used For | Implementation and Rationale |
| --- | --- | --- |
| `_normalize_wastewater_strength` | Used to select the same validation strength for SF and TF comparisons. | Uses the same accepted profile names and early input checks as TF validation. A shared validator prevents an apparently identical SF case from quietly selecting a different water composition. |
| `_create_validation_components` | Used to hold nutrient-content assumptions fixed while testing filtration settings. | Creates fresh components and applies the project content indices. It isolates concurrent or sequential validation cases instead of changing a globally shared component table. |
| `_create_validation_wastewater` | Used to separate the effect of filter settings from the effect of changing the influent. | Builds the common strength-dependent feed. Reusing the builder means differences in TF/SF results can be attributed to the units/settings rather than duplicated influent recipes. |
| `_WASTEWATER_STRENGTHS`, `_WASTEWATER_ORGANIC_FRACTIONS`, `_NONDEFAULT_COMPONENT_PROPERTIES` | Used to keep SF and TF calibration inputs synchronized when shared feed definitions are updated. | Retained aliases reference the TF definitions directly for project callers/tests. They are intentionally shared objects, not independent SF calibration tables that need to be updated twice. |

### 3.3. SF3: Shared Source Capture

**Module purpose:** `provenance.py` is used to preserve the settings and source
behind a particular filtration-performance comparison, so it can be repeated.
Use `snapshot_sources` before creating the test case. Its `_source_payloads` and
`_git_state` helpers collect the current project bytes
and repository identity; their rationale is documented under
[SVCW Section 3.4](SVCW_README.md#34-s4-shared-source-capture). SF does not define
a competing archive format or capture results automatically merely on construction.

## 4. Environment and Inputs

Reviewed environment: Python 3.14.4; QSDsan/EXPOsan installed metadata 1.6.0/1.5.3,
NumPy 2.4.6 and SciPy 1.17.1. The local editable source includes the project's
dependencies; package versions alone are not a sufficient source identifier.
Use a source manifest for the full environment record.

For low/mid/high profiles the project uses 10 MGD and 20 degC:

| Strength | COD [mg/L] | NH4-N [mg/L] | PO4-P [mg/L] |
| --- | ---: | ---: | ---: |
| low | 339 | 14 | 1.6 |
| mid | 508 | 20 | 2.4 |
| high | 1016 | 41 | 4.7 |

COD fractions and component N/P contents are the shared project settings in
[TF Sections 4-5](tricklingfilter_README.md#4-environment-and-inputs), not an
additional SF fit. `'default'` is a separate delegated feed profile. The caller
supplies any case-specific split/fouling/energy overrides through `sf_kwargs`.
No initial-condition workbook is applied by the SF factory, and no GPS-X
comparison workbook is required for these examples.

## 5. Key Assumptions and Conceptual Models

**Controlled unit comparison.** The ESAP system is `RWW -> SF -> effluent/backwash`
without upstream biological treatment. Identical project feeds are deliberately
used for TF/SF validation; the SF recipe is not fitted by altering influent P/N
contents independently from TF.

**Explicit water/solids destinations.** The example selects particulate capture
0.60 and liquid backwash fraction 0.02 to match SVCW's project configuration.
The governing project check is component-wise
`mass_feed = mass_effluent + mass_backwash` for constant splitting. These fractions
describe a selected separation/return-flow approximation, not biological nutrient
removal or a measured fouling trajectory.

**Mode-specific topology.** Empirical fouling is an alternative test configuration
with a filtered-water recycle; it is not substituted into SVCW's two-outlet
connection. No additional SF kinetic coefficients were calibrated in this wrapper.
No equivalence to a particular GPS-X filtration conceptual model is asserted.

## 6. Reproducibility, Outputs, and Limits

Record `wastewater_strength`, `model`, every `sf_kwargs` override and any dynamic
solver arguments. Use `snapshot_sources(..., settings=...)` before construction;
the common example/archive checks are in
[SVCW Sections 1 and 6](SVCW_README.md#1-simulation-example).
Static mode provides final streams, not a scope time series. For dynamic mode,
save tracked histories and raw solver arrays as well as the completed endpoint.

```bash
MPLBACKEND=Agg python -m pytest -q tests/test_esap_sand_filter.py
```

The tests protect feed-definition identity, constant-split component conservation
and dynamic recycle wiring/short integration. They do not establish a fitted
long-term fouling trajectory or validate particulate capture against samples.
Likewise, matching a phosphorus concentration after separation does not establish
chemical precipitation or biological P removal in this standalone test.
