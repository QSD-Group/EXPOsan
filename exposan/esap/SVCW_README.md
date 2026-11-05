# SVCW Treatment System

**Purpose:** Use this workflow to simulate and calibrate the integrated SVCW plant
against site data, assess whether the DW liquid matches the collected SVCW-AD
sample, and distinguish numerical/startup artifacts from genuine process effects.
This guide documents project additions, not underlying model implementations.

**Contents:** [1. Example](#1-simulation-example) |
[2. Layout](#2-layout-and-ownership) |
[3. Functions](#3-function-and-class) |
[4. Inputs](#4-environment-and-inputs) |
[5. Assumptions](#5-key-assumptions-and-conceptual-models) |
[6. Reproducibility](#6-reproducibility-outputs-and-limits)

## 1. Simulation Example

```python
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
from exposan import esap
from exposan.esap.provenance import snapshot_sources
from exposan.esap.analyses.svcw_dw_validation import plot_svcw_dw_validation

stamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%fZ')
output = Path('results') / f'svcw_{stamp}'
system_kwargs = {'ad_initial_s_ip': 24.5}
simulate_kwargs = {
    'method': 'BDF', 't_span': (0, 300), 'rtol': 1e-5, 'atol': 1e-7,
}
snapshot_sources(output / 'sources.zip', settings={
    'system_kwargs': system_kwargs, 'simulate_kwargs': simulate_kwargs,
})
system = esap.create_svcw_system(**system_kwargs)
system.simulate(**simulate_kwargs)
sol = system.scope.sol
assert sol.success, sol.message
assert np.isclose(sol.t[-1], 300), sol.t[-1]
assert np.isfinite(sol.y).all()
system.scope.export().to_pickle(output / 'tracked.pkl')
np.savez_compressed(output / 'solver_states.npz', t=sol.t, y=sol.y)
fig, axes = plot_svcw_dw_validation(
    system, save_as=output / 'dw.jpg', show=False,
)
print(system.flowsheet.unit.DW.outs[0])
```

This uses the current calibration in `systems/SVCW.py`. The initial S_IP is
24.5 mg P/L **before** optional mineral startup projection. SVCW has no ESAP
early-stop option; check the actual endpoint rather than only a flat outlet curve.
For a quick integration check, shorten the interval and endpoint assertion to
0.01 d. Full plant runs can take minutes or longer.

To plot an already-simulated system, call only `plot_svcw_dw_validation(system)`.
The standalone command `python -m exposan.esap.analyses.svcw_dw_validation`
instead constructs a new case using that script's editable argument dictionaries;
its solver tolerances can differ from this example.

## 2. Layout and Ownership

The IDs and links below map directly to the function groups in Section 3.

| ID | ESAP Source | Purpose: When to Use It | Function Group |
| --- | --- | --- | --- |
| S1 | [SVCW.py: construction](systems/SVCW.py) | Used to construct the integrated site model for calibrating plant performance against SVCW data. | [3.1](#31-s1-site-construction-and-components) |
| S2 | [SVCW.py: safeguards](systems/SVCW.py) | Used to prevent startup/numerical artifacts and inconsistent P mapping from misleading site calibration. | [3.2](#32-s2-initialization-and-runtime-safeguards) |
| S3 | [svcw_dw_validation.py](analyses/svcw_dw_validation.py) | Used to assess DW concentration fit against measured SVCW-AD ranges using two comparison panels. | [3.3](#33-s3-sample-comparison-and-plotting) |
| S4 | [provenance.py](provenance.py) | Used to preserve enough source and input evidence to repeat a reported calibration result. | [3.4](#34-s4-shared-source-capture) |
| S5 | [TF project helpers](systems/TricklingFilter_validation.py) | Used to test site TF settings without silently changing the ASR kinetic model or shared feed definitions. | [3.5](#35-s5-shared-tf-project-helpers) |

The ESAP-driven phosphorus-allocation correction in J1/J2 is a shared dependency
in `qsdsan/unit_operations/dynamic/_junction.py`; its ordinary junction machinery
is not re-documented here. Project TF transport/startup changes are likewise
shared dependencies identified in the [TF ownership map](tricklingfilter_README.md#2-layout-and-ownership).
Calibration constants remain in SVCW.py, separate from safeguards and construction.

Project data flow:

```text
RWW + return -> PC -> TF -> ASR -> FC -> SF -> SE
                               ^    |
                               +----+ RAS
PC sludge -> GT --solids--+
FC WAS    -> MT --solids--+-> M1 -> J1 -> AD -> DW -> cake
                                         |
                                         +-> biogas (AD outlet 0)
DW liquid -> J2
GT/MT overflow + J2 liquid + SF backwash -> M2 -> HD -> PC
```

For inspection: M1 outlet 0 is pre-junction AD feed; AD outlet 1 is liquid
digestate; DW outlet 0 is centrate; SF outlet 0 is the final effluent. Tracking is
attached at selected project interfaces, including AD liquid, DW liquid and J2
outlet. DW and its cake now use ADM1p components; centrate ammonium and phosphate
are `S_IN` and `S_IP` in mg N/L and mg P/L. Only the centrate is translated by J2
to mASM2d `S_NH4`/`S_PO4` before returning to M2.

## 3. Function and Class

Read **Used For** to see the calibration or diagnostic task each function serves.
The implementation column explains its role in that task. Private helpers are
normally invoked by system construction, simulation or plotting, not run alone.

### 3.1. S1: Site Construction and Components

Source: S1 in the [layout map](#2-layout-and-ownership).

**Script purpose:** `SVCW.py` is used to construct the complete site model for
simulating plant performance and testing calibration assumptions. Call
`create_svcw_system`, then simulate and compare results as in Section 1. Running
the module's construction entry point alone is not a completed calibration run.

| Function | Used For | Implementation and Rationale |
| --- | --- | --- |
| `create_svcw_system` | Used to prepare the integrated SVCW plant for simulation and calibration against the collected site data. | Applies site inputs, builds independent TF/ASR/AD models, connects the treatment/recycle trains, initializes units and attaches tracking. DW separates digestate in the ADM1p basis before J2 translates only its liquid return. This keeps cake biomass/mineral P out of the instantaneous junction conversion. Keeping these steps in one factory prevents notebook-specific wiring or calibration from becoming an undocumented second system definition. |
| `_create_svcw_components` | Used to test site nutrient-content assumptions without contaminating other calibration cases. | Creates fresh project components and applies the optional XS/XI nitrogen-content multiplier. A system-owned copy prevents a calibration trial from changing the elemental properties of a previously created case. |
| `_synchronize_adm1p_components` | Used to prevent artificial phosphorus gains or losses across the AD interfaces when project component contents are changed. | Copies shared PAO/storage/inert component bases and elemental indices before AD model compilation, returning a change ledger. Synchronization at construction avoids inconsistent phosphorus accounting that would result from changing indices after reactions were compiled. |
| `_create_svcw_wastewater` | Used to test whether a plausible site influent composition explains the measured plant and digestate concentrations. | Merges site composition/COD-fraction overrides with the project recipe and checks positive flow. One explicit feed builder separates measured inputs from calibration estimates instead of distributing adjustments among downstream units. |

### 3.2. S2: Initialization and Runtime Safeguards

Source: S2 in the [layout map](#2-layout-and-ownership).

**Purpose of this code group:** Used to stop inconsistent initial chemistry,
invalid trial-state rate evaluations or altered initial records from misleading
the site calibration. These are internal safeguards, not parameters to adjust
to force agreement with measured concentrations.

| Function or Class | Used For | Implementation and Rationale |
| --- | --- | --- |
| `_equilibrate_ad_mineral_init` | Used to reduce artificial startup phosphate/metal swings caused by inconsistent initial dissolved-ion and mineral pools. | Builds a precipitation-consistent startup state using the selected model's own reaction rows. Nested `saturation_indices` identifies active seeded minerals; sequential projection avoids introducing a second set of elemental stoichiometry, but is not a global equilibrium solver. |
| `_SVCWAnaerobicCSTR.__init__` | Used to apply the site's AD numerical safeguards consistently during an integrated calibration run. | Installs model-local wrappers around both kinetic entry points while retaining the original callables. Local wrapping avoids changing other digesters or evaluating one AD rate path with a different trial-state convention. |
| `_project_kinetic_state` | Used to prevent impossible trial concentrations from generating invalid rates without deleting material from the integrated balance. | Copies the trial state and projects liquid/gas concentrations to nonnegative values for rate evaluation only. Copying is essential: changing the solver's unknowns in place would conceal negatives and alter the balance being integrated. |
| `_projected_flex_rhos`, `_projected_rate_function` | Used to keep both AD kinetic entry points on the same nonnegative trial-state convention. | Apply that same projection before their respective kinetic calls. Both paths are covered so switching the pH/rate entry point does not reintroduce an inconsistent concentration convention. |
| `_update_state` | Used to retain honest raw-state diagnostics while publishing physical outlet values. | Publishes physical outlet values while restoring the raw solver concentrations and retaining updated flow. The hydraulic derivative can then respond to a negative trial instead of silently receiving a clipped, different state. |
| `_SVCWSystem.dynamic_run` | Used to run a plant calibration with sufficient numerical resolution of the TF nutrient/storage states. | Applies the shared TF tolerance guard for the full plant and restores the requested method afterward. This centralizes the numerical setup without requiring users to maintain a plant-length absolute-tolerance vector. |
| `_restore_initial_solution` | Used to ensure that startup plots and initial-inventory calculations refer to the actual initial condition. | Restores the saved initial sample from the actual integration y0 after the run. This preserves the starting inventory in diagnostics instead of allowing subsequent state updates to overwrite the first recorded state. |

`_SVCWAnaerobicCSTR` is used to keep AD kinetic evaluations valid during temporary
negative solver trials. `_SVCWSystem` is used to run the integrated calibration
with the shared TF accuracy guard and trustworthy initial-state records. Both
keep the site safeguards local rather than changing concentration values to fit
the measurements.

### 3.3. S3: Sample Comparison and Plotting

Source: [svcw_dw_validation.py](analyses/svcw_dw_validation.py), S3.

**Script purpose:** Used to assess calibration fit by plotting DW ammonium,
magnesium, calcium and phosphate against measured SVCW-AD concentration ranges.
Use the two-panel comparison to see which species remain outside their measured
ranges and need further investigation; it does not choose parameter updates.

| Function | Used For | Implementation and Rationale |
| --- | --- | --- |
| `measured_concentration_ranges` | Used to compare predicted ammonium and phosphate against measurements on the same mass basis. | Validates the supplied mean/range and converts measured NH4/PO4 ion masses to the requested plotting basis. The supplied +/- range is preserved without assuming a statistical confidence level, preventing a unit conversion from becoming an unrecorded uncertainty model. |
| `_tracked_dw_series` | Used to obtain the recorded DW histories needed to assess calibration fit and residual transients. | Reads the four required DW variables, accepting ADM1p S_IN/S_IP or archived mASM2d S_NH4/S_PO4 on the same elemental basis. Requires unique tracked columns and finite, ordered records. It does not reconstruct a missing history from the final stream, because that would manufacture a time series. |
| `_style_axis` | Used to make sample-comparison figures visually consistent and suitable for reporting calibration results. | Applies the requested font, tick pattern, mirrored axes and limits consistently to both panels. Formatting remains separate from data extraction so figure edits cannot change measured values or simulated concentrations. |
| `plot_svcw_dw_validation` | Used to judge how well DW ammonium, magnesium, calcium and phosphate match the measured SVCW-AD ranges. | Plots NH4 in one panel and Mg/Ca/phosphate in the other with matching measured bands at alpha=0.5. It reads an existing run without simulating again, preserving the association between the displayed trajectory and the user's calibrated case. |

The script entry point is a convenience for a fresh run: it captures sources,
simulates, checks completion and saves timestamped figure/source files. Calling
the plot function alone never substitutes that configuration for an existing case.

### 3.4. S4: Shared Source Capture

Source: [provenance.py](provenance.py), S4; also used by TF and SF workflows.

**Module purpose:** Used to preserve the evidence needed to repeat a reported
calibration: exact source files, settings, inputs and environment identity. It
supports defensible comparisons between calibration iterations, not prediction
of treatment performance.

| Function | Used For | Implementation and Rationale |
| --- | --- | --- |
| `_git_state` | Used to identify whether a calibration depends on local source changes that a commit ID alone would miss. | Records available repository HEAD and dirty/untracked status without requiring a Git checkout. A commit alone is insufficient here because some project changes are untracked; installed-source runs must still be recordable when Git is unavailable. |
| `_source_payloads` | Used to retain the actual code, guides and shared inputs needed to repeat a reported calibration result. | Collects current ESAP code/guides, focused tests and relevant shared source/input files. Exact file bytes are retained because package versions and the small index README cannot identify the actual modified code or the three system guides. |
| `snapshot_sources` | Used to make calibration figures and spreadsheets reproducible from their original code, settings and input files. | Writes a new ZIP with file hashes, explicit arguments, environment versions and optional extra inputs. It refuses overwrite and should run before construction, rather than attaching today's source to a result created earlier with unknown settings. |

The helper uses one shared archive format instead of different provenance
conventions for each validation system. It does not serialize the live system or
replace saving raw and tracked output data.

### 3.5. S5: Shared TF Project Helpers

**Purpose:** These shared helpers are used to test TF calibration within SVCW
without introducing different component or solver conventions from standalone
validation.

| Function | Used For | Implementation and Rationale |
| --- | --- | --- |
| `_create_validation_components` | Used to keep site nutrient-content assumptions consistent with the validation component definitions. | SVCW creates its own set from the shared builder, then applies site adjustments rather than changing another case's components. |
| `_create_validation_model` | Used to test the site's TF kinetic calibration independently of activated-sludge kinetics. | It constructs a separate TF model, so a TF parameter adjustment does not silently alter ASR's model. |
| `_set_tf_init_conc` | Used to initialize the site TF with its own component definitions and meaningful PAO-storage ratios. | SVCW continues to use the configurable E2P-based builder in `TricklingFilter_validation.py`; it does not load D1's stock-component profile or change when that workbook sheet is edited. |
| `_tf_tolerance_solver` | Used to preserve the validated resolution of small TF nutrient/storage states inside the larger plant. | The shared implementation in WERF `systems/D1.py` applies the same rules in D1, SVCW and standalone validation, without applying SVCW feed or component settings to WERF. |

Full implementation details are in
[TF Section 3.1](tricklingfilter_README.md#31-tf1-system-construction-and-runtime).

## 4. Environment and Inputs

Reviewed interpreter: `/opt/anaconda3/envs/algae/bin/python`, Python 3.14.4.
Recorded distribution metadata: QSDsan 1.6.0, EXPOsan 1.5.3, BioSTEAM 2.53.11,
ThermoSTEAM 0.53.5, NumPy 2.4.6, SciPy 1.17.1, pandas 3.0.3 and Matplotlib
3.10.9. The editable-source project versions do not all match installed metadata;
use the source archive and its manifest rather than reinstalling by version
numbers alone. Arial and the Guest palette are needed for the requested styling.

Project inputs are the June 2026 SVCW means, inferred ion/COD-fraction recipe,
selected calibration values in Section 5, and `../werf/data/initial_conditions.xlsx`
(E2P profiles and F1 AD seed). Keep that workbook with the run: the final initial
state includes ESAP's seeding, S_IP override and optional mineral projection.

The figure uses the supplied **SVCW-AD** sample, not the SAF sample:

| Measured Species | Ion/Metal Mass [mg/L] | Element-Basis Mean for Plotting [mg/L] |
| --- | ---: | ---: |
| NH4+ | 1380.06 +/- 137.9 | 1073.38 as N |
| PO4 | 131.1 +/- 1.3 | 42.78 as P |
| Mg | 19.3 +/- 1.6 | 19.3 as Mg |
| Ca | 64.9 +/- 2.1 | 64.9 as Ca |

Ranges are converted with the means. `concentration_basis='element'` uses N/P;
`'ion'` converts simulated concentrations to the original ion-mass basis instead.
DW stream records are read in mg/L; direct AD kinetic-state inspection requires
the correct mass/unit basis. Do not align tracked and solver rows by index alone.

## 5. Key Assumptions and Conceptual Models

### Site Parameters Selected by ESAP

Only the site's configured/calibrated inputs are summarized here, not unchanged
unit parameter sets. Calibration estimates should not be relabeled as raw samples.

| Project Choice | Current Setting | Reason/Evidence |
| --- | --- | --- |
| RWW flow, COD, NH4-N | 13.685 MGD; 600 mg COD/L; 41 mg N/L. | June 2026 plant means used to define the site feed. |
| RWW PO4-P | 2.4 mg P/L. | Calibration assumption; contemporaneous influent orthophosphate was unavailable. |
| RWW Na/K/Mg/Ca/Cl | 278.5/66.5/60/220/643.9 mg/L. | Inferred/calibrated against digestate composition, not measured raw influent. |
| COD allocation | SI=0.14, SF=0.24, SA=0.05, XI=0.235; remainder XS. | Shifts COD from particulate inert to soluble inert to match site TSS without changing assigned total COD/BOD. |
| TF profile | GPS-X flags False; mu_H=6/d, q_PHA=3, q_PP=1.5, b_PP=0.2/d. | Base mASM2d biological rates, also used by D1; TF and ASR remain separate model objects. SVCW component indices remain site-specific. |
| TF setup | Six sections/four layers; 9464 m2 area; recycle 1; dynamic thickness with H0=4e-4 m. | SVCW's physical configuration; D1 retains five layers and its WERF influent and component basis. |
| TF startup/detachment | Shared meaningful PAO seed and capacity detachment scale 0.065. | Avoids artificial guild-startup delay and retains the selected thickness/BOD fit; see TF Section 5. |
| AD initial dissolved P | 24.5 mg P/L before projection. | Replaces an inconsistent startup pool; projected startup may be lower. |
| AD initial dissolved Mg/Ca | 10 mg Mg/L and 80 mg Ca/L after optional projection. | Explicit t=0 concentrations independent of pKsp and feed-ion settings; mineral projection must not silently change them. |
| AD pH / gas-transfer setting | pH fixed at 7.81; kLa=10/d. | Measured pH and configured site gas-transfer assumption; not independently identified by one sample. |
| Apparent pKsp | CaCO3=5.8, struvite=13.0, ACP=25.8, MgCO3=5.9. | Site calibration under the selected chemistry; carbonate estimates used the 270-mM inorganic-C assumption. |
| Mineral rate overrides | CaCO3, newberyite, MgCO3, AlPO4 and FePO4 k_mmp each set to 0.24, in model-specific units. | Reduces competing rate constants relative to the retained struvite/ACP constants; intended to favor the selected phosphate-mineral pathways. |
| SF split | Particulate mass fraction 0.60 and liquid backwash fraction 0.02. | Explicit site solids-capture/return-flow approximation without adding a fouling cycle. |

### Project Safeguards and Governing Checks

**Consistent phosphorus mapping.** The project synchronizes shared component
indices before model compilation and checks `max(abs(nu @ i_P)) <= 1e-10`.
For J1, the allocation correction uses available source P to constrain the
lipid/carbohydrate split while preserving its COD pool. Its conceptual ledger is
`IP_out = TP_in - P_in_other_outlet_pools`; individual dissolved-P concentrations
need not be equal across that component translation. A residual deficit raises
an error instead of being hidden by a zero clamp.

**Dewatering before component translation.** DW uses the ADM1p particulate
flags and mass indices to separate native AD biomass and minerals. It has no
reaction or mineral-equilibrium step: dissolved S_IN, S_IP, Mg and Ca retain
their AD-feed concentrations in both outlets during dynamic simulation. Only
residual solids in the liquid outlet undergo J2's existing instantaneous
conversion. The alternative, converting all digestate before DW, would release
biomass/storage P that should leave in the cake into the liquid return first.
No junction stoichiometry or precipitation kinetics is changed. J2 can still
increase dissolved phosphate downstream of DW because its centrate contains
uncaptured particulates. Cake composition/TSS and recycled nutrient loads are
therefore not expected to reproduce results from the previous topology; retain
the matching source snapshot when comparing runs.

**Startup mineral projection.** ESAP advances the fastest active seeded mineral
along its own reaction row until saturation, then repeats. The selected target
is 0.99 for every pKsp profile to reduce switching at the precipitation
boundary. The projection uses the F1 pools with the selected initial phosphate;
afterward, dissolved Mg and Ca are fixed at 10 and 80 mg/L, respectively. These
are the actual t=0 concentrations even when all four pKsp values match ADM1p
defaults. Imposing them can leave the final startup state supersaturated;
precipitation then follows the normal dynamic model. Disable this project
projection with `equilibrate_ad_minerals=False`; it is skipped with `ad_pH=None`.
This is a startup approximation, not a new precipitation law or a forced mineral
priority during integration. Changing the Ksp-profile switch does not reset the
site's mineral rate overrides.

**Separate raw and kinetic states.** The added wrapper evaluates kinetics at
`C_kin=max(C_trial,0)` on a copy, while retaining the raw concentration for the
balance/hydraulic response. This prevents impossible rate evaluations without
pretending that outlet clipping proves every accepted solver state is positive.
The TF-specific tolerance guard is shared, not independently reimplemented.

**Sample comparison.** The project conversion is NH4-to-N = 14/18 and PO4-to-P
= 31/95, applied to uncertainty bounds as well as means. Guest-colored bands
cover the supplied range with alpha=0.5. They represent a grab-sample comparison,
not an observed 300-day confidence envelope.

No new GPS-X-equivalent AD or SF conceptual model is claimed here. The project's
shared TF conceptual reference is cited in
[TF Section 5](tricklingfilter_README.md#5-key-assumptions-and-conceptual-models).

## 6. Reproducibility, Outputs, and Limits

Save the source ZIP, explicit factory/solver settings, raw `sol.t/sol.y`, tracked
histories, completion status and figure together, as in Section 1. New snapshots
include all three guides. The DW script automatically captures sources and makes
timestamped figure/source files; its plot function alone does not save a full
trajectory or establish provenance for an older system.

To check a saved archive without extracting it:

```python
from hashlib import sha256
import json
from zipfile import ZipFile

archive_path = output / 'sources.zip'  # output from Section 1
with ZipFile(archive_path) as archive:
    manifest = json.loads(archive.read('manifest.json'))
    for name, expected in manifest['sha256'].items():
        assert sha256(archive.read(name)).hexdigest() == expected, name
```

Source capture is not a binary environment or complete dependency-data image;
use `extra_files` for additional inputs and do not edit inputs/code during the
run. Archives can contain private plant data. Earlier immutable archives in
`results/source_snapshots` precede this documentation split and remain historical
records, not proof of the source used for every earlier figure.

```bash
MPLBACKEND=Agg python -m pytest -q \
  tests/test_esap_svcw.py \
  tests/test_esap_provenance.py
```

The [review record](REVIEW_20261015.md) reports the prior 103-test suite and
bitwise-identical short pre/post-cleanup trajectories. It does not establish a
fresh 300-day fit. Remaining project limits are uncertain feed-ion estimates,
fixed-pH calibration, seed-sensitive startup projection, and possible accepted
near-zero negative states. Matching one sample or closing reaction P balances
does not alone establish whole-plant balance or parameter identifiability.
