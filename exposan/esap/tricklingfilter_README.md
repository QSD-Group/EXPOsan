# Trickling-Filter Validation

**Purpose:** Use this workflow to calibrate and validate standalone TF settings
against GPS-X results and visualize the remaining discrepancies with heatmaps.
For the TF inside the plant, use the [SVCW guide](SVCW_README.md).

**Contents:** [1. Example](#1-simulation-example) |
[2. Layout](#2-layout-and-ownership) |
[3. Functions](#3-function-and-class) |
[4. Inputs](#4-environment-and-inputs) |
[5. Assumptions](#5-key-assumptions-and-conceptual-models) |
[6. Reproducibility](#6-reproducibility-outputs-and-limits)

## 1. Simulation Example

### One Strength and Layer Count

```python
from exposan import esap

system = esap.create_trickling_filter_validation_system(
    wastewater_strength='low',
    gpsx_kinetics=False, gpsx_decay=False, gpsx_aob_proxy=False,
    asm_kwargs={'mu_H': 3.3, 'q_PHA': 6.0, 'q_PP': 2.0, 'b_PP': 0.05},
    tf_kwargs={
        'N_section': 6, 'N_biofilm_layer': 5,
        'filter_bed_surface': 9464, 'recirculation_ratio': 1,
        'dynamic_biofilm_thickness': True,
        'initial_biofilm_thickness': 4e-4,
    },
)
system.simulate(method='BDF', t_span=(0, 300), rtol=1e-3, atol=1e-2,
                stop_at_steady_state=False)
sol = system.scope.sol
assert sol.success, sol.message
assert abs(sol.t[-1] - 300) < 1e-8
TF = system.flowsheet.unit.TF
effluent = TF.outs[0]
print({
    'DO': float(effluent.iconc['S_O2']), 'BOD': effluent.BOD,
    'COD': effluent.COD, 'NH4-N': float(effluent.iconc['S_NH4']),
    'NOx-N': float(effluent.iconc['S_NO3']), 'TN': effluent.TN,
    'PO4-P': float(effluent.iconc['S_PO4']), 'TP': effluent.TP,
})  # Concentrations in mg/L of the indicated basis.
print('Section thicknesses [m]:', TF.total_biofilm_thickness)
print('Power [kW]:', system.get_energy_summary())
```

For a quick integration check, change the end time and endpoint assertion to
0.01 d. Six **sections** are locations down the filter; five **layers** resolve
biofilm depth within each section. Thickness output is one total per section.

### Nine-Case Workbook and Heatmap

```python
from exposan.esap.analyses.tf_validation_sweep import run_validation_sweep

run = run_validation_sweep(
    template_path='/absolute/path/TF_specific_GPS-X_validation.xlsx',
    output_dir='results/tf_comparison', show_plot=False,
    system_kwargs={'gpsx_kinetics': False, 'gpsx_decay': False,
                   'gpsx_aob_proxy': False},
)
print(run.workbook_path, run.heatmap_path, run.source_path)
```

Replace the template path with the actual workbook. The sweep tests low/mid/high
strength at 3/4/5 layers with six sections. Its library defaults use the four
ASM overrides in the first example. The editable launcher is
`python -m exposan.esap.analyses.tf_validation_script`.

**Profile distinction:** the launcher enables all three GPS-X flags and retains
`mu_H=3.3, q_PHA=6, q_PP=2, b_PP=0.05`; the library sweep disables the flags.
Explicit `asm_kwargs` win over flag values. Nested dictionaries merge; use
`asm_kwargs=None` to clear inherited overrides.

## 2. Layout and Ownership

The IDs below are repeated in Section 3. Function-group links connect each
file's responsibility directly to its implementation rationale.

| ID | ESAP Source | Purpose: When to Use It | Function Group |
| --- | --- | --- | --- |
| TF1 | [TricklingFilter_validation.py](systems/TricklingFilter_validation.py) | Used to generate standalone TF performance results for calibrating or validating one chosen configuration. | [3.1](#31-tf1-system-construction-and-runtime) |
| TF2 | [models.py](models.py) | Used to compare electricity demand and cost when treatment settings such as recirculation are changed. | [3.2](#32-tf2-energy-accounting) |
| TF3 | [tf_validation_sweep.py](analyses/tf_validation_sweep.py) | Used to evaluate a candidate parameter set across nine cases and save the calibration comparison. | [3.3](#33-tf3-sweep-and-workbook-export) |
| TF4 | [tf_validation_heatmap.py](analyses/tf_validation_heatmap.py) | Used to show where simulated outputs disagree with GPS-X through a comparison heatmap. | [3.4](#34-tf4-heatmap) |
| TF5 | [tf_validation_script.py](analyses/tf_validation_script.py) | Used for manual parameter calibration with heatmap illustrations comparing QSDsan and GPS-X results. | [3.5](#35-tf5-launcher) |
| TF6 | [provenance.py](provenance.py) | Used to preserve the code and inputs needed to repeat a calibration iteration. | [3.6](#36-tf6-shared-source-capture) |

ESAP-driven TF transport/startup changes reside in the shared
[TF implementation](../werf/_units.py) and
[solver helpers](../werf/systems/D1.py). The configurable startup builder lives
in [TF1](systems/TricklingFilter_validation.py) and is also used by SVCW.
D1/D1E instead read the explicit section/layer profile from workbook sheet D1.
Section 5 records the project changes only; existing biological-model internals
are outside this guide.

D1 now uses the same base biological rate profile as the current SVCW TF:
`mu_H=6`, `q_PHA=3`, `q_PP=1.5`, `b_PP=0.2`, with no additional anaerobic
heterotrophic-decay term. It deliberately retains F1's standard WERF influent and
stock component definitions, rather than SVCW's site-specific inputs. D1's
physical defaults, including five layers, are unchanged; `tf_kwargs` can override
them. D1, SVCW and standalone validation share the targeted tolerance helpers.

## 3. Function and Class

Read **Used For** to decide why a function is needed; read **Implementation and
Rationale** for how it supports that task. Private helpers normally run through
the public factories or analysis scripts, rather than being called individually.

### 3.1. TF1: System Construction and Runtime

Source: [TF1 in the layout](#2-layout-and-ownership).

**Script purpose:** Used to generate effluent, biofilm-thickness and power results
for testing one TF parameter set. Use the system factory for a single-case
calibration or stability investigation; use the launcher in Section 3.5 for the
standard nine-case heatmap comparison. Creating the system alone does not perform
a parameter fit or run its simulation.

#### Construction and Shared Inputs

| Function | Used For | Implementation and Rationale |
| --- | --- | --- |
| `create_trickling_filter_validation_system` | Used to prepare one TF case for calibrating or validating effluent quality, biofilm thickness and power demand independently of the full plant. | Builds a fresh `RWW -> TF -> effluent` case, applies project settings/startup, and attaches tracking. A standalone factory excludes plant recycle and upstream-treatment effects when comparing the TF itself. |
| `create_trickling_filter_validation_model` | Used to evaluate electricity and operating-cost tradeoffs for an already simulated TF parameter set. | Adds the project energy metrics to an existing or new case without integrating again. Creation arguments are rejected with an existing system so metric evaluation cannot silently replace a calibrated simulation. |
| `_normalize_wastewater_strength` | Used to ensure that a requested strength selects a defined, comparable validation case. | Validates case-insensitive default/low/mid/high names before construction. Sharing this check with SF prevents the same strength label from selecting different feeds across validation workflows. |
| `_create_validation_components` | Used to keep nutrient-content assumptions consistent while testing multiple calibration cases independently. | Creates fresh components, applies project content indices, and refreshes constants. Independent copies prevent later cases from modifying the elemental definitions of an earlier simulation. |
| `_create_validation_wastewater` | Used to supply the controlled influent needed to attribute output differences to TF settings rather than different feed recipes. | Builds the selected project feed at 10 MGD and 20 degC, or delegates the explicitly requested default profile. Shared construction avoids separate TF/SF tables drifting during calibration. |
| `_create_validation_model` | Used to test whether changing biological coefficients or the selected kinetic profile improves the TF comparison. | Applies selectable project kinetic profiles and then caller overrides. Nested `rates_with_anaerobic_decay` copies the rate buffer and adds only the declared heterotrophic term; a full alternative biological model is not introduced. |
| `_set_tf_init_conc` | Used to initialize SVCW and standalone validation while changing component definitions, dry density or PAO-storage seeds. | Retains the E2P-based gradients, meaningful storage ratios and dry-solids normalization. It lives here rather than in D1 because D1/D1E now read their fixed startup profiles from the workbook. |

#### Numerical Safeguards and Optional Stopping

`_ValidationBDF` and `_tf_tolerance_solver` are implemented in
[D1.py](../werf/systems/D1.py). The validation
module imports them under their existing names for compatibility. This keeps WERF
independent from ESAP while avoiding separate numerical implementations.
Default BDF runs use a numerical Jacobian without a custom sparsity pattern.
The perturbation guard and targeted tolerances do not depend on such a pattern.

| Class or Method | Used For | Implementation and Rationale |
| --- | --- | --- |
| `_create_d1_tf_model` | Used to apply the current SVCW TF biological rates in the WERF D1 scenario without importing site-specific feed assumptions. | Constructs base electron-acceptor-dependent mASM2d with explicit `mu_H`, `q_PHA`, `q_PP` and `b_PP`; removes the former GPS-X overrides and additional anaerobic decay. It uses the already-created stock WERF components. |
| `create_d1_system` | Used to evaluate the D1 treatment train with SVCW-aligned TF biology while keeping comparisons with F1 on the same influent and component basis. | Retains `default_rww`, stock components, physical defaults and TF overrides. Reads the explicit D1 profile and AD row from `initial_conditions.xlsx`, and returns `_D1System` for targeted accuracy. |
| `TricklingFilter.set_init_conc` | Used to load explicit liquid/biofilm startup profiles without regenerating biological seeds. | A two-level section/layer DataFrame loads directly on an identical grid or uses conservative overlap-weighted averages on a different uniform grid. It does not renormalize workbook concentrations; array and single-index inputs retain their existing behavior. |
| `_D1System.dynamic_run` | Used to apply the same TF numerical safeguards inside D1 as in SVCW. | Wraps the requested solver after state assembly and restores the caller's solver method even if integration fails. Only the TF state block is tightened; other units, `rtol` and `max_step` retain caller settings. |
| `TricklingFilter.get_retained_mass` | Used to include attached and suspended TF biomass in the existing WERF SRT report. | Sums selected component concentrations times their dry-mass indices and liquid/biofilm volumes, returning g dry solids. Biofilm volumes retain the ODE's thickness bounds and numerical floors; a reporting method is used instead of changing the SRT convention or any kinetic equation. |
| `_ValidationBDF._validate_jac` | Used to prevent extreme numerical trial states from disrupting a TF calibration run. | Adds nested `bounded_jacobian` to guard persistent perturbations of inactive states. Bounding numerical probes avoids extreme trial concentrations without changing reaction equations or weakening requested accuracy. |
| `_tf_tolerance_solver` | Used to improve numerical resolution of small nutrient/storage concentrations and biofilm thickness during parameter testing. | Creates local `TFToleranceSolver` after the system state layout exists. Its `__init__` tightens selected concentration/thickness tolerances and copies y0, avoiding user-assembled tolerance vectors and initial-state aliasing. |
| `_SteadyStateWindow.__init__` | Used to define what local steady state means before optionally shortening a validation run. | Validates the project window, relative criterion and scale floors, then allocates history. Stopping settings remain separate from integration tolerances because apparent stationarity and numerical accuracy are different tests. |
| `_SteadyStateWindow.update` | Used to decide whether a simulated TF has settled over a complete window rather than merely returning to the same endpoint. | Checks the full sampled range over a complete window with an interpolated start. Endpoint equality alone was rejected because an oscillation can return to its initial value while remaining unstable between endpoints. |
| `TricklingFilterValidationSystem.dynamic_run` | Used to obtain a guarded TF trajectory or an optional local steady-state result for comparison with reference performance. | Combines guarded integration, optional stopping and `steady_state_report`. Nested `monitored_state` selects states/inventories; `SteadySolver.__init__` and `step` monitor accepted states and publish the actual endpoint rather than a trial evaluation. |

`TricklingFilterValidationSystem` is used to manage a TF validation run and keep
its numerical controls and treatment/energy summaries together. It keeps reporting
with the case that produced the data, rather than requiring each notebook to
rebuild its stream relationships.

| Method | Used For | Implementation and Rationale |
| --- | --- | --- |
| `nitrogen_summary` | Used to inspect where nitrogen remains in the feed and effluent when calibrating nitrogen removal. | Returns inlet/outlet NH4-N, NOx-N, TN, other/organic N and particulate N. Explicit NOx labeling avoids presenting the lumped state as a separate nitrate/nitrite prediction. |
| `get_influent_pumping_energy`, `get_recirculation_pumping_energy`, `get_rotating_arm_energy` | Used to compare the raw-feed, recycle and distributor contributions to the power demand of candidate TF settings. | Read three project power terms separately from the current case. This makes recycle/head assumptions auditable; the historical names return kW, not time-integrated energy. |
| `get_pumping_energy` | Used to compare total pumping demand without conflating it with distributor power. | Sums the two pump powers without the distributor. An explicit subtotal avoids ambiguity about which devices are included in a reported pumping demand. |
| `get_energy_summary` | Used to report a consistent energy breakdown alongside treatment performance when choosing a TF configuration. | Returns all three named powers in one mapping. Workbook export can convert the same quantities to daily electricity without defining a competing accounting convention. |

### 3.2. TF2: Energy Accounting

Source: [models.py](models.py), corresponding to TF2 in Section 2.

**Module purpose:** Used to compare the electricity and operating-cost consequences
of candidate TF settings alongside their treatment performance. This support
module is called by the model/reporting workflow, not run as a calibration script.

| Function | Used For | Implementation and Rationale |
| --- | --- | --- |
| `_validate_efficiency` | Used to reject pump-efficiency inputs that would invalidate an energy comparison. | Requires a finite hydraulic efficiency in (0,1]. Early rejection prevents singular or physically invalid power estimates from entering an otherwise valid comparison workbook. |
| `get_influent_pumping_power` | Used to estimate the energy consequence of the chosen raw-water lift and flow assumptions. | Estimates raw-water lift power using project head/efficiency assumptions. Raw feed rather than recycle-inclusive flow is used so the inlet pump is not counted twice. |
| `get_trickling_filter_recirculation_pumping_power` | Used to evaluate the pumping penalty of increasing TF recirculation while calibrating treatment performance. | Estimates recycle demand from ratio and bed/distribution head, returning zero without recycle. A separate term isolates the energy penalty of changing this validation input. |
| `get_trickling_filter_rotating_arm_power` | Used to include distributor demand when comparing total electricity among TF configurations. | Applies the selected media-volume-specific distributor demand. Keeping this empirical estimate separate avoids hiding it in pump efficiency or total electricity. |
| `add_trickling_filter_energy_metrics` | Used to evaluate treatment-performance and electricity-cost tradeoffs together through the model metrics. | Registers power and optional daily cost. Nested `get_influent_pump_power`, `get_recirculation_pump_power`, `get_tf_rotating_arm_power`, and `get_tf_energy_cost` are callbacks reading current results, not extra solves or unused functions. |

### 3.3. TF3: Sweep and Workbook Export

Source: [tf_validation_sweep.py](analyses/tf_validation_sweep.py), TF3.

**Script purpose:** Used to evaluate one candidate calibration consistently across
low/mid/high strengths and 3/4/5 layers, then save the comparison workbook and
heatmap. It removes the manual work of constructing and reporting nine separate
cases; it does not search for the best parameter values automatically.

| Function or Class | Used For | Implementation and Rationale |
| --- | --- | --- |
| `ValidationRun` | Used to keep the workbook, heatmap, numerical results and source record associated with the same calibration iteration. | Returns paths, case metrics, variation, figure and source archive together. Associating data and artifacts avoids ambiguous notebook variables left over from different runs. |
| `_json_value` | Used to make a tested parameter set recordable and reusable in calibration metadata. | Converts NumPy values/paths to recordable settings and rejects opaque objects. A model object's string representation is not accepted as sufficient information to reproduce its configuration. |
| `_prepare_kwargs` | Used to apply one candidate parameter set consistently across all cases represented in the heatmap. | Deep-copies/merges arguments and enforces sweep-owned dimensions, fresh flowsheets and six sections. These checks prevent a complete-looking workbook whose cases do not match their labels. |
| `_TAG` | Internal helper used to address the correct spreadsheet elements when exporting calibration results. | Constructs qualified spreadsheet XML names consistently. This avoids relying on whichever namespace prefix a particular workbook uses. |
| `_worksheet_parts` | Internal helper used to find the correct strength sheet in the comparison workbook. | Resolves sheet names through workbook relationships instead of assuming internal filenames. Templates can therefore retain their own sheet numbering and layout. |
| `_cell_value` | Internal helper used to recover the reference values needed for the calibration comparison. | Reads typed cells and cached formula results. It does not evaluate formulas, so an already calculated reference workbook is required instead of attempting unreliable partial spreadsheet recalculation. |
| `_load_template` | Used to check that a reference workbook can support the planned comparison before running nine simulations. | Validates sheets, headers, reference cells and case rows before simulation. Early validation avoids spending nine integrations on an incompatible comparison template. |
| `_case_results` | Used to collect the same performance indicators for every candidate TF case. | Extracts the agreed concentrations, six thicknesses and daily electricity, rejecting invalid values. One extraction path keeps printed summaries and exported metrics consistent. |
| `_setting_rows` | Used to identify exactly which settings and solver endpoints produced a calibration workbook and heatmap. | Records effective arguments, endpoints and provenance; nested `add` flattens dictionaries. Separate readable parameter rows are preferred to a single large, difficult-to-audit serialized cell. |
| `_xml_bytes` | Internal helper used to keep exported comparison workbooks readable by Excel. | Preserves compatibility namespace declarations during XML serialization. Removing apparently unused declarations can invalidate Excel-specific content that the scientific analysis does not inspect. |
| `_set_cell` | Internal helper used to populate simulated-result cells without altering the reference presentation. | Replaces a typed cell value while retaining style attributes. Inserting simulated results therefore does not discard the reference template's presentation. |
| `_settings_xml` | Used to publish the calibration recipe inside the same workbook as its results. | Builds the wrapped, filtered settings sheet with a frozen header. The sheet is an attached record of the run, not another independent source of defaults. |
| `_save_workbook` | Used to retain a calibration iteration's results while preserving the fixed GPS-X benchmark. | Updates result/settings parts without rewriting reference rows or formula caches. Nested `result_style` adds formats without changing reference styles; whole-workbook rewriting was avoided to protect cached reference values. |
| `run_validation_sweep` | Used to evaluate one chosen parameter set across all nine strength/layer cases and produce a calibration workbook and heatmap. | Captures sources, runs nine fresh cases, checks endpoints, writes results, then plots saved values. A failed case does not yield a misleading complete workbook, although its initial source archive may remain. |

### 3.4. TF4: Heatmap

Source: [tf_validation_heatmap.py](analyses/tf_validation_heatmap.py), TF4.

**Script purpose:** Used to visualize which predicted concentrations, thicknesses
and energy demands disagree with GPS-X, so you can decide what to calibrate next.
Use it directly when you already have a result workbook and only need to inspect
or regenerate its comparison figure.

| Function | Used For | Implementation and Rationale |
| --- | --- | --- |
| `_metric_values` | Used to ensure that each heatmap cell compares matching treatment or electricity indicators. | Selects 15 metrics and sums NOx/electricity before comparison. This keeps both sources on the same reporting basis instead of comparing mismatched individual columns. |
| `calculate_percent_variation` | Used to quantify which outputs overpredict or underpredict the GPS-X reference. | Calculates the bounded signed max-normalized difference, including zero/zero handling. This was selected for the diverging heatmap and is not conventional reference-denominator percent error. |
| `compile_percent_variation` | Used to assemble all strength/layer discrepancies into the matrix needed for a calibration heatmap. | Reads saved sheets and orders the nine comparisons consistently. The workbook, not a separate in-memory array, is the figure's source of truth. |
| `_format_percent` | Used to make small and large model discrepancies readable in heatmap annotations. | Selects compact annotation precision without rounding result data. The figure remains readable while full precision stays available in the workbook. |
| `heatmap` | Used to visually locate the outputs, strengths and layer counts that need further calibration. | Draws Guest colors, labels, annotations and colorbar, then optionally exports the figure. Separating rendering from calculation prevents a style adjustment from changing the comparison metric. |
| `plot_validation_heatmap` | Used to inspect or redraw a saved calibration comparison without rerunning the nine simulations. | Connects workbook reading and rendering, returning both data and figure. A heatmap can be regenerated without rerunning the process model. |

### 3.5. TF5: Launcher

**Purpose:** [tf_validation_script.py](analyses/tf_validation_script.py) is used
for **manual parameter calibration with heatmap illustrations comparing QSDsan
and GPS-X results**. This is the main file to edit when testing a different TF
parameter set through the standard nine-case workflow.

1. Change `SYSTEM_KWARGS`: biological coefficients belong in `asm_kwargs`;
   physical/transport settings such as `detachment_rate` belong in `tf_kwargs`.
2. Run the script to simulate the nine cases and generate their workbook and
   comparison heatmap.
3. Inspect the discrepancies, choose the next parameter adjustment, and repeat.

The script delegates simulation/export to `run_validation_sweep`, so every
iteration uses the same comparison procedure. It automates the runs and heatmap
generation, not parameter selection or numerical optimization.

The launcher's `SIMULATE_KWARGS` uses BDF with `rtol=1e-5`, `atol=1e-6`,
`max_step=0.25` d, and no early steady-state termination. This tighter numerical
profile supports calibration without changing kinetics. The existing
`_tf_tolerance_solver` automatically limits `S_O2`, `S_NH4`, `S_NO3`, `S_PO4`,
`X_PHA`, `X_PP`, and `X_PAO` absolute tolerances to `1e-6` mg/L in both liquid
and biofilm states, and thickness tolerances to `1e-9` m. Tighter caller values
are preserved. No AD-specific concentration or gas tolerances are applied:
the standalone TF has neither state block. The sweep library's own defaults
remain unchanged; this stricter profile is selected by the launcher and recorded
in the output workbook and source archive. Greater numerical accuracy does not
guarantee that biological transients disappear or that 300 days is steady state.

### 3.6. TF6: Shared Source Capture

**Purpose:** `provenance.py` is used to make a calibration iteration reproducible.
Its `snapshot_sources` function preserves code, arguments and the exact reference
template before the sweep, so a later reader can repeat the comparison instead
of relying on a heatmap filename alone. Its support functions are documented under
[SVCW Section 3.4](SVCW_README.md#34-s4-shared-source-capture); a separate TF
archive format was not introduced.

## 4. Environment and Inputs

Reviewed runtime: Python 3.14.4 in `/opt/anaconda3/envs/algae`; installed metadata
QSDsan 1.6.0, EXPOsan 1.5.3, NumPy 2.4.6, SciPy 1.17.1, pandas 3.0.3,
Matplotlib 3.10.9, openpyxl 3.1.5 and Numba 0.65.1. Local source changes, not
version strings alone, define ESAP behavior. Archives record actual paths and
distribution versions. Arial is needed for exact figure typography.

Project profiles use 10 MGD, 20 degC and these concentrations in mg/L. The
explicit `default` profile is separate from low/mid/high.

| Strength | COD | NH4-N | PO4-P | K | Cl |
| --- | ---: | ---: | ---: | ---: | ---: |
| low | 339 | 14 | 1.6 | 11 | 39 |
| mid | 508 | 20 | 2.4 | 16 | 59 |
| high | 1016 | 41 | 4.7 | 32 | 118 |

COD fractions are SI=0.075, SF=0.24, SA=0.05, XI=0.30, with the remainder XS.
SF shares these feeds and the indices below. Startup reads
`../werf/data/initial_conditions.xlsx`, sheet E2P. The sweep additionally needs a
calculated GPS-X workbook with low/mid/high sheets and all 3/4/5-layer rows.
The default input path is the user's Google Drive TF_validation directory; use
an explicit path elsewhere.

The example supplies the 9464-m2/recycle-one comparison recipe. The ESAP factory
alone selects 200 m2, no recycle and H0=3e-4 m. Record overrides rather than
assuming validation, sweep and SVCW constructors select identical settings.
The validation factory's flag defaults are decay=True, kinetics=True and
AOB proxy=False; the examples therefore specify these options explicitly.

## 5. Key Assumptions and Conceptual Models

### Project Transport and Startup Changes

The section/layer layout follows the same conceptual approach as the
[GPS-X Technical Reference, Chapter 7, Trickling Filter Model](https://www.hydromantis.com/help/GPS-X/docs/8.5/Technical/7_CHAPTER7.htm).
This is a conceptual reference, not a claim of identical implementation or fit.
ESAP's changes are conservative transport/geometry, meaningful PAO seeds,
explicit ion transport and the detachment approximation below.

| Project Change | Setting and Reason |
| --- | --- |
| Startup | PAO=1% of X_H+X_PAO COD; PP/PAO=0.5*K_MAX; PHA/PAO=0.05, followed by dry-density normalization. Replaces near-zero arbitrary seeds that delayed guild establishment. |
| Ion transport | K/Mg/Na/Ca transport coefficient at 1e-5 cm2/s; Cl at 3.1e-5 cm2/s. Missing K/Mg transport coefficient previously allowed PP-associated consumption without replenishment. |
| Detachment | Capacity form, exponent 1, scale 0.065 kg dry solids/(m2 d). Selected for increasing loss near maximum thickness and the project's thickness/BOD fit; not a copied GPS-X rate equation. |
| Targeted accuracy | Concentration atol at most 1e-6 mg/L for selected TF nutrient/storage states and thickness atol at most 1e-9 m; tighter caller settings remain. Avoids large negative trace states without tightening all states equally. |
| Jacobian guard | Perturbation factors bounded at 1e-3 to avoid extreme inactive-state probes. This changes numerical evaluation, not reaction equations. |

Project bookkeeping uses `d(CV)/dt = V*dC/dt + C*dV/dt`; internal transfers
must cancel in section mass totals. Detachment scales with
`[x/max(1-x,0.001)]^p`, x=H/Hmax. Arbitrary hard resets of storage ratios were
not added because they would change conserved pools.
The meaningful seed changes the initial elemental inventory; it is not an
additional source term during integration.

### Project Component and Kinetic Overrides

| Components | i_N [g N/g COD] | i_P [g P/g COD] | Other Override |
| --- | ---: | ---: | --- |
| S_F, X_S | 0.023464 | 0.006149 | X_S i_mass=0.63 |
| S_I, X_I | 0.042021 | 0.007139 | X_I i_mass=0.63 |
| X_H, X_PAO, X_AUT | 0.060305 | 0.02369 | None |

These non-default indices were adopted for the project's Table-S1 comparison.
Keep the original table and supplied GPS-X screenshots with publication evidence;
they are not bundled source documents. The following table lists **only ESAP
overrides**, not the unchanged biological parameter set. Screenshot-derived
values are benchmark choices, not universal vendor defaults.

| ESAP Option | Parameters Changed and Purpose |
| --- | --- |
| `gpsx_kinetics` | mu_H=3.2/d; K_F=K_A_H=5 mg COD/L; eta_NO3=0.28; q_PHA=6 g COD/(g PAO-COD d); q_PP=2 g P/(g PAO-COD d); b_PP=0.05/d; eta_NO3_PAOl=0.9; b_AUT=0.17/d; eta_NO3_AUTl=0.5. Makes the directly mapped benchmark profile explicit. |
| `gpsx_decay` | b_H=0.62/d, eta_NO3_Hl=0.9, f_XI_H=0.08, f_SI=0; anaerobic reduction factor 0.6. Adds the project's smooth anaerobic **heterotrophic** decay contribution, not a PAO decay extension. |
| `gpsx_aob_proxy` | mu_AUT=0.9/d, K_NH4_AUT=0.7 mg N/L, K_O2_AUT=0.25 mg O2/L, b_AUT=0.17/d. Declares a lumped approximation rather than claiming separate nitrifier populations. |
| Calibrated four-parameter recipe | mu_H=3.3/d, q_PHA=6, q_PP=2, b_PP=0.05/d. Retains the current project fit independently of optional flag profiles. |

The added anaerobic term is proportional to
`b_H*X_H*eta_anaerobic*I_O2*I_NO3`, using smooth inhibition gates. Caller
`asm_kwargs` win over flags; an explicitly supplied TF model wins over the
constructed model. No complete MANTIS3 extension is claimed.

### Comparison and Energy Conventions

Heatmap comparison is `D=(QSDsan-GPSX)/max(abs(QSDsan),abs(GPSX))`, with zero
for zero/zero. NOx and electricity terms are summed first; the nitrite column is
zero because no separate nitrite result is produced by this project setting.

Project power accounting uses `P=rho*g*Q*head/eta` with conversion to kW.
Selected assumptions are 50+3.28 ft inlet head, bed depth+(5.5+3.28) ft recycle
head, 60% efficiency, and 0.003 kW/m3 media for the distributor. These are
accounting assumptions, not biological calibration. Daily electricity is P*24.
Optional stopping tests the full-window range, not endpoint equality; the sweep
turns stopping off for full-duration comparison.

## 6. Reproducibility, Outputs, and Limits

The sweep creates timestamped `.xlsx`, `_heatmap.tif` and `_sources.zip` files.
`run_settings` records effective arguments, archive/hash and actual end times.
References/formula caches remain in their original workbook format. Save tracked
histories and raw `sol.t/sol.y` separately when diagnosing dynamics; a heatmap
alone does not preserve a trajectory.

For a single case, call `snapshot_sources` before construction with explicit
JSON-compatible arguments. Shared instructions are in
[SVCW Section 6](SVCW_README.md#6-reproducibility-outputs-and-limits).
The [review record](REVIEW_20261015.md) documents 103 passing focused tests and
unchanged short cleanup trajectories; it is not a new 300-day validation.

```bash
MPLBACKEND=Agg python -m pytest -q tests/test_esap_trickling_filter.py
```

Project limits: short stationary windows can miss later transitions; inventory
stopping permits layer rearrangement. Kinetic proxies and
detachment fits do not establish complete GPS-X equivalence. Workbook tests use
mocked long case solves, not nine fresh 300-day calibration checks.
