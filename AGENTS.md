# Agent Context

EXPOsan is the applied systems repository. Keep case studies, scenarios, datasets, benchmark systems, and project-specific analyses here. Consume QSDsan through public APIs.

For paired QSDsan/EXPOsan architecture guidance, read the architecture skill for the agent runtime in use:

```text
../QSDsan/.codex/skills/qsdsan-exposan-architecture/SKILL.md
../QSDsan/.claude/skills/qsdsan-exposan-architecture/SKILL.md
```

If an EXPOsan unit or helper becomes broadly reusable outside one case study, promote it to QSDsan first, add QSDsan tests, then migrate EXPOsan to the public API.

## Module Template

Each EXPOsan system or system family should live in one top-level module. Use `system.py` when the module exposes one system, and `systems.py` when it exposes multiple systems or configurations. Keep existing modules stable, but use this shape for new work:

```text
exposan/<module>/
  __init__.py
  system.py or systems.py
  model.py or models.py
  _components.py           optional
  _units.py                optional, for project-specific units
  _process_settings.py     optional
  _tea.py or _lca.py       optional
  data/                    input data only
  results/                 created only when writing results
  figures/                 created only when writing figures
```

`__init__.py` should stay thin: paths, light public imports, `load()`, and lazy access errors. Build systems in `system.py` or `systems.py`; build uncertainty/sensitivity models in `model.py` or `models.py`; keep large constants in data files or focused settings modules.

Use `None` rather than mutable default dictionaries in public functions. For multi-system dispatch, prefer an explicit registry such as `SYSTEM_CREATORS = {'A': create_systemA}` over `globals()` lookup.

## Dependencies

EXPOsan's `pyproject.toml` declares only `qsdsan` as a required dependency; `scikit-learn`, `tsp_solver`, `ortools` live behind an `exposan[complete]` extra (the `ci`/`dev` extras self-reference `complete`, so CI still installs everything). Use `exposan.utils.require_package(pip_name, import_name=None, feature=None)` to lazily import an optional dependency and raise a friendly `ModuleNotFoundError` naming the exact pip package — don't invent a second helper. `hap` (routing depends on ortools/tsp_solver) requires `complete` as a whole module; its routing submodule is not independently lazy. QSDsan's own dependencies (`sympy` in particular — required unconditionally by `Process`/`CompiledProcesses`) are intentionally left untouched.

## Testing: division of labor with QSDsan

`QSDsan/tests/test_exposan.py` is a fast **construction-only smoke canary** ("does EXPOsan still wire up against this QSDsan?") — dynamic systems (`adm`, `asm`, `bsm1`, `bsm2`) call `load(simulate=False)` there and never run their ODE solve. **Full dynamic simulation coverage is EXPOsan's own job** (this repo's test suite + the weekly backstop), not QSDsan's. `cas` is steady-state (not dynamic) and always simulates on `load()` — don't add `simulate=False` to it. `metab`/`werf`/`hap` expose no `load()` at all.

Heavy stiff-ODE dynamic tests (`adm`/`asm`/`bsm1`/`bsm2`/`metab`/`pm2`/`werf`) run in `tests/conftest.py` via subprocess isolation (a `pytest_pyfunc_call` hook re-invokes `pytest <nodeid>` in a fresh process) to avoid cross-test numerical contamination discovered during BioSTEAM-master compat testing — do not remove this without re-verifying that contamination is actually gone.

**Flaky-looking dynamic-sim failures:** don't diagnose a stochastic test from a single run. Loop 6+ times and sweep `PYTHONHASHSEED=0..7` before attributing a cause — the werf I3 sim looked like a clean upstream BioSTEAM regression in a one-shot A/B but was actually seed-dependent on both versions (root cause: `flexsolve` sets `np.seterr(divide='raise', invalid='raise')` globally at import, turning harmless transient invalid/divide ops in a stiff BDF solve into hash-sensitive `FloatingPointError`s). Fixed upstream in `biosteam/_system.py::dynamic_run` via an `np.errstate` guard.

## Module-specific status (verify current state before relying on these)

- **metab**: two TEA/LCA "expected" baselines exist and are NOT interchangeable — `results/discrete_DVs/*.xlsx` (reproducible to ~0.5% by the current stack, a legitimate regression target) vs. the commented-out asserts in `tests/test_metab.py` (off by 4-9× even on the exact Feb-2024 stack that produced them; not a usable target). When debugging metab drift, confirm which baseline is being asked about first.
- **metab UASB+M**: uses ~400× sidestream recirculation for membrane degassing — a peculiar METAB-only operating point, not typical UASB. A resulting "pancake" reactor-sizing result is a legitimate study finding, not a bug (see the sizing-constraints-are-signals gotcha in `QSDsan/AGENTS.md`).
- **saf**: crude distillation made CI-robust via a deterministic Lr-ladder spec (`screen_results` replacement in `saf/systems.py`); `tests/test_saf.py` re-enabled.
- **biobinder**: crude distillation kept `simulate=False` / disabled — the ~53/47 target cut sits on a singular VLE region (heavy pseudo-components) that no column spec can reliably converge; needs a modeling change (pin the mass split, retain real reboiler/condenser cost), not a different spec. Don't re-attempt a spec-only fix without new information.
- **htl**: has (or had) a real upstream ordering bug in `systems.py` — two `ReversedSplitter` units executed before the units whose demand they reflect, one simulate-pass behind. Fixed by declaring their outs as `sys.recycle` plus a tight `sys.set_tolerance`. Verify this landed on `main` before assuming it's still open (branch was `fix-htl-recycle`, PR against QSD-Group/EXPOsan).
