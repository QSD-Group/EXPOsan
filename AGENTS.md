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
