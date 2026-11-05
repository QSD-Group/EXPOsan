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
