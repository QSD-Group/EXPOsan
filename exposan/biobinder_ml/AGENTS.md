# Module conventions

For future Python work in this module, use the following header for new files
and add it when modifying a file without a project header. Preserve existing
contributor attribution in files that already have a header.

```python
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
```

Keep incentive cash-flow treatments distinct: loans and guarantees affect
financing; grants affect capital funding and potentially tax basis; tax credits
affect tax liability or realized credit proceeds; tradable credits are market
revenue; procurement preferences affect contracted sales; avoided-disposal
contracts are service income. Do not double count tipping fees already included
in negative feedstock prices.

Test files in this module are only for system-level files and their results.
Do not create dedicated test files for individual Python scripts or helpers.

In scripts, do not leave comments or docstrings describing material as supplied
by the user, recounting conversation history, or linking to local reference
files supplied during a task. For source attribution and modeling assumptions,
include only active source URLs and relevant reference values with units as
comments beside the applicable calculation. Keep functional API documentation
and the required project/license header.

Use `exposan.biobinder_ml._distill_utils_flex` as the single shared adaptive
distillation helper. Import it by its full package path so systems and reporting
share the same trial registry. Do not add helper copies in deploy or systems.
