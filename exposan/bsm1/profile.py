#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 20:19:58 2021

@author: Yalin Li
"""

import os, sys, io
os.environ['NUMPY_EXPERIMENTAL_ARRAY_FUNCTION'] = '0'

import numpy as np
import cProfile, pstats
from pstats import SortKey
from exposan.bsm1 import bsm1, bsm1_path


# =============================================================================
# Run cProfile for systems.py, print results in the console
# =============================================================================

t = 50
t_step = 0.5

# Looks like `tuna` isn't able to show all the details
# https://github.com/nschloe/tuna
# cProfile.run(f'bsm1.simulate(t_span = (0, {t}))',
#               f"systems_{''.join(str(t).split('.'))}.prof")

# Instead using `pstats`, for t = 0.1, the stats are stored in `systems_01`
cProfile.run(f'bsm1.simulate(t_span=(0, {t}), t_eval=np.arange(0, {t+t_step}, {t_step}))',
             f"systems_{''.join(str(t).split('.'))}.prof")

p = pstats.Stats(f"systems_{''.join(str(t).split('.'))}.prof")

# # If want to load previously saved stats
# p = pstats.Stats(f"systems_{''.join(str(t).split('.'))}.prof")


# # Print all stats
# # because of the many lines that will be printed, might consider using `%clear`
# # to clear up the console before printing
# p.strip_dirs().sort_stats(-1).print_stats()

# # Sort by the functions' names and print
# p.sort_stats(SortKey.NAME)
# p.print_stats()

# Sort by cumulative time, print the N most significant ones
p.sort_stats(SortKey.CUMULATIVE).print_stats(100)

# Sort by time spent within each function, print top N
p.sort_stats(SortKey.TIME).print_stats(20)

# Sort by total time spent on the function itself, excluding subcalls
p.sort_stats('tottime').print_stats(20)


# # Dump stats, looks like it doesn't matter what the extensions are
# # (or whether there's an extension),
# # `.txt` can be used, though won't be able to open it directly
# p.dump_stats(os.path.join(dir_path, f"systems_{''.join(str(t).split('.'))}.prof"))


# %%

# =============================================================================
# Save all results to a csv file
# =============================================================================

output = io.StringIO()
# Save the old stream so we can swap back
console_strm = p.stream
p.stream = output

# Pass the outputs to `text`
p.sort_stats('tottime').print_stats()

# Add deliminator
text = output.getvalue()
text = 'ncalls' + text.split('ncalls')[-1]
text = '\n'.join([','.join(line.rstrip().split(None,5)) for line in text.split('\n')])

# Save the results to csv, here
csv_path = os.path.join(bsm1_path, f"systems_{''.join(str(t).split('.'))}.csv")
with open(csv_path, 'w+') as f:
    f.write(text)
    f.close()

# Swap back the output stream
p.stream = console_strm