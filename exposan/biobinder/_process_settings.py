#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Yalin Li <mailto.yalin.li@gmail.com>
    
    Ali Ahmad <aa3056@scarletmail.rutgers.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from exposan.saf import _process_settings
from exposan.saf._process_settings import *

__all__ = [i for i in _process_settings.__all__ if i is not 'dry_flowrate']
__all__.extend(['central_dry_flowrate', 'pilot_dry_flowrate'])

central_dry_flowrate = dry_flowrate # 110 tpd converted to kg/hr
pilot_dry_flowrate = 11.46 # kg/hr

# Salad dressing waste
HTL_yields = {
    'gas': 0.1756,
    'aqueous': 0.2925,
    'biocrude': 0.5219,
    'char': 1-0.1756-0.2925-0.5219,
    }

# https://idot.illinois.gov/doing-business/procurements/construction-services/transportation-bulletin/price-indices.html
# bitumnous, IL
price_dct['biobinder'] = 0.67