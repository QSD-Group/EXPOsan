# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
import numpy as np
from qsdsan.utils import auom
from math import pi

__all__ = ('pipe_design', 'pipe_friction_head', 'hdpe_price')

#%%
_m2in = auom('m').conversion_factor('inch')
_ft2m = auom('ft').conversion_factor('m')
_lbft2kgm = auom('lbs/ft').conversion_factor('kg/m')

_Hazen_Williams_coefficients = {
    'Stainless steel': 110,
    'HDPE': 140
    }

# https://www.piping-designer.com/index.php/datasheets/piping-datasheets/1651-pipe-hdpe-ansi-dr-11-0-ips-in
hdpe_ID_inch = np.array([
    1.340, 1.533, 1.917, 2.826, 3.633, 5.349, 6.963, 8.679,
    10.293, 11.301, 12.915, 14.532, 16.146, 17.760, 19.374, 
    20.988, 22.605, 24.219, 25.833, 27.447, 29.061, 33.906
    ])
hdpe_lb_per_ft = np.array([
    0.31, 0.41, 0.64, 1.39, 2.31, 5.00, 8.47, 13.16,
    18.51, 22.32, 29.15, 36.89, 45.54, 55.1, 65.58, 
    76.96, 89.26, 102.47, 116.58, 131.61, 147.55, 200.84
    ])

# Schedule 10 pipes @ https://amerpipe.com/wp-content/uploads/2015/10/APP-chart-v7-web.pdf
ssteel_ID_inch = np.array([
    0.307, 0.410, 0.545, 0.674, 0.884, 1.097, 1.442, 1.682, 
    2.157, 2.635, 3.260, 3.760, 4.260, 5.295, 6.357, 8.329, 
    10.42, 12.39, 13.50, 15.50, 17.50, 19.50, 21.50, 23.50
    ])
ssteel_lb_per_ft = np.array([
    0.049, 0.065, 0.065, 0.083, 0.083, 0.109, 0.109, 0.109, 
    0.109, 0.120, 0.120, 0.120, 0.120, 0.134, 0.134, 0.148, 
    0.165, 0.180, 0.250, 0.250, 0.250, 0.250, 0.250, 0.250
    ])


def pipe_design(F_vol, vmin, stainless_steel=False):
    '''Pipe sizing based on flow rate.'''
    ID = (F_vol/vmin/pi)**(1/2) * 2 * _m2in
    if stainless_steel: 
        IDs = ssteel_ID_inch
        weights = ssteel_lb_per_ft
    else:
        IDs = hdpe_ID_inch
        weights = hdpe_lb_per_ft
    ids = IDs[IDs <= ID]
    if ids.size == 0: 
        ID = IDs[0] # inch
        kg_per_m = weights[0] * _lbft2kgm
    else: 
        ID = ids[-1]
        kg_per_m = weights[IDs == ID][0] * _lbft2kgm
    return ID, kg_per_m

def pipe_friction_head(q, L, ID, stainless_steel=False):
    '''Hazen-Williams equation. 
    https://www.engineeringtoolbox.com/hazen-williams-water-d_797.html'''
    if stainless_steel: c = _Hazen_Williams_coefficients['Stainless steel']
    else: c = _Hazen_Williams_coefficients['HDPE']
    return 2.083e-3 * (100*q / c)**1.852 / ID**4.8655 * L

def hdpe_price(ID):
    '''Price in [USD/kg] as a function of inner diameter [inch],
    projection based on prices in https://hdpesupply.com/hdpe-straight-length-pipe/'''
    return 9.625*ID**(-0.368)

