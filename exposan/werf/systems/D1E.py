# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Zixuan Wang <wyatt4428@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

from .D1 import create_d1_system
from exposan.werf.utils._chp import create_chp_system

__all__ = ('create_d1e_system',)


def create_d1e_system(flowsheet=None, default_init_conds=True):
    """Create D1 with the E-suffix CHP accounting used by WERF models.

    Parameters
    ----------
    flowsheet : qsdsan.Flowsheet, optional
        Flowsheet to populate; a new D1E flowsheet is used by default.
    default_init_conds : bool, optional
        Apply D1's TF and digester initial concentrations. Defaults to True.

    Notes
    -----
    Reuses D1's treatment train, feed, components, kinetics and TF-specific
    solver safeguards unchanged. The D1E ID enables CHP energy recovery and
    operating-cost credits when werf.models.add_OPEX_metrics is called;
    no physical CHP unit or CHP capital cost is added by this factory.
    """
    return create_chp_system('D1E', create_d1_system, flowsheet, default_init_conds)
