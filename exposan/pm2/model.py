# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from scipy.interpolate import interp1d
# from warnings import warn
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter, ospath, load_data

# from qsdsan.utils import DictAttrSetter, AttrGetter, FuncGetter, \
#     get_SRT as srt, ospath, time_printer

from exposan.pm2 import (
    create_system,
    data_path,
    )

__all__ = (
    'sensitive_params',
    'import_scada_data',
    'create_model',
    )

sensitive_params = {
    'arr_e': (6842, 'K', (1000, 10000)),
    'K_P': (1.0, 'g P/m^3', (0.0001, 10)),
    'K_STO': (1.566, 'g COD/g COD', (0.0001, 20)),
    'exponent': (4, '', (1, 10)),
    'm_ATP': (15.835, 'g ATP/g COD/d', (1, 50)),
    'q_CH': (0.594, 'g COD/g COD/d', (0.001, 10)),
    'V_NH': (0.254, 'g N/g COD/d', (0.0001, 5)),
    'V_P': (0.016, 'g P/g COD/d', (0.0001, 5))
    }

#%%
def import_scada_data():

    file = ospath.join(data_path, 'condti_scada_result.xlsx')
    result_scada = load_data(file, sheet=0, index_col=None)

    t_scada = result_scada.t_stamp.to_numpy()

    vss_scada = result_scada.VSS.to_numpy()
    snh_scada = result_scada.S_NH.to_numpy()
    sp_scada = result_scada.S_P.to_numpy()
    # sp_scada = np.maximum(result_scada.S_P.to_numpy(), 0.015)

    return t_scada, vss_scada, snh_scada, sp_scada

#%%

def create_model(system=None, flowsheet=None):
    # sys = create_system(flowsheet)
    sys = system or create_system(flowsheet)

    stream = sys.flowsheet.stream
    DYINF, PHO, ME, INT, TE, RETEN, RE, CE, CEN, ALG, RAA = stream.DYINF, stream.PHO, stream.ME, stream.INT, stream.TE, stream.RETEN, stream.RE, stream.CE, stream.CEN, stream.ALG, stream.RAA

    unit = sys.flowsheet.unit
    SE, MIX, PBR1, PBR2, PBR3, PBR4, PBR5, PBR6, PBR7, PBR8, PBR9, PBR10, PBR11, PBR12, PBR13, PBR14, PBR15, PBR16, PBR17, PBR18, PBR19, PBR20, MEM, MEV, POST_MEM, CENT, RET = unit.SE, unit.MIX, unit.PBR1, unit.PBR2, unit.PBR3, unit.PBR4, unit.PBR5, unit.PBR6, unit.PBR7, unit.PBR8, unit.PBR9, unit.PBR10, unit.PBR11, unit.PBR12, unit.PBR13, unit.PBR14, unit.PBR15, unit.PBR16, unit.PBR17, unit.PBR18, unit.PBR19, unit.PBR20, unit.MEM, unit.MEV, unit.POST_MEM, unit.CENT, unit.RET

    model = qs.Model(system=sys, exception_hook='warn')
    param = model.parameter
    metric = model.metric

    pm2 = MIX.suspended_growth_model

    cmps = MIX.components

    ##### General model with all uncertain variables #####

    # Add Uncertainty Parameters

    for k, v in sensitive_params.items():
        b, units, bounds = v
        D = shape.Uniform(lower=bounds[0], upper=bounds[1])
        param(setter=DictAttrSetter(pm2.rate_function, '_params', k),
              name=k, element=RET, kind='coupled', units=units, distribution=D, bounds=bounds)   # distribution=D (optional)

    ##### Add universal evaluation metrics #####

    idx_vss = cmps.indices(['X_ALG', 'X_CH', 'X_LI', 'X_N_ALG', 'X_P_ALG'])
    idx_snh = cmps.indices(['S_NH'])
    idx_sp = cmps.indices(['S_P'])

    imass = cmps.i_mass[idx_vss]

    # import_scada_data
    scada_data = import_scada_data()

    t_scada, vss_scada, snh_scada, sp_scada = scada_data

    vss_range = max(vss_scada) - min(vss_scada)
    snh_range = max(snh_scada) - min(snh_scada)
    sp_range = max(sp_scada) - min(sp_scada)

    @metric(element='Error')
    def NRMSE_VSS():

        t_simul = RET.scope.time_series
        result_simul_vss = RET.scope.record[:,idx_vss]

        vss = np.sum(result_simul_vss * imass, axis = 1)
        # vss = result_simul[:,0] * cmps.X_ALG.i_mass + result_simul[0:,1] * cmps.X_CH.i_mass + result_simul[0:,2] * cmps.X_LI.i_mass + result_simul[0:,3] + result_simul[0:,4]

        f = interp1d(t_simul, vss)
        vss_simul = f(t_scada)

        rmse_vss = (sum((vss_scada - vss_simul)**2)/len(t_scada))**0.5
        nrmse_vss = rmse_vss/vss_range
        return nrmse_vss

    @metric(element='Error')
    def NRMSE_SNH():
        t_simul = TE.scope.time_series
        result_simul_snh = TE.scope.record[:,idx_snh].flatten()

        f = interp1d(t_simul, result_simul_snh)
        snh_simul = f(t_scada)

        rmse_snh = (sum((snh_scada - snh_simul)**2)/len(t_scada))**0.5
        nrmse_snh = rmse_snh/snh_range
        return nrmse_snh

    @metric(element='Error')
    def NRMSE_SP():
        t_simul = TE.scope.time_series
        result_simul_sp = TE.scope.record[:,idx_sp].flatten()

        f = interp1d(t_simul, result_simul_sp)
        sp_simul = f(t_scada)

        rmse_sp = (sum((sp_scada - sp_simul)**2)/len(t_scada))**0.5
        nrmse_sp = rmse_sp/sp_range
        return nrmse_sp

    return model