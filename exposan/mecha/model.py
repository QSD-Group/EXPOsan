# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ga-Yeong Kim <gayeong1225@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np, qsdsan as qs
from scipy.interpolate import interp1d
from chaospy import distributions as shape
from qsdsan.utils import DictAttrSetter, ospath, load_data
from datetime import datetime

from exposan.mecha import (
    create_system,
    data_path,
    )

__all__ = (
    'sensitive_params',
    'import_measured_data',
    'create_model',
    )

modified_asm1_kwargs = dict(
    Y_A=0.24, Y_H=0.67, f_P=0.08, i_XB=0.08, i_XP=0.06,
    mu_H=6.0, K_S=20.0, K_O_H=0.2, K_NO=0.5, b_H=0.62,
    eta_g=0.8, eta_h=0.4, k_h=3.0, K_X=0.03, mu_A=0.8,
    K_NH=1.0, b_A=0.1, K_O_A=0.4, k_a=0.08, fr_SS_COD=0.75,
    path=None,
    )

sensitive_params = {
    'Y_A': (0.24, 'g COD/g N', (0.1, 1)),
    'Y_H': (0.67, 'g COD/g COD', (0.1, 1)),
    # 'f_P': (0.08, '', (0.008, 0.8)),
    # 'i_XB': (0.08, 'g N/g COD', (0.008, 0.8)),
    # 'i_XP': (0.06, 'g N/g COD', (0.006,0.6)),
    'mu_H': (6.0, 'd^(-1)', (1, 10)),
    'K_S': (20.0, 'g COD/m^3', (10, 50)),
    'K_O_H': (0.2, 'g O2/m^3', (0.02, 2)),
    'K_NO': (0.5, 'g N/m^3', (0.05, 5)),
    'b_H': (0.62, 'd^(-1)', (0.062, 6.2)),
    'eta_g': (0.8, '', (0.1, 1)),
    'eta_h': (0.4, '', (0.1, 1)),
    'k_h': (3.0, 'd^(-1)', (0.3, 30)),
    'K_X': (0.03, 'g COD/g COD', (0.003, 0.3)),
    'mu_A': (0.8, 'd^(-1)', (0.08, 8)),
    'K_NH': (1.0, 'g N/m^3', (0.1, 10)),
    'b_A': (0.1, 'd^(-1)', (0.01, 1)),
    'K_O_A': (0.4, 'g O2/m^3', (0.04, 4)),
    'k_a': (0.08, 'm^3/g COD/d', (0.008, 0.8)),
    # 'fr_SS_COD': (0.75, 'g SS/g COD', (,)),
    }

# KLa_params = {
#     'FAC KLa': (30, '', (10, 50)),
#     'AER1 KLa': (60, '', (10, 100)),
#     'AER2 KLa': (60, '', (10, 100)),
#     'AER3 KLa': (60, '', (10, 100)),
#     }

#%%

def import_measured_data():

    file = ospath.join(data_path, 'train_val_test_online_dataset_calibration.xlsx')

    measured = load_data(file, sheet='train_val_test_online_dataset')
    # measured = load_data(file, sheet='train_val_test_online_dataset', index_col=None)

    day0 = datetime.strptime('5-5-2021 11:57 AM', '%m-%d-%Y %I:%M %p')
    calc_day = lambda t_stamp: (t_stamp-day0).total_seconds()/60/60/24

    t_stamp = measured.index
    t_measured = calc_day(t_stamp).to_numpy()

    nh_measured = measured.NH4.to_numpy()
    no_measured = measured.NO3.to_numpy()

    return t_measured, nh_measured, no_measured

#%%

def create_model(system=None, flowsheet=None):
    sys = create_system(asm_kwargs=modified_asm1_kwargs)

    stream = sys.flowsheet.stream
    DYINF, effluent, WAS, INT, RAS = stream.Dynamic_influent, stream.effluent, stream.WAS, stream.INT, stream.RAS

    unit = sys.flowsheet.unit
    WW, ANO, FAC, AER1, AER2, AER3, C1 = unit.Waste_Water, unit.ANO, unit.FAC, unit.AER1, unit.AER2, unit.AER3, unit.C1

    model = qs.Model(system=sys, exception_hook='warn')
    param = model.parameter
    metric = model.metric

    asm1 = ANO.suspended_growth_model
    cmps = ANO.components

    aerFac = FAC.aeration
    aer1 = AER1.aeration
    aer2 = AER2.aeration
    aer3 = AER3.aeration

    ##### General model with all uncertain variables #####

    # Add Uncertainty Parameters

    for k, v in sensitive_params.items():
        b, units, bounds = v
        D = shape.Uniform(lower=bounds[0], upper=bounds[1])
        param(setter=DictAttrSetter(asm1, '_parameters', keys=(k,)),
              name=k, element=ANO, kind='coupled', units=units, distribution=D, bounds=bounds)   # distribution=D (optional)
        # param(setter=DictAttrSetter(asm1.rate_function, '_params', k),
        #       name=k, element=RET, kind='coupled', units=units, distribution=D, bounds=bounds)   # distribution=D (optional)

    # for k, v in KLa_params.items():
    #     b, units, bounds = v
    #     D = shape.Uniform(lower=bounds[0], upper=bounds[1])

    #     if k == 'FAC KLa':
    #         @param(name=k, element=FAC, kind='coupled', units='',
    #                baseline=b, distribution=D)
    #         def set_FAC_KLa(i):
    #             aerFac.KLa = i

    #     elif k == 'AER1 KLa':
    #         @param(name=k, element=AER1, kind='coupled', units='',
    #                baseline=b, distribution=D)
    #         def set_AER1_KLa(i):
    #             aer1.KLa = i

    #     elif k == 'AER2 KLa':
    #         @param(name=k, element=AER2, kind='coupled', units='',
    #                baseline=b, distribution=D)
    #         def set_AER2_KLa(i):
    #             aer2.KLa = i

    #     else:
    #         @param(name=k, element=AER3, kind='coupled', units='',
    #                baseline=b, distribution=D)
    #         def set_AER3_KLa(i):
    #             aer3.KLa = i

    # # Aeration
    # b = aerFac.KLa
    # D = shape.Uniform(lower=10, upper=50)
    # @param(name='FAC KLa', element=FAC, kind='coupled', units='',
    #        baseline=b, distribution=D)
    # def set_Fac_KLa(i):
    #     aerFac.KLa = i

    # b = aer1.KLa
    # D = shape.Uniform(lower=10, upper=100)
    # @param(name='AER1 KLa', element=AER1, kind='coupled', units='',
    #        baseline=b, distribution=D)
    # def set_AER1_KLa(i):
    #     aer1.KLa = i

    # b = aer2.KLa
    # D = shape.Uniform(lower=10, upper=100)
    # @param(name='AER2 KLa', element=AER2, kind='coupled', units='',
    #        baseline=b, distribution=D)
    # def set_AER2_KLa(i):
    #     aer2.KLa = i

    # b = aer3.KLa
    # D = shape.Uniform(lower=10, upper=100)
    # @param(name='AER3 KLa', element=AER3, kind='coupled', units='',
    #        baseline=b, distribution=D)
    # def set_AER3_KLa(i):
    #     aer3.KLa = i

    ##### Add universal evaluation metrics #####

    idx_snh = cmps.indices(['S_NH'])
    idx_sno = cmps.indices(['S_NO'])

    measured_data = import_measured_data()
    t_measured, nh_measured, no_measured = measured_data

    @metric(element='Error')
    def RMSE():
        t_predicted = effluent.scope.time_series
        nh_predicted = effluent.scope.record[:,idx_snh].flatten()
        no_predicted = effluent.scope.record[:,idx_sno].flatten()

        f_nh = interp1d(t_predicted, nh_predicted)
        nh_predicted_intp = f_nh(t_measured)

        f_no = interp1d(t_predicted, no_predicted)
        no_predicted_intp = f_no(t_measured)

        eq_predicted = 30 * nh_predicted_intp + 10 * no_predicted_intp
        eq_measured = 30 * nh_measured + 10 * no_measured

        rmse = (sum((eq_predicted - eq_measured)**2)/len(t_measured))**0.5

        return rmse

    return model