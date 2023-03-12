# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
from qsdsan.utils import ospath, load_data
from exposan.metab import create_system, create_model, run_model, data_path
import qsdsan as qs

#%%
def run_discrete_DVs(samples_path):
    qs.PowerUtility.price = 0.0913
    dct = load_data(samples_path, sheet=None)
    for i in (
            # 'UASB', 
              # 'FB', 
              'PB',
             ):
        for n, j in (
                (1,'P'), 
                (2,'P'), 
                (2,'M'), 
                (2,'H'),
                ):
            sys = create_system(n_stages=n, reactor_type=i, gas_extraction=j)
            print(sys.ID)
            mdl = create_model(sys, kind='DV', exception_hook='warn')
            sample = dct[i+'_'+j].to_numpy()
            run_model(mdl, sample)

#%%
# def plot_clusters()


#%%
if __name__ == '__main__':
    path = ospath.join(data_path, 'analysis_framework.xlsx')
    run_discrete_DVs(path)
