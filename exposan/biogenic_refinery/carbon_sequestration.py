#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd, qsdsan as qs
import models as m
from chaospy import distributions as shape
from qsdsan import Model, Metric
from exposan import biogenic_refinery as br
from exposan.biogenic_refinery import create_system, results_path, run_uncertainty

# Filter out warnings related to uptime ratio
import warnings
warnings.filterwarnings('ignore', message='uptime_ratio')

# input data for changing ash content and temperature
input_dct = {
    'LitData': {'f_ash_content': ('triangle', 17, 38.45, 58.7),
                'pyrolysis_temp': ('uniform', 400, 600, 800)},
    'SysA': {'f_ash_content': ('uniform', 46.08, 51.2, 56.32),
             'pyrolysis_temp': ('uniform', 400, 600, 800)},
    'SysB': {'f_ash_content': ('uniform', 34.02, 37.8, 41.58),
             'pyrolysis_temp': ('uniform', 400, 600, 800)}
    }

# =============================================================================
# Functions to create models (only working with sysA and sysB here)
# =============================================================================

def create_modelA(country_specific=False, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysA = br.create_system('A', flowsheet=flowsheet)
    modelA = Model(sysA, **model_kwargs)
    return modelA

def create_modelB(country_specific=False, **model_kwargs):
    flowsheet = model_kwargs.pop('flowsheet', None)
    sysB = br.create_system('B', flowsheet=flowsheet)
    modelB = Model(sysB, **model_kwargs)
    return modelB

# =============================================================================
# Function to run models for biochar/carbon sequestration metrics
# =============================================================================

def run_CS(dct):
    results = []    
    for sys_name, sys_dct in dct.items():
        
        # create model according to system
        if sys_name == 'SysB':
            model = create_modelB(country_specific=False)
            sys = model.system
            u = sys.flowsheet.unit
            pyrolysis = u.B11
        else:
            model = create_modelA(country_specific=False)
            sys = model.system
            u = sys.flowsheet.unit            
            pyrolysis = u.A8
                    
        # create metrics related to biochar production
        CS_metrics = [
            pyrolysis.yield_db,
            pyrolysis.outs[0].F_mass * 24,
            pyrolysis.CS,
            pyrolysis.outs[0].imass['C'] * 24
            ]
        model.metrics = [
            Metric('Biochar Yield', CS_metrics[0], '% db'),
            Metric('Mass Biochar Produced', CS_metrics[1], 'kg/day'),
            Metric('Carbon Sequestration Potential', CS_metrics[2], '%'),
            Metric('Mass Carbon Sequestered', CS_metrics[3], 'kg/day')
            ]
            
        # set parameters
        param = model.parameter
        
        # feedstock ash content
        kind, low_val, peak_val, max_val = sys_dct['f_ash_content']
        b = peak_val
        if kind == 'triangle':
            D = shape.Triangle(lower=low_val, midpoint=peak_val, upper=max_val)
        else:
            D = shape.Uniform(lower=low_val,upper=max_val)
        @param(name='Feedstock ash content', element=pyrolysis, kind='coupled', units='%',
               baseline=b, distribution=D)
        def set_ash_content(i):
            pyrolysis.f_ash_content = i
            
        # pyrolysis temperature
        kind, low_val, peak_val, max_val = sys_dct['pyrolysis_temp']
        b = peak_val
        if kind == 'triangle':
            D = shape.Triangle(lower=low_val, midpoint=peak_val, upper=max_val)
        else:
            D = shape.Uniform(lower=low_val,upper=max_val)
        @param(name='Pyrolysis temperature', element=pyrolysis, kind='coupled', units='deg C',
               baseline=b, distribution=D)
        def set_pyrolysis_temp(i):
            pyrolysis.pyrolysis_temp = i

        results[sys_name] = m.run_uncertainty(model)
        del sys, model
        
    return results

results = run_CS(dct=input_dct)
