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

from qsdsan.utils import ospath, time_printer
from exposan.pm2 import (
    results_path,
    create_model,
    sensitive_params,
    )

import numpy as np, pandas as pd
from scipy.optimize import minimize, basinhopping, shgo

__all__ = ('cali_setup', 'optimizer', 'objective_function')

#%%

mdl = create_model()

def cali_setup():
    params = []
    baseline = np.array([])
    boundary = []

    for k, v in sensitive_params.items():
        b, units, bounds = v
        params.append(k)
        baseline = np.append(baseline, b)
        boundary.append(bounds)
    boundary = tuple(boundary)

    return params, baseline, boundary

sense_params, opt_params, bnds = cali_setup()         # start from init_guess

'''
sense_params: list of sensitive parameters
opt_params: initial guess of sensitive parameters
bnds: min & max of sensitive parameters
'''

#%%
def optimizer():

    opt = shgo(objective_function, bounds=bnds, iters=2, minimizer_kwargs={'method':'SLSQP', 'ftol':1e-3})   #2_with gygutfeeling_init->error

    # opt = shgo(objective_function, bounds=bnds, iters=1, minimizer_kwargs={'method':'SLSQP'})   #1_with gygutfeeling_init


    # opt = shgo(objective_function, bounds=bnds, iters=1, minimizer_kwargs={'method':'SLSQP', 'ftol':1e-3}, sampling_method='simplicial')
    # opt = basinhopping(objective_function, opt_params, niter=3, minimizer_kwargs={'method':'SLSQP'}, niter_success=5)

    # opt = minimize(objective_function, opt_params, method='SLSQP', bounds=bnds, tol=1e-3)

    # opt = minimize(objective_function, opt_params, method='SLSQP', bounds=bnds, tol=1e-4)

    # tol=1e-6 warning, fail
    # tol=1e-5 warning, fail
    # tol=1e-4 warning, success

    opt_as_series = pd.Series(opt)
    opt_as_series.to_excel(excel_writer=(ospath.join(results_path, 'calibration_result_newbase_conti_iter2.xlsx')))

    # opt_as_series.to_excel(excel_writer=(ospath.join(results_path, 'calibration_result_exclude.xlsx')))

    # scipy.optimize.minimize(fun, x0, args=(), method=None, jac=None, hess=None, hessp=None, bounds=None,\
    #                         constraints=(), tol=None, callback=None, options=None)

    # scipy.optimize.shgo(func, bounds, args=(), constraints=None, n=None,
    #                     iters=1, callback=None, minimizer_kwargs=None, options=None, sampling_method='simplicial')

    # scipy.optimize.basinhopping(func, x0, niter=100, T=1.0, stepsize=0.5, minimizer_kwargs=None, take_step=None,
    #                             accept_test=None, callback=None, interval=50, disp=False, niter_success=None, seed=None, *, target_accept_rate=0.5, stepwise_factor=0.9)

#%%
'''
<scipy.optimize.minimize>
fun: The objective function to be minimized.
x0: Initial guess. Array of real elements of size (n,), where n is the number of independent variables.
args: Extra arguments passed to the objective function and its derivatives (fun, jac and hess functions).
method: Type of solver. Should be one of ‘Nelder-Mead’,‘Powell’,‘CG’,‘BFGS’,‘Newton-CG’,‘L-BFGS-B’,‘TNC’,‘COBYLA’,‘SLSQP',‘trust-constr’,‘dogleg’,‘trust-ncg',‘trust-exact’,‘trust-krylov’
jac: Method for computing the gradient vector. Only for CG, BFGS, Newton-CG, L-BFGS-B, TNC, SLSQP, dogleg, trust-ncg, trust-krylov, trust-exact and trust-constr.
hessp: Hessian of objective function times an arbitrary vector p. Only for Newton-CG, trust-ncg, trust-krylov, trust-constr.
bounds: Bounds on variables for Nelder-Mead, L-BFGS-B, TNC, SLSQP, Powell, and trust-constr methods.
constraints: Constraints definition. Only for COBYLA, SLSQP and trust-constr. Constraints for ‘trust-constr’ are defined as a single object or a list of objects specifying constraints to the optimization problem.
tol: Tolerance for termination.
callback: Called after each iteration.
res: The optimization result represented as a OptimizeResult object.

<scipy.optimize.OptimizeResult>
x: The solution of the optimization
success: Whether or not the optimizer exited successfully.
status: Termination status of the optimizer. Its value depends on the underlying solver. Refer to message for details.
message: Description of the cause of the termination.
fun, jac, hess, hess_inv: Values of objective function, Jacobian, Hessian or its inverse (if available). The Hessians may be approximations, see the documentation of the function in question.
nfev, njev, nhev: Number of evaluations of the objective functions and of its Jacobian and Hessian.
nit: Number of iterations performed by the optimizer.
'''
#%%
@time_printer
def objective_function(opt_params, *args):
    
    # for p,v in zip(mdl.parameters, opt_params):
    #     p.setter(v)
    
    # sys = mdl.system
    # sys.simulate(t_span=(0, 50), t_eval = np.arange(0, 51, 1), method='RK23', state_reset_hook='reset_cache', print_t=True)

    # sys.simulate(t_span=(0, 50), t_eval = np.arange(0, 51, 1), method='RK23', state_reset_hook='reset_cache', print_t=True)

    try:
        mdl._update_state(opt_params, t_span=(0, 50), t_eval = np.arange(0, 51, 1), method='RK23', state_reset_hook='reset_cache', print_t=True)
        # mdl._update_state(opt_params, t_span=(0, 7), t_eval = np.arange(0, 7.01, 0.01), method='BDF', state_reset_hook='reset_cache')

    except:
        return 0.5

    out = [metric() for metric in mdl.metrics]
    obj = np.average(out)
    
    print('succeed')

    return obj

optimizer()