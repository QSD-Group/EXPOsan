===========================================================================================
pm2_ecorecover: Calibration & Validation of Phototrophic-Mixotrophic Process Model (PM2)
===========================================================================================

Summary
-------
This module is used to calibrate and validate the Phototrophic-Mixotrophic Process Model (PM2; ``qsdsan.processes.PM2``) using the continuous online monitoring data from the EcoRecover plant. Note that the package `optuna <https://optuna.org/>`_ is needed to run the calibration script and the package `scikit-fda <https://pypi.org/project/scikit-fda/>`_ is needed to run the data cleaning script.

.. figure:: ./readme_figures/pm2_ecorecover.svg

   *Layout of the system used for the calibration and validation of PM2 implemented in QSDsan*

- ``system.py`` constructs EcoRecover system above.
- ``model.py`` creates a model for uncertainty and sensitivity analysis.
- ``calibration.py`` calibrates model parameters.
- ``data_cleaning.py`` cleans raw online sensor data.


Load the system
---------------
.. code-block:: python

   >>> from exposan import pm2_ecorecover
   >>> pm2_ecorecover.load()
   >>> sys = pm2_ecorecover.sys
   >>> sys.simulate(t_span=(0,3), method='RK23')
   >>> PBR = pm2_ecorecover.PBR20
   >>> fig, axis = PBR.scope.plot_time_series(('S_P'))
   >>> fig

.. figure:: ./readme_figures/pm2_ecorecover_phosphorus.png

    *Phosphorus removal in the EcoRecover system*

