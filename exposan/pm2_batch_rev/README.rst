=====================================================================================
pm2_batch: Calibration of Phototrophic-Mixotrophic Process Model (PM2)
=====================================================================================

Summary
-------
This module is used to calibrate the Phototrophic-Mixotrophic Process Model (PM2; ``qsdsan.processes.PM2``) using batch experiment results. Note that the package `optuna <https://optuna.org/>`_ is needed to run the calibration script.

.. figure:: ./readme_figures/pm2_batch.svg

   *Layout of the system used for the calibration of PM2 implemented in QSDsan*

- ``system.py`` constructs the batch system above.
- ``model.py`` creates a model for uncertainty and sensitivity analysis.
- ``analyses.py`` runs uncertainty and sensitivity analysis and plots the results.
- ``calibration.py`` calibrates kinetic model parameters.
- ``_batch_may_unit_cali.py`` is used to check the calibration restuls.


Load the system
---------------
.. code-block:: python

   >>> from exposan import pm2_batch
   >>> pm2_batch.load()
   >>> sys = pm2_batch.sys
   >>> sys.simulate(t_span=(0,1), method='RK23')
   >>> PBR = pm2_batch.PBR
   >>> fig, axis = PBR.scope.plot_time_series(('S_P'))
   >>> fig

.. figure:: ./readme_figures/pm2_batch_phosphorus.png

    *Phosphorus removal in the batch reactor setting*

