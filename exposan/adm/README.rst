==================================================
adm: Validation of Anaerobic Digestion Model No. 1
==================================================

Summary
-------
This module is used to validate the implementation of Anaerobic Digestion Model No. 1 (ADM1) [1]_ in QSDsan using a continuous stirred-tank reactor (CSTR). The implemented algorithms (i.e., ``qsdsan.processes.ADM1`` and  ``qsdsan.sanunits.AnaerobicCSTR`` classes) has been verified against the MATLAB/Simulink implementation with algebraic solutions for H2 and pH developed by Flores-Alsina et al. [2]_ , [3]_.

 .. figure:: ./figures/ADM1_cstr.png

    *Layout of the system used for the validation of ADM1 implemented in QSDsan.*

- ``system.py`` constructs the system above.
- ``model.py`` creates a model with identical parameter settings and influent composition as in Rosen and Jeppsson [3]_ and performs simulations of the system with 100 randomly sampled initial conditions until the steady state is reached.
- ``analyses.py`` processes and plots the simulated time-series data of the system state variables with steady-state values reported in Rosen and Jeppsson [3]_.

References
----------
.. [1] IWA Task Group for Mathematical Modelling of Anaerobic Digestion Processes. Anaerobic Digestion Model No.1 (ADM1); IWA Publishing, 2005. `<https://doi.org/10.2166/9781780403052>`_
.. [2] wwtmodels. Anaerobic-Digestion-Models; Github (Accessed July 20, 2022) `<https://github.com/wwtmodels/Anaerobic-Digestion-Models>`_
.. [3] Rosen and Jeppsson. Aspects on ADM1 Implementation within the BSM2 Framework; Lund, 2006. `<https://www.iea.lth.se/publications/reports/lth-iea-7224.pdf>`_
