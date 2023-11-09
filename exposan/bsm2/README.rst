======================================
bsm2: Benchmark Simulation Model No. 2
======================================



BELOW IS OUTDATED


==============================================================
interface: Demo System for Dynamic Simulation with ASM and ADM
==============================================================

Summary
-------
This modules contains a simple system to demonstrate dynamic simulation with more than one process models by including interface models (``ADMtoASM`` and ``ASMtoADM`` for anaerobic digestion model no.1 and activated sludge model no.1 [1]_) .

The system contains the entire Benchmark Simulation Model No. 1 (BSM1) [2]_ configuration as modeled in the `bsm1 <https://github.com/QSD-Group/EXPOsan/tree/main/exposan/bsm1>`_ module. But the waste activated sludge from the clarifier is sent to an anaerobic digester (AD), and the effluent from the AD is recycled into the first activated sludge bioreactor.

.. figure:: ./readme_figures/interface.svg

    *Interface system layout: BSM1 with an anaerobic digester*


Load the system
---------------
.. code-block:: python

	>>> from exposan import interface as inter
	>>> inter.load()
	>>> sys = inter.sys
	>>> sys.simulate(method='BDF', t_span=(0, 15)) # the default 'RK45' method can't solve it
	>>> u = sys.flowsheet.unit # unit registry
	>>> # Check the states of any one or more components
	>>> fig, axis = u.AD1.scope.plot_time_series(u.AD1.components.IDs[:3])
	>>> fig


.. figure:: ./readme_figures/AD_select_states.png

    *Select component states from the anaerobic digester in the interface system*


References
----------
.. [1] Nopens et al., An ASM/ADM Model Interface for Dynamic Plant-Wide Simulation. Water Res. 2009, 43, 1913–1923. `<https://doi.org/10.1016/j.watres.2009.01.012>`_
.. [2] Alex et al., Benchmark simulation model no. 1 (BSM1). Report by the IWA Taskgroup on benchmarking of control strategies for WWTPs (2008): 19-20. `<http://iwa-mia.org/benchmarking/#BSM1>`_