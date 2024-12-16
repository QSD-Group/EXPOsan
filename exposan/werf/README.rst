=================================================================================
werf: 18 Benchmark Models for Water Resource Recovery Facilities (WRRF) in the US
=================================================================================

Summary
-------
This modules contains 18 benchmark models of typical WRRF configurations in the US [1]_, [2]_. Plant-wide dynamic simulations are enabled with all configurations.

** To be completed... **


Load the system
---------------
.. code-block:: python

	>>> from exposan.werf import create_g1_system
	>>> sys = create_g1_system()
	>>> # sys.simulate(method='BDF', t_span=(0, 300))
	>>> # Check the states of any one or more components
	>>> # fig, axis = bsm2.AD1.scope.plot_time_series(u.AD1.components.IDs[:3])
	>>> # fig


.. figure:: ./figures/G1.png

    *Current module diagram.*


References
----------
.. [1] Tarallo, S.; Shaw, A.; Kohl, P.; Eschborn, R. A Guide to Net-Zero Energy Solutions for Water Resource Recovery Facilities; ENER1C12; The Water Environment Research Foundation, 2015. `<https://iwaponline.com/ebooks/book/293/>`_.
.. [2] El Abbadi, S. H.; Feng, J.; Hodson, A. R.; Amouamouha, M.; Busse, M. M.; Polcuch, C.; Zhou, P.; Macknick, J.; Guest, J. S.; Stokes-Draut, J. R.; Dunn, J. B. Benchmarking Greenhouse Gas Emissions from U.S. Wastewater Treatment for Targeted Reduction. 2024. `<https://doi.org/10.31223/X5VQ59>`_
