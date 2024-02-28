========================================================
pou_disinfection: Point-Of-Use Disinfection Technologies
========================================================

Summary
-------
Systems included in this module are developed for the manuscript Elijah et al. [1]_ for four point-of-use (POU) disinfection technologies.

Four system configurations are included in this module, each corresponds to one technology.

.. figure:: ./readme_figures/sysA.svg

    *Technology A (sysA): POU chlorination*

.. figure:: ./readme_figures/sysB.svg

    *Technology B (sysB): Silver nanoparticle-enabled ceramic water filter (AgNP CWF)*

.. figure:: ./readme_figures/sysC.svg

    *Technology C (sysC): Low-pressure UV mercury lamp system*

.. figure:: ./readme_figures/sysD.svg

    *Technology D (sysD): UV LED disinfection lamp system*


Loading systems
---------------
.. code-block:: python

    >>> # Import and load the system
    >>> from exposan import pou_disinfection as pou
    >>> pou.load()
    >>> # Quick look at the systems
    >>> pou.sysA.show() # doctest: +ELLIPSIS
	System: sysA
	ins...
	[0] A_naclo  
	    phase: 'l', T: 298.15 K, P: 101325 Pa
	    flow (kmol/hr): NaClO  1.55e-08
	[1] A_cl_bottle  
	    phase: 's', T: 298.15 K, P: 101325 Pa
	    flow (kmol/hr): Polyethylene  3.95e-09
	outs...
	[0] treated_water  
	    phase: 'l', T: 298.15 K, P: 101325 Pa
	    flow (kmol/hr): Mg     1.9e-06
	                    Ca     1.15e-06
	                    Ecoli  0.308


Techno-economic analysis and life cycle assessment (TEA/LCA)
------------------------------------------------------------
TEA and LCA are performed through :class:`TEA` and :class:`LCA` objects and named as "teaA", "lcaA", etc.

.. code-block:: python

    >>> # Check the TEA results for `sysA`
    >>> pou.teaA.show() # doctest: +SKIP
	TEA: sysA
	NPV  : -404 USD at 5.0% discount rate
    >>> # There are also handy functions to allow you quickly see important results
    >>> pou.print_summaries((pou.sysA, pou.sysB, pou.sysC, pou.sysD)) # doctest: +SKIP

	---------- Summary for sysA ----------


	Total cost: 0.09 USD/cap/yr.

	Total GWP: 0.13 kg CO2-eq/cap/year.

	---------- Summary for sysB ----------


	Total cost: 0.74 USD/cap/yr.

	Total GWP: 0.07 kg CO2-eq/cap/year.

	---------- Summary for sysC ----------


	Total cost: 5.50 USD/cap/yr.

	Total GWP: 2.74 kg CO2-eq/cap/year.

	---------- Summary for sysD ----------


	Total cost: 14.36 USD/cap/yr.

	Total GWP: 1.16 kg CO2-eq/cap/year.
    >>> # You can also save system reports (will be saved to your current working directory)
    >>> # pou.sysA.save_report('sysA.xlsx')


Uncertainty and sensitivity analyses
------------------------------------
Uncertainty and sensitivity analyses are performed through :class:`Model` objects.


References
----------
.. [1] Elijah et al., Assessing the relative sustainability of point-of-use water disinfection technologies for off-grid communities. In Prep._