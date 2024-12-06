==============================
saf: Sustainable Aviation Fuel
==============================

Summary
-------
This module includes a hydrothermal liquefaction (HTL)-based system for the production of sustainable aviation fuel (SAF) and valuable coproducts (hydrogen and fertilizers) from wet organic wastes (manuscript to be submitted [1]_).

Two system configurations are included in the module describing the three scenarios discussed in the manuscript, but the system diagram looks identical (the electrochemical [EC] unit is a placeholder that does nothing in the baseline scenario).

.. figure:: ./readme_figures/sys.svg

    *System diagram.*

Loading systems
---------------
.. code-block:: python

    >>> # Import and load the system
    >>> from exposan import saf
    >>> saf.load(configuration='baseline') # 'baseline', 'EC', 'EC-Future', 'no PSA'
    >>> # Quick look at the systems
    >>> saf.sys.show() # doctest: +ELLIPSIS
    System: sys
    ins...
    [0] makeup_H2  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2  49.7
    [1] recycled_H2  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [2] EC_replacement_surrogate  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [3] H2_HC  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2  23.5
    [4] HCcatalyst_in  
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): HCcatalyst  0.119
    [5] natural_gas  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CH4  25.5
    [6] air  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  60.2
                        N2  226
    [7] H2_HT  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2  26.3
    [8] HTcatalyst_in  
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): HTcatalyst  0.216
    [9] -  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [10] -  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow: 0
    [11] -  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow: 0
    [12] -  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): H2O  223
    [13] feedstock  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): Lipids         9.58
                         Proteins       1.32
                         Carbohydrates  9.79
                         Ash            284
                         H2O            803
    [14] transportation_surrogate  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): Lipids         9.58
                         Proteins       1.32
                         Carbohydrates  9.79
                         Ash            284
                         H2O            803
    [15] feedstock_water  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): H2O  223
    outs...
    [0] process_H2  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2  49.7
    [1] excess_H2  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [2] EC_H2  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [3] recovered_N  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [4] recovered_P  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [5] recovered_K  
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [6] ww_to_disposal  
        phase: 'l', T: 333.05 K, P: 101325 Pa
        flow (kmol/hr): HTLaqueous  48.4
                        1E2PYDIN    0.697
                        ETHYLBEN    0.281
                        4M-PHYNO    0.119
                        4EPHYNOL    0.0214
                        INDOLE      0.000343
                        7MINDOLE    3.97e-05
            1.03e+03
    [7] HCcatalyst_out  
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): HCcatalyst  0.119
    [8] gas_emissions  
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  48.5
                        CO2  62.1
                        N2   226
    [9] solids_to_disposal  
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): HTLchar  284
    [10] HTcatalyst_out  
         phase: 's', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): HTcatalyst  0.216
    [11] s18  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow: 0
    [12] s19  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): H2O  223
    [13] s20  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow: 0
    [14] mixed_fuel  
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): C14H30  3.93
                         C21H44  1.21
                         C8H18   5.09
    >>> # Results
    >>> saf.simulate_and_print(saf.sys) # doctest: +ELLIPSIS
    Fuel properties
    ---------------
    gasoline: 47.86 MJ/kg, 2.77 kg/gal, 212.82 GGE/hr.
    jet: 47.35 MJ/kg, 2.87 kg/gal, 279.48 GGE/hr.
    diesel: 47.10 MJ/kg, 2.99 kg/gal, 130.74 GGE/hr.
    Minimum selling price of all fuel is $3.96/GGE.
    NPV is 1 USD
    AOC is 417,393 USD/yr
    sales is 16,355,081 USD/yr
    net_earnings is 12,590,774 USD/yr
    Global warming potential of all fuel is -5.39 kg CO2e/GGE.

More settings can be changed in the `systems.py` script, the /analyses_ directory includes two sensitivity analyses (with regard to plant size and biocrude yield).


References
----------
.. [1] Si et al., In Prep.