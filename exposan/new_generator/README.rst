============
NEWgenerator
============

Summary
-------
Systems included in this module are based on the work described in Watabe et al., [1]_ for the NEWgenerator system designed by researchers at the University of South Florida [2]_. Please note that because detailed information of the unit operations is protected by a non-disclosure agreement, only the system configuration is available on the public repository. Please contact the corresponding authors of the referenced papers if further information is needed.

Two system configurations are included in this module.

.. figure:: ./figures/sysA.png

    *Scenario A (sysA): original NEWgenerator with front-end and photovoltaic electricity*

.. figure:: ./figures/sysB.png

    *Scenario B (sysB): original NEWgenerator with front-end and grid-tied electricity*


Loading systems
---------------
.. code-block:: python

    >>> # Import and load the system
    >>> from exposan import new_generator as ng
    >>> ng.load()
    >>> # Quick look at the systems
    >>> ng.sysA.show() # doctest: +ELLIPSIS
    System: sysA
    ins...
    [0] toilet_paper
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Tissue  0.000282
    [1] flushing_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  0.0231
    [2] cleansing_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow: 0
    [3] desiccant
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): WoodAsh  0.00633
    [4] air
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  0.00138
                        N2  0.0052
    [5] LPG
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): LPG  0.0136
    [6] zeolite
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Zeolite  2.81e-05
    [7] GAC
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): GAC  0.00152
    [8] NaCl1
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NaCl  0.000234
    [9] NaOH
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): SodiumHydroxide  1.71e-05
    [10] NaCl
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NaCl  8.18e-05
    outs...
    [0] A11_in
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NH3     0.000163
                        NonNH3  0.000219
                        P       7.69e-05
                        K       0.000185
                        H2O     2.07
    [1] CH4
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CH4         0.000142
                        SolubleCH4  0.000343
    [2] N2O
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): N2O  1.11e-08
    [3] sol_N
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NonNH3  0.000173
    [4] sol_P
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): P  5.14e-05
    [5] sol_K
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): K  5.42e-05
    [6] A_sol_non_fertilizers
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NH3         0.000471
                        Mg          7.71e-05
                        Ca          0.000227
                        H2O         0.00668
                        OtherSS     0.395
                        Tissue      0.0282
                        WoodAsh     0.633
                        ...
    [7] used
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  0.00138
                        CO2  0.000691
    [8] wasted
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CH4  0.000565
    [9] spent_zeolite
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NH3      0.000277
                        Zeolite  2.81e-05
    [10] spent_GAC
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): GAC  0.00152
    [11] conc_NH3
        phase: 's', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NH3  0.000325
    >>> # You can look at individual units
    >>> ng.A3.show() # doctest: +ELLIPSIS
    NEWgeneratorAnMBR: A3
    ins...
    [0] mixed_waste  from  MURT-A2
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): NH3      21.5
                     NonNH3   6.51
                     P        4.26
                     K        9.37
                     Mg       1.88
                     Ca       9.08
                     H2O      4.81e+04
                     OtherSS  390
                     Tissue   28.2
                     WoodAsh  633
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 3171.6 mg/L
         TN         : 564.9 mg/L
         TP         : 85.9 mg/L
         TK         : 188.9 mg/L
    outs...
    [0] A3_treated  to  NEWgeneratorIonExchange-A5
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): NH3         13
                     NonNH3      3.94
                     P           2.67
                     K           7.25
                     H2O         3.72e+04
                     SolubleCH4  5.51
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 570.9 mg/L
         TC         : 110.3 mg/L
         TOC        : 110.3 mg/L
         TN         : 453.5 mg/L
         TP         : 71.4 mg/L
         TK         : 193.9 mg/L
    [1] sludge  to  SludgePasteurization-A4
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (g/hr): NH3         8.02
                     NonNH3      2.42
                     P           1.59
                     K           2.12
                     Mg          1.88
                     Ca          9.08
                     H2O         1.09e+04
                     OtherSS     395
                     Tissue      28.2
                     WoodAsh     633
                     SolubleCH4  1.61
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 570.9 mg/L
         TC         : 98.4 mg/L
         TOC        : 98.4 mg/L
         TN         : 851.5 mg/L
         TP         : 129.7 mg/L
         TK         : 173.0 mg/L
    [2] biogas  to  SludgePasteurization-A4
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (g/hr): CH4  22.4
        WasteStream-specific properties: None for non-liquid waste streams
    [3] A3_CH4  to  Mixer-A12
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow: 0
        WasteStream-specific properties: None for non-liquid waste streams
    [4] A3_N20  to  Mixer-A13
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (g/hr): N2O  0.000491
        WasteStream-specific properties: None for non-liquid waste streams


Techno-economic analysis and life cycle assessment (TEA/LCA)
------------------------------------------------------------
TEA and LCA are performed through :class:`TEA` and :class:`LCA` objects and named as "teaA", "lcaA", etc.

.. code-block:: python

    >>> # Check the TEA results for `sysA`
    >>> ng.teaA.show() # doctest: +SKIP
    SimpleTEA: sysA
    NPV  : -169,400 USD at 5.0% discount rate
    >>> # There are also handy functions to allow you quickly see important results
    >>> ng.print_summaries(ng.sysA) # doctest: +SKIP

    ---------- Summary for sysA ----------


    Total N recovery: 52.5 %.

    Total P recovery: 37.3 %.

    Total K recovery: 22.6 %.

    Total cost: 109.77 USD/cap/yr.

    Net emission: 80.34 kg CO2-eq/cap/yr.

    Net ecosystems damage: 0.19 points/cap/yr.

    Net health damage: 12.39 points/cap/yr.

    Net resources damage: 9.75 points/cap/yr.
    >>> # You can also save system reports
    >>> # ng.sysA.save_report('full_file_path.xlsx')


Uncertainty and sensitivity analyses
------------------------------------
Uncertainty and sensitivity analyses are performed through :class:`Model` objects. Two models (``modelA`` and ``modelB`` for ``sysA`` and ``sysB``, respectively) have been added with all the uncertainty parameters and ranges used in Watabe et al. [1]_.

You can make changes (e.g., add or remove parameters, change uncertainty ranges) to these models or construct your own ones.

.. code-block:: python

    >>> # Run the default model for `sysA`
    >>> modelA = ng.create_model('A')
    >>> # Try use larger samples, here is just to get a quick demo result
    >>> # The result will be automatically saved under the results folder
    >>> ng.run_uncertainty(modelA, N=10) # doctest: +ELLIPSIS
    function `run_uncertainty`
    Total ...


Contextual analysis
-------------------
The ``country_specific.py`` module contains scripts to evaluate the system with contextual parameters that are specific to the deployed country.


References
----------
.. [1] Watabe et al., Prioritizing Research, Development, and Deployment Opportunities for the NEWgeneratorTM Non-Sewered Sanitation System through Quantitative Sustainable Design. In Prep.
.. [2] Shyu et al., The NEWgeneratorTM Non-Sewered Sanitation System: Long-Term Field Testing at an Informal Settlement Community in EThekwini Municipality, South Africa. Journal of Environmental Management 2021, 296, 112921. https://doi.org/10.1016/j.jenvman.2021.112921.