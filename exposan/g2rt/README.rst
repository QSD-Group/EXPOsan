========================================================================================
The potential of thermomechanical and thermochemical processes to enable sustainable household sanitation
========================================================================================

Summary
-------
Systems included in this module are based on the work described in Wang et al. [1]_, which follows the HRT design described in [2]_ and [3]_.

Two system configurations are included in this module:

- **Scenario A (sysA):** dry-based pasteurization mechanical dewatering (PMD) HRT
- **Scenario B (sysB):** wet-based supercritical water oxidation (SCWO) HRT

Loading systems
---------------
.. code-block:: python

    >>> # Import and load the module
    >>> from exposan import g2rt
    >>> g2rt.load()

    >>> # System flow diagram
    >>> g2rt.sysA.diagram()

.. figure:: ./readme_figures/PMD.png
   :alt: Process model diagram for sysA

.. code-block:: python

    >>> g2rt.sysB.diagram()

.. figure:: ./readme_figures/SCWO.png
   :alt: Process model diagram for sysB

.. code-block:: python

    >>> # Quick look at the PMD HRT system (sysA)
    >>> g2rt.sysA.show()  # doctest: +ELLIPSIS
    System: sysA
    Highest convergence error among components in recycle
    stream A10-0 after 1 loops:
    - flow rate   5.67e-04 kmol/hr (19%)
    - temperature 0.00e+00 K (0%)
    ins...
    [0] Wood_pellets
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O         0.000127
                        WoodPellet  0.0387
    [1] Air
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): O2  0.000639
                        N2  0.0024
    [2] toilet_paper
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): Tissue  1.56e-06
    [3] tap_water
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): H2O  0.0194
    outs...
    [0] hot_gas
        phase: 'g', T: 375 K, P: 101325 Pa
        flow (kmol/hr): P        1.02e-06
                        K        1.65e-06
                        H2O      0.000255
                        OtherSS  0.00111
                        N2       0.00242
                        CO2      0.000638
    [1] NO
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NO  7.39e-06
    [2] SO2
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): SO2  1.69e-07
    [3] CH4
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): CH4  1.52e-06
    [4] N2O
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): N2O  8.86e-07
    [5] NH3
        phase: 'g', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NH3  5.67e-05
    [6] sol_N
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NonNH3  1.39e-06
    [7] sol_P
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): P  6.85e-06
    [8] sol_K
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): K  1.21e-05
    [9] A_sol_non_fertilizers
        phase: 'l', T: 298.15 K, P: 101325 Pa
        flow (kmol/hr): NH3      9.27e-07
                        Mg       4.6e-06
                        Ca       1.35e-05
                        H2O      7.61e-06
                        OtherSS  0.000196
                        WoodAsh  3.26e-06
    [10] hot_gas
         phase: 'g', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): H2O  0.0074
    [11] A5_reject_discharge
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow: 0
    [12] A6_effluent
         phase: 'l', T: 298.15 K, P: 101325 Pa
         flow: 0
    [13] H2O_vapor1
         phase: 'g', T: 298.15 K, P: 101325 Pa
         flow (kmol/hr): P    3.43e-07
                         K    5.97e-07
                         H2O  0.0324

    >>> # Individual units
    >>> g2rt.A10.show()  # doctest: +ELLIPSIS
    VolumeReductionFilterPress: A10
    ins...
    [0] A9_pasteurized  from  VRpasteurization-A9
    phase: 'l', T: 363.15 K, P: 101325 Pa
    flow (g/hr): NH3      0.0487
                    NonNH3   0.0118
                    P        0.00771
                    K        0.017
                    Mg       0.00339
                    Ca       0.0164
                    H2O      52.2
                    OtherSS  0.712
                    Tissue   1.78
                    sCOD     0.0795
                    xCOD     8.16
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 238.1 mg/L
         TC         : 44362.3 mg/L
         TOC        : 44362.3 mg/L
         TN         : 930.3 mg/L
         TP         : 118.6 mg/L
         TK         : 260.9 mg/L
         TSS        : 152946.1 mg/L
    outs...
    [0] filter_press_liquid  to  G2RTLiquidsTank-A13
    phase: 'l', T: 298.15 K, P: 101325 Pa
    flow (g/hr): NH3      0.0382
                    NonNH3   0.00925
                    P        0.00606
                    K        0.0133
                    Mg       0.00266
                    Ca       0.0129
                    H2O      39.5
                    OtherSS  0.559
                    Tissue   0.178
                    sCOD     0.0625
                    xCOD     0.816
        WasteStream-specific properties:
         pH         : 7.0
         COD        : 21225.2 mg/L
         TC         : 7428.8 mg/L
         TOC        : 7428.8 mg/L
         TN         : 1147.8 mg/L
         TP         : 146.4 mg/L
         TK         : 321.9 mg/L
         TSS        : 24022.4 mg/L
    [1] A10_pressed_solid_cake  to  VRdryingtunnel-A11
    phase: 's', T: 298.15 K, P: 101325 Pa
    flow (g/hr): NH3      0.0104
                    NonNH3   0.00252
                    P        0.00165
                    K        0.00364
                    Mg       0.000728
                    Ca       0.00353
                    H2O      12.6
                    OtherSS  0.153
                    Tissue   1.6
                    sCOD     0.0171
                    xCOD     7.34

Techno-economic analysis and life cycle assessment (TEA/LCA)
------------------------------------------------------------
TEA and LCA are performed through :class:`TEA` and :class:`LCA` objects and named ``teaA`` and ``lcaA`` for PMD HRT, and ``teaB`` and ``lcaB`` for SCWO HRT.

.. code-block:: python

    >>> # Check TEA results
    >>> g2rt.teaA.show()
    TEA: sysA
    ...

    >>> g2rt.teaB.show()
    TEA: sysB
    NPV  : -43,832 USD at 5.0% discount rate

    >>> # Handy functions to quickly see important results
    >>> g2rt.print_summaries(g2rt.sysA)  # doctest: +SKIP

    ---------- Summary for sysA ----------

    Total N recovery: 0.9 %.

    Total P recovery: 79.2 %.

    Total K recovery: 80.1 %.

    Total H2O recovery: 60.0 %.

    Total cost: 534.56 USD/cap/yr.

    Net emission: 396.74 kg CO2-eq/cap/yr.

    Net ecosystems damage: 4.56 points/cap/yr.

    Net health damage: 32.73 points/cap/yr.

    Net resources damage: 19.98 points/cap/yr.

    >>> # You can also save system reports
    >>> # g2rt.sysA.save_report('full_file_path.xlsx')

Uncertainty and sensitivity analyses
------------------------------------
Uncertainty and sensitivity analyses are performed through :class:`Model` objects. ``modelA`` and ``modelB`` (for ``sysA`` and ``sysB``, respectively) have been added with all the uncertainty parameters, distributions, and ranges used in Wang et al. [1]_.

.. code-block:: python

    >>> # Run the default model for sysA
    >>> modelA = g2rt.models.create_model('A')
    >>> # The result will be automatically saved under the results folder
    >>> g2rt.models.run_uncertainty(modelA, N=100)  # doctest: +ELLIPSIS
    function `run_uncertainty`
    Total ...

Contextual analysis
-------------------
The ``location_specific.py`` module contains scripts to evaluate the system with contextual parameters that are specific to the deployed country.

.. code-block:: python

    >>> # Run PMD HRT across 77 countries
    >>> g2rt.run_multiple_cities(system_IDs='A', N=10, note='test')

References
----------
.. [1] Wang et al., *The potential of thermomechanical and thermochemical processes to enable sustainable household sanitation*. In prep.
.. [2] Yee et al., *Volume reduction non-sewered single unit toilet system*.
       `WO2023288326A1 <https://patents.google.com/patent/WO2023288326A1/en?oq=WO2023288326A1>`_
.. [3] Yee et al., *Water oxidation non-sewered single unit toilet system*.
       `WO2023288331A1 <https://patents.google.com/patent/WO2023288331A1/en?oq=WO2023288331A1>`_
