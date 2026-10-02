==========================================================================================================
biobinder_ml: Process and Systems Modeling of HTL Fuels and Bioproducts from Wet Organic Wastes
==========================================================================================================

Summary
-------
``biobinder_ml`` is an EXPOsan module for process and systems analysis of hydrothermal liquefaction (HTL) of wet organic wastes. The module links HTL yield prediction with process simulation, techno economic analysis (TEA), life cycle assessment (LCA), and uncertainty analysis for biobinder and biofuel and sustainable aviation fuel (SAF) pathways.

The representative feedstocks used in the system analysis are food waste, green waste, manure, and wastewater sludge. Feedstock compositions are defined in ``_feedstocks.py``. The module supports centralized HTL with centralized upgrading (c HTL) and decentralized HTL with centralized upgrading (d HTL).

Main system scripts
-------------------
``Dist_flex.py``
    Main biobinder and biofuel system model. It creates the HTL system, applies the machine learning yield prediction and distillation specifications, and attaches TEA and LCA objects.

``SAF_FLEX.py``
    Main SAF system model.

``_feedstocks.py``
    Representative wet organic waste compositions used by the system models.

``_rf_htl_predictor.py`` and ``_ml_features.py``
    Utilities for applying the trained random forest HTL yield model to system inputs.

``_distill_utils.py``
    Distillation feasibility and product recovery utilities.

``_components.py``, ``_process_settings.py``, ``_modified_hx.py``, and ``_units.py``
    Components, process settings, heat exchange settings, and unit operations used by the system models.

Techno economic metrics
-----------------------
The system TEA provides internal rate of return (IRR), net present value (NPV), annual revenue, annual operating cost, fixed capital investment, and net earnings. ``Dist_flex.py`` also reports EBITDA, EBITDA margin, and EBITDA relative to total capital investment.

Tax and depreciation are represented through the BioSTEAM/QSDsan TEA cash flow calculation and can be inspected through the TEA cash flow table. Debt schedules, debt service coverage ratio (DSCR), and levered equity returns are not separately modeled in this module.

Life cycle assessment
---------------------
The system models attach QSDsan LCA objects for evaluating environmental impacts. The baseline biobinder system reports global warming potential (GWP). Additional TRACI analyses used in the associated study are analysis workflows rather than required system construction functions.

Uncertainty and sensitivity analysis
------------------------------------
The associated analysis workflows use Monte Carlo simulation to evaluate scenario distributions and use SHAP and rank correlation methods to examine economic drivers. Tipping fee analyses evaluate the waste economics associated with target investment performance. These analyses build on the core system models and are not required to create or simulate a baseline system.

Basic usage
-----------
A biobinder system can be created for a representative feedstock using::

    from exposan.biobinder_ml.Dist_flex import create_system, simulate_and_print

    sys = create_system(
        feedstock_id='food',
        decentralized_HTL=False,
        decentralized_upgrading=False,
        skip_EC=True,
        generate_H2=False,
        EC_config=None,
    )
    simulate_and_print(sys)

Set ``feedstock_id`` to ``'food'``, ``'green'``, ``'manure'``, or ``'sludge'`` for the representative feedstocks. Set ``decentralized_HTL=True`` and ``decentralized_upgrading=False`` for the d HTL configuration.

Testing
-------
The representative feedstock test creates and simulates the baseline biobinder system for food waste, green waste, manure, and wastewater sludge. It checks that the principal TEA and LCA results are finite and prints the calculated IRR, NPV, revenue, EBITDA, and GWP values. This provides a compact regression and smoke test for the feedstock specific system calculations.

Notes
-----
The trained machine learning model and other required data files must be available in the module data or results locations expected by the system scripts. Large Monte Carlo outputs, manuscript figure generation scripts, and postprocessing outputs are not required for basic system simulation.
