#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Edit the settings below, then run all nine TF cases and their heatmap.

Run: python -m exposan.esap.analyses.tf_validation_script
The spreadsheet, heatmap, and source archive are saved together in OUTPUT_DIR. The original
TF_specific_GPS-X_validation.xlsx remains the source of the GPS-X rows.
"""

from pathlib import Path

from exposan.esap.analyses.tf_validation_sweep import run_validation_sweep

# %% System creation arguments (strength and layer count are swept automatically)
SYSTEM_KWARGS = {
    'gpsx_kinetics': False,
    'gpsx_decay': False,
    'gpsx_aob_proxy': False,
    # Make the currently inherited overrides explicit; this preserves results.
    # Use None here for the GPS-X profile without calibrated ASM overrides.
    # 'asm_kwargs': {'mu_H': 3.3, 'q_PHA': 6.0, 'q_PP': 2.0, 'b_PP': 0.05},
    'tf_kwargs': {
        'N_section': 6,
        'filter_bed_surface': 9464,
        'recirculation_ratio': 1,
        'dynamic_biofilm_thickness': True,
        'initial_biofilm_thickness': 4e-4,
        'detachment_rate': 0.07,
    },
}

# %% Simulation arguments
# TF dynamic_run already applies nutrient/storage atol <= 1e-6 mg/L and
# thickness atol <= 1e-9 m after assembling states, preserving tighter values.
# This standalone system has no AD liquid or gas states to scale separately.
SIMULATE_KWARGS = {
    'method': 'BDF',
    't_span': (0, 300),
    'rtol': 1e-5,
    'atol': 1e-6,
    'max_step': 0.25,
    'stop_at_steady_state': False,
}

# %% Source workbook and output folder
OUTPUT_DIR = Path(
    '/Users/zixuanwang/Library/CloudStorage/GoogleDrive-'
    'wyatt4428@gmail.com/My Drive/4- Research/430- UIUC/434. Projects/'
    '43412. ARPA_E_N_Mg_recovery/process_models/TF_validation/'
)
TEMPLATE_PATH = OUTPUT_DIR / 'TF_specific_GPS-X_validation.xlsx'

# %% Run the sweep, populate the workbook, and plot from the saved values
if __name__ == '__main__':
    validation_run = run_validation_sweep(
        system_kwargs=SYSTEM_KWARGS,
        simulate_kwargs=SIMULATE_KWARGS,
        template_path=TEMPLATE_PATH,
        output_dir=OUTPUT_DIR,
        show_plot=True,
    )
    results = validation_run.results
