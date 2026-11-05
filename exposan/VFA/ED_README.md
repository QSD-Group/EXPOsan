# Isolated VFA electrodialysis prototype

Use `_ed_components.create_ed_components`, `_electrodialysis.Electrodialysis`,
and `ed_example.create_ed_feeds`. The new feed construction uses VFA anion
formulae C3H5O2, C4H7O2, C6H11O2 with separately tracked Na+, K+, Cl-.
The original neutral components and RedoxED flowsheet remain unchanged.
No upstream pH adjustment, dissociation, or reagent demand is inferred.

Run with the local editable repositories on the Python path if not installed:

```sh
PYTHONPATH=/Users/blues/Documents/QSDsan:/Users/blues/Documents/EXPOsan /opt/anaconda3/envs/py312/bin/python -m exposan.VFA.ed_example
PYTHONPATH=/Users/blues/Documents/QSDsan:/Users/blues/Documents/EXPOsan /opt/anaconda3/envs/py312/bin/python -m exposan.VFA.test_electrodialysis
```

## Equations and assumptions

One cell pair, two well-mixed liquid channels, fixed current density i,
width b and depth h. Unknown transfers T are mol/s, positive diluate to
concentrate. The unknown length L gives projected area a=bL, installed CEM+AEM
area 2a, and current I=ia. For each ion, T=aJ; recovery is T_target/N_D,in.
Channel outlets are N_D,in-T and N_C,in+T. Outlet stream properties determine
volume flow and concentrations; water flow stays in its original channel.

Migration is (t_CEM-t_AEM)*i/(zF), using signed charge. Back diffusion is
sum(permeability_m)*(c_surface,C-c_surface,D). The flux is migration minus
back diffusion, already summed over both membranes (no factor of two).
Each membrane's transport fractions sum to one. Net charge transfer must be
zero; independent diffusion coefficients that violate this are rejected.
No artificial sodium correction is made. TODO: coupled multicomponent
zero-current diffusion model for unequal ionic permeabilities.

The provisional limiting-current law is i_lim=k_lim*c_eq,D, with c_eq the
positive charge-equivalent concentration. Surface concentrations are
c_D*(1-i/i_lim) and c_C*(1+i/i_lim). This uses the symmetric-electrolyte
polarization approximation for mixed VFAs. TODO: validate this approximation
and calibrate k_lim for the specified channel flow/geometry. Fixed width,
depth, and feeds define velocity; no hydraulic optimization is performed.

Conductivity=sum(lambda_j*c_j), including all tracked ions. Liquid areal
resistance is h/kappa_D+h/kappa_C. Each membrane's provisional resistance is
r_base+coefficient/mean(c_eq,D,c_eq,C). Voltage=i*r_total and power=voltage*I.
No separate Nernst-layer ohmic term is added to the full-depth channel
resistance. No nonohmic/electrode terms, thermal balance, pressure drop,
pumping, costing, or water transport are included. Output electrical power
is the cell-pair ohmic contribution, not complete apparatus consumption.

All example coefficients, conductivities, geometry, recovery, feed flows,
solute pseudo-properties and membrane resistances are illustrative assumptions.
TODO: replace them with data; the example is verification, not calibration.
Sizing uses SLSQP with a hard nonlinear inequality
`(1 - 1e-6)*i_lim >= i`, evaluated from each trial outlet composition,
even when the surface-polarization correction is disabled. The margin
represents a numerical strict inequality, not an engineering safety factor.
Surface clipping only protects intermediate infeasible trials. Both the
constraint and transport/recovery residuals are checked after solving.

Positive finite parameters, nonnegative feeds, electroneutrality, sub-limiting
current, bounded transfers and actual scaled solver residuals are checked.

## Sources

Equations are adapted to well-mixed channels from WaterTAP ED1D:
https://watertap.readthedocs.io/en/latest/technical_reference/unit_models/electrodialysis_1D.html
https://github.com/watertap-org/watertap/blob/main/watertap/unit_models/electrodialysis_1D.py
Specifically mass-transfer, concentration-polarization and limiting-current
relationships. The inverse-concentration membrane-resistance closure is a
project assumption, not claimed as the WaterTAP constitutive equation.
The multi-ion limitation is discussed in the WaterTAP documentation.
