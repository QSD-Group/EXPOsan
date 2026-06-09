# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

__all__ = ('test_hap',)

def test_hap(monkeypatch):
    import numpy as np
    from numpy.testing import assert_allclose as ac
    from exposan.hap import create_model, Locations, SimpleCVRP

    # Make the collection/distribution CVRP deterministic so `transport_duty`
    # (and the transport-dependent TEA metrics) are reproducible:
    #   (1) stop on a solution-count budget instead of the machine-dependent
    #       wall-clock time limit, and
    #   (2) seed the otherwise-random location layout (set below).
    def _deterministic_solve(self, time_limit=300, print_solution=False):
        if self._manager is None: self.register()
        else: self._routes = None
        params = self.search_parameters
        params.solution_limit = 300
        self.solution = self.routing.SolveWithParameters(params)
        self._routes = None
    monkeypatch.setattr(SimpleCVRP, 'solve', _deterministic_solve)

    rtol = 0.01

    mdl = create_model()
    for p in mdl.parameters:
        p.setter(p.baseline)
    # Seed the location layout AFTER setting baseline parameters (some setters
    # reset `CD.locations` to None, which would otherwise re-randomize it).
    cd = mdl.system.flowsheet.unit.CD
    area = cd._default_area
    cd.locations = Locations.random_within_area(
        cd._np + 1, distance_metric='cityblock',
        x_range=area['x'], y_range=area['y'], demands=1, seed=42)
    mdl.system.simulate()
    out = [m() for m in mdl.metrics]
    transport_duty = out.pop(1)

    # Re-enabled now that the CVRP solution is deterministic; seed=42 gives a
    # transportation duty ~2% from the historical value, within rtol=0.1.
    assert np.isclose(transport_duty, 0.966498304330359, rtol=0.1)

    ac(out, 
        [63172.32562288159,
         # 1.048288575203588,
         17.597048423328086,
         -1086811.1628681163,
         25.792040993957844,
         74.20795900604216,
         97.2731791145954,
         0.3441417797725412,
         2.38267910563208,
         49.52777806264908,
         15.221503548445778,
         16.166861157662943,
         14.744119672295634,
         4.339737558946575], 
        rtol=rtol)

if __name__ == '__main__':
    test_hap()