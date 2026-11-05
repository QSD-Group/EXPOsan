"""Focused checks, executable with py312 without requiring pytest."""
import unittest
import numpy as np
from .ed_example import create_example


class TestElectrodialysis(unittest.TestCase):
    def test_conservation_recovery_and_power(self):
        ed = create_example()
        for ID, charge in ed.charges.items():
            self.assertAlmostEqual(ed.components[ID].i_charge * ed.components[ID].MW, charge)
        d, c = ed.ins
        do, co = ed.outs
        np.testing.assert_allclose(d.mol + c.mol, do.mol + co.mol, atol=1e-10)
        for inlet, out in zip(ed.ins, ed.outs):
            self.assertEqual(inlet.imol['H2O'], out.imol['H2O'])
            self.assertAlmostEqual(float(np.asarray(out.imol[ed.ions]) @ ed.z), 0., places=9)
        self.assertAlmostEqual((d.imol['Propionate'] - do.imol['Propionate']) / d.imol['Propionate'], 0.5, places=7)
        r = ed.design_results
        self.assertAlmostEqual(r['Installed membrane area'], 2 * r['Projected area'])
        self.assertAlmostEqual(r['Current'], ed.current_density * r['Projected area'])
        self.assertAlmostEqual(ed.power_utility.rate, r['Voltage'] * r['Current'] / 1000.)
        first = ed.power_utility.rate
        ed.simulate()
        self.assertAlmostEqual(ed.power_utility.rate, first)
        self.assertLess(ed.solver_residual, 1e-7)
        self.assertLessEqual(ed.current_density, (1. - 1e-6) * r['Limiting current density'] + 1e-9)

    def test_back_diffusion_and_polarization(self):
        ed = create_example()
        length = ed.design_results['Length']
        ed.permeabilities[:] = 0.
        ed.simulate()
        self.assertLess(ed.design_results['Length'], length)
        # Ideal migration gives an independent Faraday sizing target.
        expected_transfer = ed.recovery * ed.ins[0].imol[ed.target] / 3.6
        flux = 0.25 * ed.current_density / 96485.33212
        self.assertAlmostEqual(ed.design_results['Projected area'], expected_transfer / flux, places=6)
        _, bulk = ed._bulk_properties(np.array(list(ed.transfers.values())))
        surface, _, ratio = ed._polarization(bulk)
        self.assertGreater(ratio, 0.)
        self.assertTrue(np.all(surface[0] <= bulk[0]))
        self.assertTrue(np.all(surface[1] >= bulk[1]))
        ed.membrane_resistance_coefficient *= 2
        self.assertGreater(ed._resistance(bulk), ed.design_results['Areal resistance'])

    def test_current_constraint_without_surface_correction(self):
        ed = create_example(False)
        ed.polarization = False
        ed.limiting_coefficient = 1e-6
        with self.assertRaises(RuntimeError):
            ed.simulate()

    def test_infeasibility(self):
        ed = create_example(False)
        ed.max_length = 1e-6
        with self.assertRaises(RuntimeError):
            ed.simulate()
        ed = create_example(False)
        ed.limiting_coefficient = 1e-6
        with self.assertRaises(RuntimeError):
            ed.simulate()
        ed = create_example(False)
        ed.permeabilities[:, 0] *= 10
        with self.assertRaisesRegex(ValueError, 'net charge transfer'):
            ed.simulate()
        ed = create_example(False)
        ed.ins[0].imol['Na'] *= 2
        with self.assertRaisesRegex(ValueError, 'charge balanced'):
            ed.simulate()


if __name__ == '__main__':
    unittest.main()
