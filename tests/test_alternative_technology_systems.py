# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""System-level screening balances, financial reconciliation and native LCA."""
import importlib.util
from pathlib import Path
import sys
import unittest
from dataclasses import replace
from math import isclose

# Load the standalone technology package without importing the HTL/ML stack.
root = Path(__file__).resolve().parents[1]
path = root/'exposan/biobinder_ml/deploy/technologies'
spec = importlib.util.spec_from_file_location('technologies', path/'__init__.py',
                                             submodule_search_locations=[str(path)])
tech = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = tech
spec.loader.exec_module(tech)


class TechnologySystems(unittest.TestCase):
    def test_ad_mass_carbon_and_capacity(self):
        for name in tech.FEEDSTOCKS:
            for route in ('chp','rng'):
                with self.subTest(feedstock=name, route=route):
                    c = tech.ADConfig(route=route)
                    r = tech.create_ad_system(50,feedstock=name,config=c)
                    s = r.streams_kg_h
                    incoming = sum(s['feedstock'].values())+sum(s['makeup_water'].values())
                    outgoing = sum(sum(s[k].values()) for k in ('raw_biogas','digestate_cake','centrate'))
                    self.assertAlmostEqual(incoming,outgoing)
                    f = r.feedstock
                    carbon = 50_000/24*f.total_solids*f.vs_fraction*f.carbon_in_vs
                    gas_c = s['raw_biogas']['CH4']*12.011/16.043+s['raw_biogas']['CO2']*12.011/44.01
                    self.assertAlmostEqual(carbon,gas_c+r.design['retained_carbon_kg_h'])
                    sized = tech.create_ad_system(biogas_nm3_h=r.design['raw_biogas_nm3_h'],feedstock=name,config=c)
                    self.assertAlmostEqual(sized.capacity_tpd,50)
                    self.assertGreaterEqual(r.tea()['annual_revenues']['electricity'],0)
                    r.to_json()

    def test_compost_mass_and_carbon(self):
        for name in tech.FEEDSTOCKS:
            r = tech.create_composting_system(50,feedstock=name)
            s = r.streams_kg_h
            incoming = sum(sum(s[k].values()) for k in ('feedstock','bulking_agent','makeup_water','oxygen'))
            outgoing = sum(sum(s[k].values()) for k in ('compost','emissions'))
            self.assertAlmostEqual(incoming,outgoing,places=7)
            self.assertGreater(r.land_acres,0)
            self.assertGreaterEqual(r.design['retained_carbon_kg_h'],0)
            r.to_json()

    def test_reference_cost_and_escalation(self):
        from technologies.ad import REFERENCE_TPD, REFERENCE_USD
        e = tech.Economics(target_cepci=567.3)
        r = tech.create_ad_system(REFERENCE_TPD,economics=e)
        self.assertAlmostEqual(r.tea()['installed_cost'],REFERENCE_USD)
        scaled = tech.create_ad_system(REFERENCE_TPD,economics=replace(e,target_cepci=1134.6))
        self.assertAlmostEqual(scaled.tea()['installed_cost'],2*REFERENCE_USD)
        less_uptime = tech.create_ad_system(REFERENCE_TPD,economics=replace(e,operating_days=200))
        self.assertAlmostEqual(r.tea()['installed_cost'],less_uptime.tea()['installed_cost'])
        for factory in (tech.create_ad_system,tech.create_composting_system):
            small,big = factory(50),factory(100)
            ratio = big.tea()['installed_cost']/small.tea()['installed_cost']
            self.assertTrue(1 < ratio < 2)

    def test_financial_and_lca_reconciliation(self):
        for factory in (tech.create_ad_system,tech.create_composting_system):
            r = factory(50)
            t = r.tea()
            self.assertAlmostEqual(t['annual_revenue']-t['annual_opex'],t['EBITDA'])
            self.assertEqual(t['annual_revenues']['gate_fee'],40*r.annual_tonnes)
            self.assertEqual(t['annual_costs']['feedstock'],0)
            self.assertAlmostEqual(t['cashflows'][-1]-t['cashflows'][1],t['land']+t['working_capital'])
            with self.assertRaises(ValueError):
                r.lca({})
            factors = {k:0. for k in r.inventory_annual}
            factors.update(methane_kg=27.,n2o_kg=273.,electricity_import_kwh=.4,electricity_export_kwh=.4)
            gross = r.lca(factors)
            net = r.lca(factors,avoided_disposal_kg_co2e_t=100)
            self.assertAlmostEqual(gross['annual_kg_co2e']-net['annual_kg_co2e'],100*r.annual_tonnes)
            self.assertAlmostEqual(sum(gross['breakdown'].values()),gross['annual_kg_co2e'])
            zero_quality = factory(50,feedstock=replace(tech.FEEDSTOCKS['food'],marketable_fraction=0))
            self.assertGreater(zero_quality.inventory_annual['residue_kg'],0)
            self.assertEqual(zero_quality.design['marketable_product_kg_year'],0)

    def test_native_results(self):
        try:
            import qsdsan
        except ImportError as error:
            self.skipTest(str(error))
        cases = [tech.create_ad_system(50),tech.create_composting_system(50),
                 tech.create_ad_system(50,config=tech.ADConfig(route='rng')),
                 tech.create_ad_system(50,config=tech.ADConfig(route='rng',heat_mj_t_slurry=10000))]
        for r in cases:
            factors = {k:1. for k in r.inventory_annual}
            sys_,lca,streams = r.to_qsdsan(factors)
            impact = next(iter(lca.get_total_impacts(annual=True).values()))
            self.assertTrue(isclose(impact,r.lca(factors)['annual_kg_co2e'],rel_tol=1e-8))
            self.assertAlmostEqual(sum(s.F_mass for s in sys_.feeds),sum(s.F_mass for s in sys_.products))
            self.assertAlmostEqual(sys_.installed_equipment_cost,r.tea()['installed_cost'])


if __name__ == '__main__':
    unittest.main()
