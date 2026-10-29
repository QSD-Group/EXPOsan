# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Convert an empirical result to native streams, utilities, equipment and LCA."""
from uuid import uuid4


def to_qsdsan(result, factors, *, ID=None, avoided_disposal_kg_co2e_t=0.):
    """Return ``(system, lca, streams)`` without changing the global thermo.

    This adapter is a fixed design snapshot. Rebuild after changing inputs.
    ``result.tea()`` owns the screening cash flows, including land; the native
    unit exposes installed costs for use with a project-specific BioSTEAM TEA.
    Impact quantities include the lifetime exactly once. No existing registries
    are cleared. All CFs follow Result.lca's explicit units and sign convention.
    """
    import thermosteam as tmo
    import biosteam as bst
    import qsdsan as qs
    result.lca(factors, avoided_disposal_kg_co2e_t=avoided_disposal_kg_co2e_t)
    tag = ID or 'alternative_'+uuid4().hex[:10]
    chemicals = [tmo.Chemical(i) for i in ('Water','CH4','CO2','O2','N2O','NH3','N2')]
    for name in ('VS','Ash','DryMatter'):
        chemical = tmo.Chemical.blank(name, phase='s', MW=1.)
        chemical.default()
        chemicals.append(chemical)
    chemicals = tmo.Chemicals(chemicals)
    chemicals.compile()
    thermo = tmo.Thermo(chemicals)
    maps = dict(result.streams_kg_h)
    # Raw biogas is an intermediate, not a second final product.
    if 'raw_biogas' in maps:
        maps.pop('raw_biogas')
        cake = maps.pop('digestate_cake')
        maps['digestate_product'] = {k:v*result.feedstock.marketable_fraction for k,v in cake.items()}
        maps['residue'] = {k:v*(1-result.feedstock.marketable_fraction) for k,v in cake.items()}
        # Combustion oxygen and water close the aggregate AD boundary.
        raw_ch4 = result.streams_kg_h['raw_biogas']['CH4']
        emitted = maps['methane_emission']['CH4']
        rng = maps['rng']['CH4']
        burned = raw_ch4-emitted-rng
        maps['oxygen'] = {'O2':burned*(44.01+2*18.01528-16.043)/16.043}
        maps['combustion_water'] = {'Water':burned*2*18.01528/16.043}
    else:
        compost = maps.pop('compost')
        maps['compost_product'] = {k:v*result.feedstock.marketable_fraction for k,v in compost.items()}
        maps['residue'] = {k:v*(1-result.feedstock.marketable_fraction) for k,v in compost.items()}
    input_names = {'feedstock','makeup_water','bulking_agent','oxygen','inoculum'}
    streams = {}
    for name, composition in maps.items():
        phase = 'g' if name in ('oxygen','rng','methane_emission','biogenic_co2','emissions','combustion_water') else 'l'
        streams[name] = qs.SanStream(tag+'_'+name, thermo=thermo, phase=phase,
                                     units='kg/hr', **composition)
    ins = [s for k,s in streams.items() if k in input_names]
    outs = [s for k,s in streams.items() if k not in input_names]

    class EmpiricalProcess(qs.SanUnit):
        _N_ins = len(ins)
        _N_outs = len(outs)
        _units = {}

        def _run(self):
            # A fixed-composition design must be rebuilt when its inlet changes.
            # Do not overwrite inlet edits and silently simulate another case.
            for stream in self.ins:
                name = stream.ID[len(tag)+1:]
                for chemical in chemicals.IDs:
                    if abs(stream.imass[chemical]-maps[name].get(chemical,0.))>1e-7:
                        raise ValueError('Changed inlet: rebuild the composition-dependent alternative system')
            for stream in self.outs:
                name = stream.ID[len(tag)+1:]
                stream.empty()
                for chemical, mass in maps[name].items():
                    stream.imass[chemical] = mass

        def _design(self):
            self.design_results.update({k:v for k,v in result.design.items() if isinstance(v,(int,float))})
            for item in result.equipment:
                self.design_results[item['name']] = item['size']
                self._units[item['name']] = {
                    'wet t/day':'tonne/day','wet mixture t/year at 330 days':'tonne/year',
                    'Nm3/h':'m3/h','Nm3 raw biogas/h':'m3/h','site':''
                    }.get(item['size_unit'],item['size_unit'])

        def _cost(self):
            self.installed_costs.clear()
            self.cost_basis_audit = []
            for item in result.equipment:
                name = item['name']
                if 'reference_installed_usd' in item:
                    cost = (item['reference_installed_usd']*
                        (self.design_results[name]/item['reference_size'])**item['exponent']*
                        result.economics.target_cepci/item['reference_cepci'])
                    basis = 'Scaled bundled installed capital; purchase/BM split unavailable'
                else:
                    cost = item['installed_usd']
                    basis = 'Site scenario in target-year USD; purchase/BM split unavailable'
                # Installed-only correlations are not FOB equipment quotations.
                # BioSTEAM retains these entries when no purchase costs are supplied.
                self.installed_costs[name] = cost
                self.cost_basis_audit.append(dict(**item,native_installed_usd=cost,
                    purchase_cost_usd=None,bare_module_factor=None,cost_basis=basis))
            inv, h = result.inventory_annual, result.hours
            self.power_utility.consumption = inv.get('electricity_import_kwh',0.)/h
            self.power_utility.production = -inv.get('electricity_export_kwh',0.)/h
            # External heat supply; CHP's internally recovered heat is already netted.
            heat = inv.get('heat_import_mj',0.)/h*1000
            if heat:
                self.add_heat_utility(heat, 298.15, 333.15,heat_transfer_efficiency=1.)
            export = -inv.get('heat_export_mj',0.)/h*1000
            if export:
                utility = self.add_heat_utility(export,298.15,333.15,heat_transfer_efficiency=1.)
                utility.reverse()

    unit = EmpiricalProcess(tag+'_process', ins=ins, outs=outs, thermo=thermo,init_with='SanStream')
    system = bst.System(tag, path=(unit,))
    system.operating_hours = result.hours
    system.simulate()
    expected = sum(item['installed_usd'] for item in result.equipment)
    if abs(system.installed_equipment_cost-expected)>max(1e-6,expected*1e-9):
        raise ValueError('Native installed costs do not reconcile with the source correlations')
    if abs(unit.mass_balance_error())>max(1e-6,unit.F_mass_in*1e-8):
        raise ValueError('Alternative native boundary fails total mass closure')
    indicator = qs.ImpactIndicator(tag+'_GWP', unit='kg CO2-eq')
    lca = qs.LCA(system, lifetime=result.economics.lifetime, indicators=(indicator,),
                 uptime_ratio=result.economics.operating_days/365, simulate_system=False)
    # Explicit inventory items retain diesel, disposal, water and signed exports,
    # alongside utilities, without automatically adding a second utility impact.
    quantities = dict(result.inventory_annual)
    quantities['avoided_disposal_t'] = -result.annual_tonnes
    all_factors = dict(factors, avoided_disposal_t=avoided_disposal_kg_co2e_t)
    for name, quantity in quantities.items():
        if not quantity:
            continue
        activity_unit = ('kWh' if name.endswith('_kwh') else
                         'MJ' if name.endswith('_mj') else
                         'L' if name.endswith('_l') else
                         'tonne' if name.endswith('_t') else 'kg')
        item = qs.ImpactItem(tag+'_'+name, functional_unit=activity_unit,
                             **{indicator.ID:all_factors[name]})
        lca.add_other_item(item, quantity*result.economics.lifetime)
    return system, lca, streams
