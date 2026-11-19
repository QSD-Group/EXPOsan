# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Native thermochemical systems with common cash flows and auditable markets."""
from dataclasses import dataclass, asdict
from pathlib import Path
from uuid import uuid4
import ast
import json
import warnings
import re
import numpy as np
from ._common import Feedstock, Economics, Result, positive
from ._thermal_basis import ThermalConfig, ProductMarket, default_product_markets
from ._thermal_chemistry import create_thermo, atom_vector, ELEMENTS, ATOMIC, characterize_product


def _native_cost_metadata(unit, name):
    """Retain the installed native correlation's own reference index and sizing."""
    module = 'heat_exchange'
    index = 567.
    reference = 'Seider et al., Product and Process Design Principles, Chapter 16'
    size_basis = dict(unit._units)
    bounds = getattr(unit,'_bounds',{})
    cost_items = getattr(unit,'cost_items',{})
    if name in cost_items:
        item = cost_items[name]
        index = float(item.CE)
        module = 'drying'
        reference = 'Kwiatkowski et al. 2006; native DrumDryer area correlation'
    elif hasattr(unit,'purchase_cost_algorithm'):
        index = float(unit.purchase_cost_algorithm.CE)
        module = 'design_tools/tank_design'
        reference = 'Native selected tank purchase-cost algorithm'
    elif type(unit).__name__ == 'Pump':
        module = '_pump'
    elif type(unit).__name__ in ('SolidsCentrifuge','BoundedCentrifuge'):
        module = 'solids_separation'
    elif 'Compressor' in type(unit).__name__:
        module = 'compressor'
    elif type(unit).__name__ in ('Conversion','GasLiquidSeparator'):
        module = 'design_tools/pressure_vessel'
    return dict(currency='USD',usd_per_currency=1.,reference_cepci=index,
        reference_cost_year=2007 if abs(index-525.2)<.3 else 2013 if abs(index-567)<.01 else None,
        cost_year_note='Native denominator retained exactly; 567 is the rounded 2013 calibration index, not 567.3 annual CEPCI' if index == 567 else 'Native reference index',
        source='https://biosteam.readthedocs.io/en/latest/_modules/biosteam/units/'+module+'.html',
        reference=reference,size_units=size_basis,native_bounds=bounds)


def _cost_owner(unit, name):
    """Resolve auxiliary-prefixed costs to their native sizing/cost owner."""
    if ' - ' in name:
        label, child_name = name.split(' - ',1)
        for auxiliary_name, child in unit.get_auxiliary_units_with_names():
            if auxiliary_name.replace('_',' ').capitalize() == label:
                return _cost_owner(child,child_name)
    return unit,name


def _capital_factors():
    keys = ('warehouse','site_development','additional_piping','proratable_costs',
            'field_expenses','construction','contingency','other_indirect_costs')
    values = {}
    root = Path(__file__).resolve().parents[3]
    for filename, variable in (('htl/_tea.py','default_kwargs'),('saf/_process_settings.py','tea_kwargs')):
        tree = ast.parse((root/filename).read_text(encoding='utf-8-sig'))
        for node in ast.walk(tree):
            if not isinstance(node,ast.Assign) or not any(isinstance(n,ast.Name) and n.id == variable for n in node.targets):
                continue
            v = node.value
            pairs = zip(v.keys,v.values) if isinstance(v,ast.Dict) else ((ast.Constant(k.arg),k.value) for k in v.keywords)
            for key,value in pairs:
                name = ast.literal_eval(key)
                if name in keys:
                    values[name] = ast.literal_eval(value)
    if set(values) != set(keys):
        raise ValueError('HTL capital policy changed; provide an explicitly reviewed capital basis')
    return values


def _market(stream, market, hours, target, year):
    quality = characterize_product(stream,target)
    if market is None:
        return dict(status='UNPRICED', revenue_usd_year=0., marketed_kg_year=0.,
                    reason='No quality-qualified plant-gate sale price supplied',quality=quality)
    if isinstance(market,dict): market = ProductMarket(**market)
    if market.target_year != year:
        raise ValueError('Product price target year must match equipment/operating cost year')
    failures = []
    for key, actual, sense in (
        ('min_purity_mol',quality['purity_mol'],'min'),
        ('max_water_wt',quality['water_wt'],'max'),
        ('max_ash_dry_wt',quality['ash_dry_wt'],'max'),
        ('max_nitrogen_dry_wt',quality['ultimate_dry_wt']['N'],'max'),
        ('max_sulfur_dry_wt',quality['ultimate_dry_wt']['S'],'max'),
        ('max_co_mol',quality['dry_gas_mol'].get('CO',0.),'max'),
        ('min_hhv_mj_kg',quality['hhv_mj_kg'],'min')):
        limit = getattr(market,key)
        if limit is not None and (actual is None or (actual < limit if sense == 'min' else actual > limit)):
            failures.append(key)
    proxy = market.basis != 'qualified_contract'
    if not market.qualified and not proxy:
        failures.append('unverified buyer acceptance/contaminants')
    mass = float(stream.F_mass)*hours*market.saleable_fraction
    quantity = quality['hhv_mj_h']*hours*market.saleable_fraction if market.unit.endswith('MJ-HHV') else mass
    return dict(status='QUALITY_REJECTED' if failures else 'PROXY_SCENARIO' if proxy else 'QUALIFIED_SCENARIO',
                revenue_usd_year=0. if failures else quantity*market.price*market.usd_per_currency*market.year_factor,
                marketed_kg_year=0. if failures else mass, failures=failures,
                qualification='conditional price/marketability proxy; no buyer acceptance inferred' if proxy else 'documentary scenario',
                price_basis=asdict(market), quality=quality)


@dataclass
class ThermalResult(Result):
    """Material/economic snapshot retaining its actual native equipment system."""
    native_system: object = None
    native_streams: dict = None
    config: object = None

    def _check_snapshot(self):
        feed = self.native_streams['feedstock']
        expected = self.capacity_tpd*1000/24
        if abs(feed.F_mass-expected) > max(1e-8,expected*1e-9):
            raise ValueError('Native feed changed: rebuild the Result before economic/inventory reporting')
        if any(abs(float(feed.imass[c.ID])-self.streams_kg_h['feedstock'].get(c.ID,0.)) > max(1e-8,expected*1e-9) for c in feed.chemicals):
            raise ValueError('Native feed composition changed: rebuild the Result before reporting')
        installed = sum(i['installed_usd'] for i in self.equipment)
        if abs(installed-self.native_system.installed_equipment_cost) > max(1e-6,installed*1e-9):
            raise ValueError('Native costs changed: rebuild Result to update its financial snapshot')

    def tea(self, *, native_system=None, return_tea=False):
        """HTL-aligned capital factors and the existing BioSTEAM cash-flow engine."""
        from ._finance import cashflow_analysis
        self._check_snapshot()
        if native_system is not None and native_system is not self.native_system:
            raise ValueError('Use the native system owned by this result')
        e, inv, q = self.economics, self.inventory_annual, self.annual_tonnes
        installed = float(self.native_system.installed_equipment_cost)
        factors = _capital_factors()
        direct_factor = 1+sum(factors[k] for k in ('warehouse','site_development','additional_piping'))
        indirect_factor = 1+sum(factors[k] for k in ('proratable_costs','field_expenses','construction','contingency','other_indirect_costs'))
        direct = installed*direct_factor
        fixed = direct*indirect_factor
        land = self.land_acres*e.land_usd_acre if e.purchase_land else 0.
        costs = dict(maintenance=installed*e.maintenance_fraction,
            insurance=fixed*e.insurance_fraction,labor=e.labor_usd_year*(1+e.labor_burden),
            other=q*e.other_opex_usd_t,feedstock=max(e.feedstock_price_usd_t,0)*q,
            electricity=inv.get('electricity_import_kwh',0)*e.electricity_buy_usd_kwh,
            heat=inv.get('heat_import_mj',0)*e.heat_buy_usd_mj,
            cooling=inv.get('cooling_mj',0)*self.design['cooling_price_usd_mj'],
            water=inv.get('water_kg',0)/1000*e.water_usd_t,
            residue=inv.get('residue_kg',0)/1000*e.residue_disposal_usd_t,
            liquid=inv.get('liquid_kg',0)/1000*e.liquid_disposal_usd_t,
            land_lease=0. if e.purchase_land else self.land_acres*e.land_lease_usd_acre_year,
            **self.extra_costs)
        revenues = dict(gate_fee=max(-e.feedstock_price_usd_t,0)*q,**self.revenues)
        foc = sum(costs[k] for k in ('maintenance','insurance','labor','land_lease'))
        opex,sales = sum(costs.values()),sum(revenues.values())
        finance,tea = cashflow_analysis(e,fixed,land,foc,opex-foc,sales,
                                       native_system=self.native_system,return_tea=True)
        wc = fixed*e.working_capital_fraction
        r,n = e.discount_rate,e.lifetime
        crf = r/(1-(1+r)**-n) if r else 1/n
        annualized = (fixed+land+wc-(land+wc)/(1+r)**n)*crf
        analysis = dict(installed_cost=installed,direct_permanent_investment=direct,
            fixed_capital=fixed,land=land,working_capital=wc,total_investment=fixed+land+wc,
            capital_factors=factors,capital_policy='HTL ISBL direct and indirect factors; no OSBL unit classification in this boundary',
            annual_costs=costs,annual_revenues=revenues,annual_opex=opex,annual_revenue=sales,
            EBITDA=sales-opex,annualized_capital=annualized,
            levelized_net_cost_usd_t=(annualized+opex-sales)/q,
            net_operating_cost_usd_t=(opex-sales)/q,**finance,
            cost_year=self.config.cost_year,target_cepci=e.target_cepci,
            comparison_ready=False,open_gaps=list(self.notes))
        return (analysis,tea) if return_tea else analysis

    def to_json(self):
        self._check_snapshot()
        fields = ('technology','capacity_tpd','equipment','streams_kg_h','inventory_annual',
                  'design','land_acres','revenues','extra_costs','notes')
        data = {name:getattr(self,name) for name in fields}
        data.update(feedstock=asdict(self.feedstock),economics=asdict(self.economics),
                    config=asdict(self.config),tea=self.tea())
        return json.dumps(data,indent=2,allow_nan=False)

    def to_qsdsan(self, factors, *, ID=None, avoided_disposal_kg_co2e_t=0.):
        """Attach QSDsan inventory LCA to the already-built native system."""
        self._check_snapshot()
        self.lca(factors,avoided_disposal_kg_co2e_t=avoided_disposal_kg_co2e_t)
        import qsdsan as qs
        tag = ID or self.native_system.ID
        indicator = qs.ImpactIndicator(tag+'_GWP',unit='kg CO2-eq')
        lca = qs.LCA(self.native_system,lifetime=self.economics.lifetime,indicators=(indicator,),
            uptime_ratio=self.economics.operating_days/365,simulate_system=False)
        quantities = dict(self.inventory_annual,avoided_disposal_t=-self.annual_tonnes)
        cf = dict(factors,avoided_disposal_t=avoided_disposal_kg_co2e_t)
        for key,value in quantities.items():
            if value:
                unit = 'kWh' if key.endswith('_kwh') else 'MJ' if key.endswith('_mj') else 'tonne' if key.endswith('_t') else 'kg'
                item = qs.ImpactItem(tag+'_'+key,functional_unit=unit,**{indicator.ID:cf[key]})
                lca.add_other_item(item,value*self.economics.lifetime)
        return self.native_system,lca,self.native_streams


def create_system(technology, capacity_tpd, *, feedstock='food', config=None,
                  economics=None, characterization=None, markets=None, ID=None,
                  cooling_price_usd_mj=0.):
    """Build a native network; capacity is wet metric t/operating day.

    Residual liquid and solids leave this boundary for paid offsite management,
    consistently with AD/composting. Product-use credits are not presumed.
    """
    import biosteam as bst
    from ._thermal_units import (Conversion, WaterHX, WaterRecuperator, SaltSeparator,
        GasLiquidSeparator, ProductSeparator, BundledCompressor, HotCyclone, OffgasOxidizer,
        BoundedCentrifuge, MinimumReciprocatingCompressor, OilCondenser)
    positive('capacity_tpd',capacity_tpd)
    positive('cooling_price_usd_mj',cooling_price_usd_mj,True)
    c = config or ThermalConfig()
    e = economics or Economics.from_htl()
    if e.engineering_fraction or e.contingency_fraction or e.site_development_usd_acre:
        raise ValueError('Thermochemical TEA already applies HTL capital factors; additional generic capital adders would double count')
    if technology not in ('fast_pyrolysis','htc','htg_methane','htg_hydrogen'):
        raise ValueError('Unknown thermochemical technology')
    is_htg = technology.startswith('htg_')
    mode = 'htg' if is_htg else technology
    allowed_costs = {'pyrolyzer','cyclone'} if mode == 'fast_pyrolysis' else {'reactor','salt_separator','gas_purification'} if is_htg else {'reactor'}
    if set(c.cost_overrides)-allowed_costs:
        raise ValueError(f'Unknown/inapplicable cost overrides: {sorted(set(c.cost_overrides)-allowed_costs)}')
    if technology == 'fast_pyrolysis' and not 673.15 <= c.temperature_k <= 873.15:
        raise ValueError('Fast-pyrolysis scenario temperature must be 400-600 C')
    if technology == 'htc' and (c.temperature_k >= 647.096 or c.pressure_pa < 1e6):
        raise ValueError('HTC requires subcritical hot compressed water')
    if is_htg and (c.temperature_k <= 647.096 or c.pressure_pa <= 22.064e6):
        raise ValueError('HTG factories represent supercritical water service')
    f = Feedstock.from_name(feedstock,characterization=characterization) if isinstance(feedstock,str) else feedstock
    if not isinstance(f,Feedstock):
        raise TypeError('feedstock must be a named feed or Feedstock')
    thermo,composition,products,atoms = create_thermo(f,c,mode)
    if mode == 'htc' and c.pressure_pa <= float(thermo.chemicals.Water.Psat(c.temperature_k)):
        raise ValueError('HTC pressure must keep bulk water below its boiling point')
    tag = ID or technology+'_'+uuid4().hex[:8]
    streams,units = {},[]
    def stream(name,**kwargs):
        result = bst.Stream(tag+'_'+name,thermo=thermo,**kwargs)
        streams[name] = result
        return result
    def unit(cls,name,ins,outs,**kwargs):
        result = cls(tag+'_'+name,ins=ins,outs=outs,thermo=thermo,**kwargs)
        units.append(result)
        return result
    wet = capacity_tpd*1000/24
    feed = stream('feedstock',units='kg/hr',**{k:v*wet for k,v in composition.items()})
    received = stream('received')
    # https://biosteam.readthedocs.io/en/latest/API/units/tank.html : native
    # V=Q*tau/fill sizing. 8 h receiving storage and 80% fill are site scenarios.
    unit(bst.StorageTank,'receiving',feed,received,tau=8.,V_wf=.8,vessel_material='Stainless steel')
    makeup = stream('makeup_water',units='kg/hr',Water=0.)
    h = e.operating_days*24
    outlets = []
    product_names = []
    if mode == 'fast_pyrolysis':
        cake,centrate = stream('dewatered_feed'),stream('centrate')
        solid_ids = tuple(k for k in composition if k != 'Water')
        unit(BoundedCentrifuge,'feed_dewatering',received,(cake,centrate),
            split={k:c.solids_capture for k in solid_ids},solids=solid_ids,
            moisture_content=c.cake_moisture,strict_moisture_content=True)
        dried,drygas,dryemissions = stream('dried_feed'),stream('dryer_exhaust',phase='g'),stream('dryer_emissions',phase='g')
        air,ng = stream('dryer_air',phase='g'),stream('dryer_natural_gas',phase='g')
        unit(bst.DrumDryer,'feed_drying',(cake,air,ng),(dried,drygas,dryemissions),
            split=0.,moisture_content=c.dryer_moisture,utility_agent='Steam',P=101325.,
            gas_composition=(('N2',.79),('O2',.21)))
        preheated = stream('heated_feed')
        unit(WaterHX,'feed_heater',dried,preheated,T=c.temperature_k,rigorous=False)
        hot = stream('reactor_effluent')
        reactor = unit(Conversion,'pyrolyzer',preheated,hot,technology=mode,
                       config=c,composition=composition,products=products)
        hotchar,vapour = stream('hot_char'),stream('pyrolysis_vapour',phase='g')
        unit(HotCyclone,'hot_cyclone',hot,(hotchar,vapour),config=c,dry_feed=dried)
        char = stream('biochar',phase='s')
        unit(WaterHX,'char_cooler',hotchar,char,T=313.15,rigorous=False)
        cooled = stream('cooled_effluent')
        # https://www.psi.ch/en/cpe/hydropilot : cooling before phase separation.
        # 40 C is an engineering condensation/utility target, not a yield fit.
        unit(OilCondenser,'condensation',vapour,cooled,T=313.15,rigorous=False,
             latent_kj_kg=products.get('PyrolysisOil',{}).get('latent_kj_kg',0.))
        gas,oil = stream('raw_gas',phase='g'),stream('bio_oil')
        unit(GasLiquidSeparator,'gas_liquid_separator',cooled,(gas,oil),T=313.15,P=101325.,Q=0.)
        outlets.extend((centrate,drygas,dryemissions,char,oil))
        product_names.extend(('bio_oil','biochar'))
    else:
        dry = wet*(1-composition['Water'])
        makeup.imass['Water'] = max(dry/c.slurry_solids-wet,0.)
        slurry = stream('slurry')
        # https://biosteam.readthedocs.io/en/latest/API/units/tank.html
        # Native stirred-tank sizing;15 min hold/80% fill are engineering scenarios.
        # https://github.com/QSD-Group/QSDsan/blob/main/qsdsan/sanunits/_reactor.py
        # 0.0985 kW/m3 corresponds to the cited0.5 hp/1000 gal agitation heuristic.
        unit(bst.MixTank,'slurry_preparation',(received,makeup),slurry,tau=.25,V_wf=.8,
             kW_per_m3=.0985,vessel_material='Stainless steel')
        pumped = stream('pressurized_slurry')
        # Native pump selection is limited to 20000 ft head for a plunger pump.
        # Series pressure stages avoid silently selecting a centrifugal pump
        # above that envelope. Each stage retains native power and cost sizing.
        stages = max(int(np.ceil((c.pressure_pa-101325.)/1e7)),1)
        inlet = slurry
        for j in range(stages):
            outlet = pumped if j == stages-1 else stream(f'pump_stage_{j+1}_out')
            pressure = 101325.+(c.pressure_pa-101325.)*(j+1)/stages
            unit(bst.Pump,f'high_pressure_pump_{j+1}',inlet,outlet,P=pressure,material='Stainless steel')
            inlet = outlet
        preheated,hot,recovered = stream('preheated_slurry'),stream('reactor_effluent'),stream('heat_recovered_effluent')
        recuperator = unit(WaterRecuperator,'recuperator',(pumped,hot),(preheated,recovered),
            dT=20.,U=.5,ft=1.,material='Stainless steel/stainless steel')
        heated = stream('heated_slurry')
        unit(WaterHX,'trim_heater',preheated,heated,T=c.temperature_k,rigorous=False,
             furnace_pressure=c.pressure_pa/6894.757293)
        if is_htg:
            desalted,brine = stream('desalted_slurry'),stream('salt_brine')
            unit(SaltSeparator,'salt_separator',heated,(desalted,brine),config=c)
            reactor_feed = desalted
            outlets.append(brine)
        else:
            reactor_feed = heated
        reactor = unit(Conversion,'conversion_reactor',reactor_feed,hot,
            technology=mode,config=c,composition=composition,products=products)
        cooled = stream('cooled_effluent')
        unit(WaterHX,'effluent_cooler',recovered,cooled,T=313.15,rigorous=False)
        depressurized = stream('depressurized_effluent')
        unit(bst.IsenthalpicValve,'letdown',cooled,depressurized,P=c.separation_pressure_pa if is_htg else 101325.,vle=False)
        gas,liquid = stream('raw_gas',phase='g'),stream('process_water')
        unit(GasLiquidSeparator,'gas_liquid_separator',depressurized,(gas,liquid),
            T=313.15,P=c.separation_pressure_pa if is_htg else 101325.,Q=0.)
        if is_htg:
            route = technology[4:]
            target = 'CH4' if route == 'methane' else 'H2'
            purified,tail = stream('purified_gas',phase='g'),stream('tail_gas',phase='g')
            unit(ProductSeparator,'purification',gas,(purified,tail),target=target,config=c,route=route)
            final = stream(route,phase='g')
            if c.product_pressure_pa > c.separation_pressure_pa:
                compressor_cls = BundledCompressor if route == 'methane' else MinimumReciprocatingCompressor
                # https://biosteam.readthedocs.io/en/latest/_modules/biosteam/units/compressor.html
                # Reciprocating/electric drive avoids an implicit steam turbine;
                # eta=.75 is a common scenario, not a vendor hydrogen guarantee.
                unit(compressor_cls,'product_compressor',purified,final,P=c.product_pressure_pa,
                    eta=.75,compressor_type='Reciprocating',driver='Electric motor',material='Stainless steel')
            else:
                unit(bst.IsenthalpicValve,'product_letdown',purified,final,P=c.product_pressure_pa,vle=False)
            gas = tail
            outlets.extend((liquid,final))
            product_names.append(route)
        else:
            cake,water = stream('hydrochar_cake'),stream('wastewater')
            unit(BoundedCentrifuge,'hydrochar_dewatering',liquid,(cake,water),
                split={'HydrocharMatrix':c.solids_capture,'Ash':c.solids_capture},
                solids=('HydrocharMatrix','Ash'),moisture_content=c.cake_moisture,
                strict_moisture_content=True)
            final,drygas,dryemissions = stream('hydrochar'),stream('dryer_exhaust',phase='g'),stream('dryer_emissions',phase='g')
            air,ng = stream('dryer_air',phase='g'),stream('dryer_natural_gas',phase='g')
            unit(bst.DrumDryer,'hydrochar_drying',(cake,air,ng),(final,drygas,dryemissions),split=0.,
                moisture_content=c.product_moisture,utility_agent='Steam',P=101325.,
                gas_composition=(('N2',.79),('O2',.21)))
            outlets.extend((water,final,drygas,dryemissions))
            product_names.append('hydrochar')
    stack,air,ng = stream('stack_gas',phase='g'),stream('oxidizer_air',phase='g'),stream('oxidizer_natural_gas',phase='g')
    unit(OffgasOxidizer,'offgas_oxidizer',(gas,air,ng),stack,duty_per_kg=0.)
    outlets.append(stack)
    system = bst.System.from_units(tag,units=units)
    system.operating_hours = h
    old_ce = bst.CE
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            bst.CE = e.target_cepci
            system.simulate()
        finally:
            bst.CE = old_ce
    feeds = tuple(system.feeds)
    boundary_outs = tuple(system.products)
    mass_in = sum(float(s.F_mass) for s in feeds)
    mass_out = sum(float(s.F_mass) for s in boundary_outs)
    def elemental(s):
        return sum((float(s.imol[k.ID])*atom_vector(k) for k in thermo.chemicals),np.zeros(5))*ATOMIC
    element_in = sum((elemental(s) for s in feeds),np.zeros(5))
    element_out = sum((elemental(s) for s in boundary_outs),np.zeros(5))
    residual = element_out-element_in
    ash_residual = sum(float(s.imass['Ash']) for s in boundary_outs)-sum(float(s.imass['Ash']) for s in feeds)
    if abs(mass_out-mass_in) > max(1e-6,mass_in*1e-8) or np.max(np.abs(residual)) > max(1e-6,np.max(element_in)*1e-8) or abs(ash_residual) > 1e-6:
        raise ValueError(f'Native boundary fails mass/element/ash closure: {mass_out-mass_in}, {residual}, {ash_residual}')
    equipment = []
    for u in system.cost_units:
        for name,installed in u.installed_costs.items():
            owner,item_name = _cost_owner(u,name)
            cost_audit = getattr(owner,'cost_audit',[])
            native_basis = None if cost_audit else _native_cost_metadata(owner,item_name)
            equipment.append(dict(name=u.ID+': '+name,unit_class=type(u).__name__,
                sizing_owner_class=type(owner).__name__,
                installed_usd=float(installed),purchase_usd=float(u.purchase_costs[name]) if name in u.purchase_costs else None,
                baseline_purchase_usd=float(u.baseline_purchase_costs[name]) if name in u.baseline_purchase_costs else None,
                bare_module_factor=float(owner.F_BM.get(item_name,1.)),material_factor=float(owner.F_M.get(item_name,1.)),
                pressure_factor=float(owner.F_P.get(item_name,1.)),design=dict(owner.design_results),
                native_correlation=native_basis,
                currency='USD',target_year=c.cost_year,target_cepci=e.target_cepci,
                source_cost_audit=cost_audit,
                source=cost_audit[0]['source'] if cost_audit else native_basis['source']))
    duties = [hu for u in system.units for hu in u.heat_utilities]
    power = sum(float(u.power_utility.rate) for u in system.units)
    markets = default_product_markets(product_names,c.cost_year,f.marketable_fraction) if markets is None else markets
    unknown_markets = set(markets)-set(product_names)
    if unknown_markets:
        raise ValueError(f'Unknown product markets: {sorted(unknown_markets)}')
    target = 'CH4' if technology == 'htg_methane' else 'H2' if technology == 'htg_hydrogen' else None
    market_results = {k:_market(streams[k],markets.get(k),h,target,c.cost_year) for k in product_names}
    for k,market in markets.items():
        settings = asdict(market) if isinstance(market,ProductMarket) else market
        if settings.get('max_sulfur_dry_wt') is not None and not c.ultimate_daf:
            market_results[k].update(status='QUALITY_REJECTED',revenue_usd_year=0.,marketed_kg_year=0.)
            market_results[k].setdefault('failures',[]).append('feed sulfur uncharacterized')
    # Charged residual outlets include unmarketed solids/oil; gases are oxidized.
    liquid_names = ('centrate','wastewater') if mode != 'htg' else ('process_water','salt_brine')
    liquid_mass = sum(float(streams[k].F_mass) for k in liquid_names if k in streams)
    unmarketed = {k:max(float(streams[k].F_mass)-market_results[k]['marketed_kg_year']/h,0.) for k in product_names}
    liquid_mass += sum(v for k,v in unmarketed.items() if k == 'bio_oil')
    residue = sum(v for k,v in unmarketed.items() if k in ('biochar','hydrochar'))
    heating = sum(max(float(hu.duty),0.) for hu in duties)/1000
    cooling = -sum(min(float(hu.duty),0.) for hu in duties)/1000
    inventory = dict(electricity_import_kwh=max(power,0.)*h,heat_import_mj=heating*h,
        cooling_mj=cooling*h,water_kg=float(makeup.F_mass)*h,residue_kg=residue*h,
        liquid_kg=liquid_mass*h,biogenic_co2_kg=float(stack.imass['CO2'])*h,
        methane_emission_kg=float(stack.imass['CH4'])*h,
        sulfur_dioxide_kg=float(stack.imass['SO2'])*h)
    native_warnings = list(dict.fromkeys(re.sub(r'\x1b\[[0-9;]*m','',str(w.message)) for w in caught))
    gaps = ['Conversion is an elemental-retention/stoichiometric scenario, not calibrated waste-specific kinetics or equilibrium.',
        'Oil/char thermophysical properties and HHV are estimates until measured; no molecular speciation or product certification is inferred.',
        'Defined biochemical waste is the receiving boundary; contaminant screening, particle milling, conveyance and site-specific safety equipment remain additional scope.',
        'Wastewater/residual treatment is an offsite boundary with shared scenario fees; environmental factors and permit acceptance are required.',
        'Labor and nonenergy operating prices remain shared scenarios without independently verified plant staffing/contracts.',
        'Offgas destruction uses native combustion and a Seider fired-heater capital proxy; emissions-control vendor design is unresolved.',
        'Cooling water tariff is unresolved unless explicitly supplied; utility inventory is retained.',
        'Default product sales use conditional HHV/literature price proxies, not qualified contracts; markets={} removes all product credits.',
        'No climate ranking is available without complete factors and aligned downstream boundaries.']
    if mode != 'fast_pyrolysis':
        gaps.append('Native vessel cost uses a fixed 15000 psi stress, 0.85 joint efficiency and 1/8 inch corrosion allowance; temperature-dependent ASME design and high-pressure vendor quotes are unresolved.')
    else:
        gaps.append('Pyrolyzer cost year is adopted as 2007 from ORNL attribution; source planning appendix contains a different projected year. Six-tenths scaling is an assumption.')
        gaps.append('Hot cyclone uses specified capture and a throughput cost anchor; particle-size/velocity design, feeder/grinder and solids conveyance require site-specific engineering.')
    if is_htg:
        gaps.extend(['High-pressure hydrocyclone is a solids-removal cost proxy, not a validated salt-precipitation design.',
            'Catalyst/sulfur-guard inventory, replacement and upstream catalyst protection require a waste-specific design; gas impurity splits are specifications.',
            'Purification recoveries/purity are specified; equipment correlations do not prove gas-grid or fuel-cell qualification.'])
    if not c.ultimate_daf:
        gaps.append('Feed sulfur is uncharacterized by the shared biochemical surrogates, not established absent. Supply measured CHONS ultimate_daf to model sulfur.')
    gaps.extend(native_warnings)
    design = dict(conversion_basis=c.conversion_basis,conversion_source=c.conversion_source,
        incoming_wet_kg_h=wet,incoming_dry_kg_h=wet*(1-composition['Water']),
        effective_wet_composition=f.composition,reactor=dict(reactor.design_results),
        product_quality={k:v['quality'] for k,v in market_results.items()},product_markets=market_results,
        raw_gas_quality=characterize_product(streams['raw_gas'],target),
        mass_balance_kg_h=mass_out-mass_in,
        element_balance_kg_h={k:float(v) for k,v in zip(ELEMENTS,residual)},
        ash_balance_kg_h=ash_residual,
        boundary_inputs=[s.ID for s in feeds],boundary_outputs=[s.ID for s in boundary_outs],
        marketable_product_kg_year=sum(v['marketed_kg_year'] for v in market_results.values()),
        cooling_price_usd_mj=cooling_price_usd_mj,
        heat_recovered_mj_h=float(recuperator.total_heat_transfer)/1000 if mode != 'fast_pyrolysis' else 0.,
        property_package=dict(streams='ThermoSTEAM',database='thermo/chemicals',
            feed_surrogates={k:dict(CAS=thermo.chemicals[k].CAS,formula=thermo.chemicals[k].formula,
                resolved_name=getattr(thermo.chemicals[k],'common_name',k)) for k in ('Lipids','Proteins','Carbohydrates')},
            water='CoolProp IAPWS real-fluid correction for hydrothermal duties/density',
            gases='native chemical species; empirical gas/liquid partition',
            oil_char='feed-dependent elemental pseudocomponents; assumed density/Cp; Channiwala-Parikh HHV estimate',
            oil_latent_basis={k:v for k,v in products.get('PyrolysisOil',{}).items() if k not in ('atoms_per_kg_wet','mass_per_kg_wet')}),
        native_warnings=native_warnings,comparison_ready=False)
    # https://biosteam.readthedocs.io/en/latest/API/units/index.html
    # Native fans/mixers/centrifuges account for electricity without routing all
    # dissipated work into stream enthalpy. Retain the real residual separately
    # from the exact material checks; do not replace it with an asserted zero.
    from ._thermal_units import sensible_h
    energy_in = sum(sensible_h(s)+s.Hf for s in feeds)
    energy_out = sum(sensible_h(s)+s.Hf for s in boundary_outs)
    process_heat = sum(float(hu.unit_duty) for hu in duties)
    design['energy_audit_kj_h'] = dict(material_enthalpy_in=energy_in,
        material_enthalpy_out=energy_out,process_heat=process_heat,
        electricity=power*3600.,residual=energy_out-energy_in-process_heat-power*3600.,
        status='unreconciled native mechanical dissipation/property approximations; not an exact enthalpy closure claim')
    stream_maps = {k:{chem.ID:float(s.imass[chem.ID]) for chem in thermo.chemicals if s.imass[chem.ID]}
                   for k,s in streams.items()}
    return ThermalResult(technology,f,e,capacity_tpd,equipment,stream_maps,inventory,design,c.land_acres,
        revenues={k:v['revenue_usd_year'] for k,v in market_results.items()},
        notes=gaps,native_system=system,native_streams=streams,config=c)
