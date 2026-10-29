# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Cornell-corrected composition-based AD with CHP/RNG and digestate handling."""
from math import log
from dataclasses import dataclass
from ._common import Feedstock, Economics, Equipment, Result, FEEDSTOCKS, positive, fraction

NREL = 'https://www.nrel.gov/docs/fy13osti/57082.pdf'
RNG_US_REFERENCE = 'https://www.epa.gov/system/files/documents/2024-01/pdh_full.pdf'
WSU = 'https://csanr.wsu.edu/anaerobic-digestion-systems/enterprise-budget-calculator/'


@dataclass(frozen=True)
class ADConfig:
    conversion_model: str = 'cornell'
    biodegradability: float | None = None
    reactor_conversion: float | None = None
    route: str = 'chp'
    hrt_days: float = 22.
    headspace_fraction: float = .20
    digester_solids: float = .12
    methane_volume_fraction: float = .60
    methane_leakage: float = .02
    combustion_slip: float = .01
    rng_recovery: float = .97
    electrical_efficiency: float = .35
    thermal_efficiency: float = .40
    boiler_efficiency: float = .85
    electricity_kwh_t: float = 35.
    heat_mj_t_slurry: float = 100.
    upgrading_kwh_nm3_biogas: float = .25
    cake_solids: float = .30
    solids_capture: float = .95
    land_acres_at_reference: float = 3.
    rng_usd_mj: float = .012
    # EPA Table 4-10: compression/treatment, 2020 USD per scfm inlet gas.
    # https://www.epa.gov/system/files/documents/2024-01/pdh_full.pdf
    rng_usd_per_scfm_at_1000: float = 8300.
    rng_usd_per_scfm_at_6000: float = 6200.
    interconnection_usd: float = 0.
    scaling_exponent: float = .60

    def __post_init__(self):
        if self.conversion_model not in ('cornell','empirical'):
            raise ValueError('conversion_model must be cornell or empirical')
        for key in ('biodegradability','reactor_conversion'):
            if getattr(self,key) is not None:
                fraction(key,getattr(self,key))
        for key in ('rng_usd_per_scfm_at_1000','rng_usd_per_scfm_at_6000'):
            positive(key,getattr(self,key))
        if self.route not in ('chp', 'rng'):
            raise ValueError('route must be chp or rng')
        for key in ('headspace_fraction', 'methane_leakage', 'combustion_slip',
                    'rng_recovery', 'electrical_efficiency', 'thermal_efficiency',
                    'solids_capture', 'scaling_exponent'):
            fraction(key, getattr(self, key))
        for key in ('digester_solids', 'methane_volume_fraction', 'cake_solids', 'boiler_efficiency'):
            positive(key, getattr(self, key))
            fraction(key, getattr(self, key))
        if self.methane_volume_fraction == 1:
            raise ValueError('Biogas must include CO2')
        if self.combustion_slip == 1:
            raise ValueError('combustion_slip must be less than one')
        if self.electrical_efficiency+self.thermal_efficiency > 1:
            raise ValueError('CHP efficiencies sum to more than one')
        for key in ('hrt_days', 'land_acres_at_reference'):
            positive(key, getattr(self, key))
        for key in ('electricity_kwh_t', 'heat_mj_t_slurry', 'upgrading_kwh_nm3_biogas',
                    'rng_usd_mj', 'interconnection_usd'):
            positive(key, getattr(self, key), True)


# https://www.nrel.gov/docs/fy13osti/57082.pdf Table 13: 7,000 US short t/yr,
# $561/(short t/yr), Jan 2013 report. 2013 annual CEPCI=567.3 is an
# approximate harmonization of heterogeneous project costs in Table 12.
REFERENCE_TPD = 7000 * .90718474 / 330
REFERENCE_USD = 7000 * 561
# Engineering allocation of the bundled installed plant, NOT item quotations.
# Includes CHP, civil work, installation, engineering and startup once.
CAPITAL_SHARES = (
    ('Reception, scale and organics collection pit', .08, 'feed'),
    ('Depackager, grinder and slurry pumps', .10, 'feed'),
    ('Digester concrete vessel', .26, 'volume'),
    ('Insulation, coating and roof/gas holder', .10, 'volume'),
    ('Mixers, heat exchangers, water tank and controls', .10, 'volume'),
    ('Gas piping, H2S scrubbing and safety flare', .08, 'gas'),
    ('Digestate separator and covered storage', .08, 'slurry'),
    ('CHP engine, gas cooling and switchgear', .12, 'power'),
    ('Civil works, electrical installation and startup', .08, 'feed'),
)


def create_system(capacity_tpd=None, *, biogas_nm3_h=None, feedstock='food',
                  config=None, economics=None, characterization=None):
    """Return a screening Result; specify exactly one wet feed or raw biogas capacity.

    Biogas sizing is based on incoming VS and achieved methane yield, without
    changing feedstock composition. Residual liquids require a treatment outlet.
    """
    c, e = config or ADConfig(), economics or Economics.from_htl()
    f = (Feedstock.from_name(feedstock,characterization=characterization)
         if isinstance(feedstock,str) and c.conversion_model == 'cornell' else
         FEEDSTOCKS[feedstock] if isinstance(feedstock, str) else feedstock)
    if (capacity_tpd is None) == (biogas_nm3_h is None):
        raise ValueError('Specify exactly one of capacity_tpd or biogas_nm3_h')
    conversion = None
    if c.conversion_model == 'cornell':
        from ._cornell_ad import predict
        if not f.composition:
            raise ValueError('Cornell AD requires a characterized Feedstock.from_name(...)')
        reactor_conversion = (c.reactor_conversion if c.reactor_conversion is not None else
                              f.ad_reactor_conversion if f.ad_reactor_conversion is not None else 1.)
        conversion = predict(f.composition,
            c.biodegradability if c.biodegradability is not None else f.ad_biodegradability,
            reactor_conversion)
        conversion['biodegradability_provenance'] = (
            dict(status='explicit ADConfig override',source='',basis='explicit fraction')
            if c.biodegradability is not None else dict(f.ad_biodegradability_metadata))
        conversion['yield_basis'] = ('corrected BMP screening; no reactor kinetic derating'
            if c.reactor_conversion is None and f.ad_reactor_conversion is None else
            'explicit reactor conversion applied to corrected BMP')
        gas_per_t = 1000*(conversion['methane_nm3_kg_wet']+conversion['co2_nm3_kg_wet'])
    else:
        gas_per_t = 1000*f.total_solids*f.vs_fraction*f.methane_yield/c.methane_volume_fraction
    positive('feedstock biogas yield', gas_per_t, True)
    if biogas_nm3_h is not None:
        positive('feedstock biogas yield for reverse sizing',gas_per_t)
        positive('biogas_nm3_h', biogas_nm3_h)
        capacity_tpd = biogas_nm3_h*24/gas_per_t
    positive('capacity_tpd', capacity_tpd)
    wet = capacity_tpd*1000/24
    dry, water = wet*f.total_solids, wet*(1-f.total_solids)
    vs = dry*f.vs_fraction
    added_water = max(dry/c.digester_solids-wet, 0.)
    slurry = wet+added_water
    methane_nm3_h = (wet*conversion['methane_nm3_kg_wet'] if conversion else vs*f.methane_yield)
    biogas = capacity_tpd*gas_per_t/24
    # Normal densities from molecular weights / 22.414 Nm3/kmol.
    methane = methane_nm3_h*16.043/22.414
    co2 = (biogas-methane_nm3_h)*44.01/22.414
    carbon_gas = methane*12.011/16.043 + co2*12.011/44.01
    incoming_carbon = wet*conversion['carbon_kg_kg_wet'] if conversion else vs*f.carbon_in_vs
    reaction_water = wet*conversion['reaction_water_kg_kg_wet'] if conversion else 0.
    ammonia = wet*conversion['ammonia_kg_kg_wet'] if conversion else 0.
    if methane+co2 > vs+reaction_water+1e-8 or carbon_gas > incoming_carbon+1e-8:
        raise ValueError('BMP/biogas composition exceeds available VS mass or carbon')
    remaining_dry = dry+reaction_water-methane-co2
    available_water = water+added_water-reaction_water
    if available_water < 0:
        raise ValueError('AD stoichiometry consumes more water than available')
    cake_dry = remaining_dry*c.solids_capture
    cake_water = cake_dry*(1/c.cake_solids-1)
    if cake_water > available_water:
        raise ValueError('Cake moisture requires more water than available; lower cake moisture')
    cake = cake_dry+cake_water
    liquid = slurry-methane-co2-cake
    marketable = cake*f.marketable_fraction
    leak = methane*c.methane_leakage
    captured = methane-leak
    # Methane LHV = 50 MJ/kg, screening approximation.
    heat_demand = slurry/1000*c.heat_mj_t_slurry
    power_demand = capacity_tpd/24*c.electricity_kwh_t
    if c.route == 'chp':
        burned = captured*(1-c.combustion_slip)
        slip, rng = captured*c.combustion_slip, 0.
        # Fixed-efficiency CHP screening, methane LHV = 50 MJ/kg.
        generation = burned*50*c.electrical_efficiency/3.6
        recovered_heat = burned*50*c.thermal_efficiency
        heat_import = max(heat_demand-recovered_heat, 0.)
        heat_export = max(recovered_heat-heat_demand, 0.)
    else:
        # Internal biogas boiler supplies process heat before upgrading.
        boiler_ch4 = min(captured, heat_demand/(50*c.boiler_efficiency*(1-c.combustion_slip)))
        burned = boiler_ch4*(1-c.combustion_slip)
        tail = (captured-boiler_ch4)*(1-c.rng_recovery)
        slip = boiler_ch4*c.combustion_slip+tail
        rng = (captured-boiler_ch4)*c.rng_recovery
        generation, heat_export = 0., 0.
        heat_import = max(heat_demand-burned*50*c.boiler_efficiency, 0.)
        power_demand += biogas*c.upgrading_kwh_nm3_biogas
    net_power = power_demand-generation
    h = e.operating_days*24
    inventory = dict(electricity_import_kwh=max(net_power, 0)*h,
                     electricity_export_kwh=-max(-net_power, 0)*h,
                     heat_import_mj=heat_import*h, heat_export_mj=-heat_export*h,
                     water_kg=added_water*h, methane_kg=(leak+slip)*h,
                     biogenic_co2_kg=(co2+burned*44.01/16.043)*h,
                     residue_kg=(cake-marketable)*h, liquid_kg=liquid*h)
    volume = slurry/1000*24*c.hrt_days*(1+c.headspace_fraction)
    # Fixed reference design, so changing feed properties changes equipment sizes.
    ref_f, ref_c = FEEDSTOCKS['food'], ADConfig()
    ref_wet = REFERENCE_TPD*1000/24
    ref_slurry = ref_wet*ref_f.total_solids/ref_c.digester_solids
    ref_gas = ref_wet*ref_f.total_solids*ref_f.vs_fraction*ref_f.methane_yield/.60
    ref_volume = ref_slurry/1000*24*ref_c.hrt_days*1.2
    ref_power = ref_gas*.60*16.043/22.414*.98*.99*50*.35/3.6
    drivers = dict(feed=(capacity_tpd, REFERENCE_TPD, 'wet t/day'),
                   slurry=(slurry, ref_slurry, 'kg/h'), volume=(volume, ref_volume, 'm3'),
                   gas=(biogas, ref_gas, 'Nm3/h'), power=(generation, ref_power, 'kW'))
    equipment = []
    for name, share, driver in CAPITAL_SHARES:
        if driver == 'power' and c.route == 'rng':
            continue
        size, ref, unit = drivers[driver]
        if size:
            equipment.append(Equipment(name, ref, unit, REFERENCE_USD*share,
                             567.3, exponent=c.scaling_exponent,
                             source=NREL+'; assumed equipment allocation').scale(size, e.target_cepci))
    rng_cost_basis = None
    if c.route == 'rng':
        # Inlet to upgrading excludes leaked methane and methane used by the
        # internal boiler. All modeled CO2 remains in the upgrading inlet.
        upgrading_nm3_h = (captured-boiler_ch4)*22.414/16.043+co2*22.414/44.01
        # Nm3: 0 C, 1 atm. Assume scf: 60 F, 1 atm (dry ideal gas).
        # 1 scfm = 60*0.028316846592*(273.15/288.7055556) Nm3/h.
        nm3_h_per_scfm = 60*.028316846592*273.15/288.7055556
        upgrading_scfm = upgrading_nm3_h/nm3_h_per_scfm
        # Log-log curve through EPA's two endpoints, not the generic 0.6 rule.
        ref_cost = 1000*c.rng_usd_per_scfm_at_1000
        large_cost = 6000*c.rng_usd_per_scfm_at_6000
        exponent = log(large_cost/ref_cost)/log(6.)
        extrapolated = not 1000 <= upgrading_scfm <= 6000
        if extrapolated:
            import warnings
            warnings.warn(f'EPA RNG upgrading cost extrapolated: {upgrading_scfm:g} scfm outside 1000-6000 scfm; AD gas is an LFG proxy', stacklevel=2)
        if upgrading_scfm > 0:
            equipment.append(Equipment('RNG gas compression and treatment',1000.,'scfm upgrading inlet',
                ref_cost,596.2,exponent=exponent,
                source=RNG_US_REFERENCE+' Table 4-10, 2020 USD; LFG proxy for AD gas; excludes pipeline/interconnection').scale(upgrading_scfm,e.target_cepci))
        rng_cost_basis = dict(source=RNG_US_REFERENCE,table='4-10',
            inlet_nm3_h=upgrading_nm3_h,inlet_scfm=upgrading_scfm,
            nm3_h_per_scfm=nm3_h_per_scfm,
            standard_conditions='Nm3 at 273.15 K; scf assumed at 288.7055556 K (60 F); both 1 atm, dry ideal gas',
            cost_2020_usd=ref_cost*(upgrading_scfm/1000)**exponent,
            scaling_exponent=exponent,reference_scfm_range=[1000.,6000.],
            extrapolated=extrapolated,feedstock_proxy='Landfill gas cost proxy applied to AD biogas',
            excluded_scope='Pipeline, interconnection; gas-specific cleanup and boiler scope require validation')
    if c.interconnection_usd:
        equipment.append(dict(name='Site-specific interconnection', size=1., size_unit='site',
                              installed_usd=c.interconnection_usd,
                              source='Scenario cost in target-year USD'))
    streams = dict(feedstock={'Water':water, 'VS':vs, 'Ash':dry-vs},
                   makeup_water={'Water':added_water},
                   raw_biogas={'CH4':methane, 'CO2':co2},
                   digestate_cake={'Water':cake_water, 'DryMatter':cake_dry},
                   centrate={'Water':available_water-cake_water,
                             'DryMatter':remaining_dry-cake_dry},
                   rng={'CH4':rng}, methane_emission={'CH4':leak+slip},
                   biogenic_co2={'CO2':co2+burned*44.01/16.043})
    return Result('AD-'+c.route, f, e, capacity_tpd, equipment, streams, inventory,
                  dict(raw_biogas_nm3_h=biogas, methane_nm3_h=methane_nm3_h,
                       digester_volume_m3=volume, slurry_kg_h=slurry,
                       gross_power_kw=generation, power_demand_kw=power_demand,
                       capital_basis=dict(
                           common_source=NREL,
                           reference_capacity_short_t_year=7000.,
                           reference_bundled_installed_usd=REFERENCE_USD,
                           reference_cepci=567.3,
                           allocation='Assumed equipment shares and scaling drivers; not NREL item quotations',
                           chp_cost_source='12% assumed share of NREL bundled installed capital' if c.route == 'chp' else None,
                           rng_upgrading=rng_cost_basis,
                           native_chp_cost_added=False,
                           status=('NREL-calibrated bundled AD-CHP screening allocation; fixed-efficiency CHP energy model'
                               if c.route == 'chp' else
                               'Hybrid screening estimate: NREL bundle minus assumed 12% CHP share plus EPA compression/treatment proxy; boiler and pipeline/interconnection scope require validation')),
                       chp_model='fixed_efficiency_methane_LHV' if c.route == 'chp' else None,
                       recovered_heat_mj_h=recovered_heat if c.route == 'chp' else 0.,
                       heat_demand_mj_h=heat_demand,
                       retained_carbon_kg_h=incoming_carbon-carbon_gas,
                       conversion_model=c.conversion_model,
                       conversion=conversion,
                       raw_methane_volume_fraction=methane_nm3_h/biogas if biogas else None,
                       retained_ammonia_kg_h=ammonia,
                       marketable_product_kg_year=marketable*h,
                       exported_energy_mj_year=rng*50*h+max(-net_power,0)*3.6*h+heat_export*h),
                  c.land_acres_at_reference*(capacity_tpd/REFERENCE_TPD)**.7,
                  revenues=dict(digestate=marketable*h/1000*f.product_price_usd_t,
                                rng=rng*50*h*c.rng_usd_mj),
                  notes=['Installed equipment allocation is assumed, calibrated to NREL bundled plant cost.',
                         'Named-feed biodegradability defaults are cited literature proxies. Without explicit reactor derating, yields represent corrected BMP screening, not validated full-scale performance.',
                         'Wet AD screening: high-solids/woody feeds require a different reactor/cost basis.',
                         'Digestate storage/application emissions and product-use effects require additional inventory.'])


def adbc_dairy_costs(wet_cow_equivalents, *, target_cepci=564.1,
                     flush_separation=False, codigest_volume_fraction=0.):
    """WSU ADBC v2-48 dairy correlations, separate from the NREL cost basis.

    Returns itemized vessel/infrastructure capital only, excluding CHP/CNG,
    fiber/nutrient/water recovery, land, and taxes. No automatic waste-to-cow mapping.
    Uniform CEPCI escalation approximates the workbook's component indices.
    """
    positive('wet_cow_equivalents', wet_cow_equivalents)
    positive('target_cepci', target_cepci)
    fraction('codigest_volume_fraction', codigest_volume_fraction)
    n = wet_cow_equivalents
    if not 500 <= n <= 10000:
        import warnings
        warnings.warn('Dairy correlation extrapolation; confirm project cost with an engineer', stacklevel=2)
    # https://csanr.wsu.edu/anaerobic-digestion-systems/enterprise-budget-calculator/
    # AnaerobicDigester AN5:AU17, 2013 basis; branch at 2500 WCE.
    data = (
        ('Pre-digester separation/thickening',289512.258663,205553.885700,-6.581411664,27.001937521),
        ('Freight, rentals and mobilization',140551.874860,39007.928206,-40.094296788,.523281874),
        ('Excavation, driveway, building and trenching',-257375.370948,68845.371642,172.933669342,42.445372306),
        ('Concrete vessel',514554.503805,297598.448402,58.647810122,145.430232283),
        ('Insulation and coating',-22644.349315,89094.747039,101.545530822,56.849892280),
        ('Hollow core roof panels',57470.136986,1682.090812,36.105914384,58.421132853),
        ('Piping and H2S scrubbing',416764.243114,90046.964333,-75.210964446,55.475947066),
        ('Heat exchangers, tank, controls and pumps',563612.549080,43065.871137,-67.315030675,140.903640502),
        ('Permitting, engineering, testing and bond',108189.534146,-176291.377205,-22.975646341,90.816718199),
        ('Startup and support',170603.414634,-1786.316144,-46.687134146,22.268758165),
        ('Off-farm organics collection pit',28744.292237,61440.302067,15.784817352,2.706413420),
    )
    output = []
    for i,(name,b1,b2,m1,m2) in enumerate(data):
        cost = b1+m1*n if n < 2500 else b2+m2*n
        if i == 0 and not flush_separation:
            continue
        if i == 10 and not codigest_volume_fraction:
            continue
        if 1 <= i <= 9:
            cost *= 1+codigest_volume_fraction
        if cost < 0:
            raise ValueError(f'ADBC correlation produces negative {name} cost at {n} WCE')
        output.append(dict(name=name, size=n, size_unit='WCE',
                           installed_usd=cost*target_cepci/564.1, source=WSU))
    return output


if __name__ == '__main__':
    print(create_system(50).to_json())
