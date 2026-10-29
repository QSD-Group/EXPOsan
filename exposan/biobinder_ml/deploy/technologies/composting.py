# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    Ali Ahmad <aliahmad1331@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''
"""Composition-based PaCTEA conversion with compost plant sizing and TEA."""
from dataclasses import dataclass, asdict, field
from math import isfinite
from ._common import Feedstock, Economics, Equipment, Result, FEEDSTOCKS, positive, fraction
from ._pactea import PaCTEAConfig, SUBSTRATES, resolve_carbohydrates, simulate as simulate_pactea


@dataclass(frozen=True)
class CompostConfig:
    conversion_model: str = 'pactea'
    pactea: PaCTEAConfig = field(default_factory=PaCTEAConfig)
    bulking_kg_per_kg_feed: float = .125
    bulking_solids: float = .60
    bulking_vs_fraction: float = .95
    bulking_carbon_in_vs: float = .50
    bulking_usd_t: float = 15.
    vs_destroyed: float = .50
    product_moisture: float = .40
    active_days: float = 28.
    curing_days: float = 42.
    bulk_density_kg_m3: float = 750.
    pile_height_m: float = 2.
    usable_site_fraction: float = .35
    minimum_acres: float = 1.
    electricity_kwh_t_mixture: float = 10.
    diesel_l_t_mixture: float = 2.
    # https://www.ipcc-nggip.iges.or.jp/public/2006gl/pdf/5_Volume5/V5_4_Ch4_Bio_Treat.pdf
    # Table 4.1 wet basis: 4 g CH4/kg, 0.24 g N2O/kg treated.
    methane_kg_kg_mixture: float = .004
    n2o_kg_kg_mixture: float = .00024
    reference_installed_usd: float = 8_900_000.
    reference_capacity_t_year: float = 40_000*.90718474
    reference_cepci: float = 596.2
    scaling_exponent: float = .60

    def __post_init__(self):
        if self.conversion_model not in ('pactea','empirical'):
            raise ValueError('conversion_model must be pactea or empirical')
        if not isinstance(self.pactea,PaCTEAConfig):
            raise TypeError('pactea must be a PaCTEAConfig')
        for key in ('bulking_solids', 'bulking_vs_fraction', 'bulking_carbon_in_vs',
                    'vs_destroyed', 'product_moisture', 'usable_site_fraction',
                    'methane_kg_kg_mixture', 'n2o_kg_kg_mixture', 'scaling_exponent'):
            fraction(key, getattr(self, key))
        for key in ('bulk_density_kg_m3', 'pile_height_m', 'usable_site_fraction',
                    'reference_capacity_t_year', 'reference_cepci'):
            positive(key, getattr(self, key))
        if self.product_moisture == 1:
            raise ValueError('Compost moisture must be less than one')
        for key,value in asdict(self).items():
            if key in ('conversion_model','pactea'):
                continue
            positive(key, value, True)


def create_system(capacity_tpd, *, feedstock='food', config=None, economics=None,
                  characterization=None, bulking_composition=None):
    """Return a Result. Capacity excludes separately purchased bulking agent.

    Equipment and land handle the full feed + bulking mixture. PaCTEA predicts
    conversion and nitrogen species; empirical mode is an explicit legacy option.
    """
    positive('capacity_tpd', capacity_tpd)
    c, e = config or CompostConfig(), economics or Economics.from_htl()
    f = (Feedstock.from_name(feedstock,characterization=characterization)
         if isinstance(feedstock,str) and c.conversion_model == 'pactea' else
         FEEDSTOCKS[feedstock] if isinstance(feedstock, str) else feedstock)
    wet = capacity_tpd*1000/24
    bulking = wet*c.bulking_kg_per_kg_feed
    mixture = wet+bulking
    dry = wet*f.total_solids+bulking*c.bulking_solids
    water = mixture-dry
    vs = wet*f.total_solids*f.vs_fraction+bulking*c.bulking_solids*c.bulking_vs_fraction
    carbon = (wet*f.total_solids*f.vs_fraction*f.carbon_in_vs
              +bulking*c.bulking_solids*c.bulking_vs_fraction*c.bulking_carbon_in_vs)
    if c.conversion_model == 'empirical':
        destroyed = vs*c.vs_destroyed
        carbon_lost = carbon*c.vs_destroyed
        methane = mixture*c.methane_kg_kg_mixture
        n2o = mixture*c.n2o_kg_kg_mixture
        ch4_carbon = methane*12.011/16.043
        if ch4_carbon > carbon_lost:
            raise ValueError('Methane emission exceeds degraded carbon; revise emission factor or conversion')
        co2 = (carbon_lost-ch4_carbon)*44.01/12.011
        # Pseudo-substrate mass balance: noncarbon degraded mass forms water;
        # oxygen addition supplies CO2 and N2O oxygen. No microbial kinetics.
        nitrogen_loss = n2o*28.014/44.013
        product_dry = dry-destroyed-nitrogen_loss
        reaction_water = destroyed-carbon_lost-(methane-ch4_carbon)
        if min(product_dry, reaction_water) < 0:
            raise ValueError('Conversion/emission factors violate the dry-mass balance')
        oxygen = co2*(44.01-12.011)/44.01+n2o*(44.013-28.014)/44.013
        product_water = product_dry*c.product_moisture/(1-c.product_moisture)
        makeup = max(product_water-water-reaction_water, 0.)
        evaporated = max(water+reaction_water-product_water, 0.)
    conversion = None
    inoculum = ammonia = nitrogen_gas = 0.
    if c.conversion_model == 'pactea':
        if not f.composition:
            raise ValueError('PaCTEA requires a characterized Feedstock.from_name(...)')
        feed_composition,feed_resolution = resolve_carbohydrates(f.composition,f.carbohydrate_fractions)
        if bulking and bulking_composition is None:
            raise ValueError('Supply measured/cited bulking_composition in wet mass fractions for PaCTEA')
        b,bulking_resolution = resolve_carbohydrates(bulking_composition or {k:0. for k in f.composition})
        if bulking and set((*SUBSTRATES,'Water','Ash'))-set(b):
            raise ValueError('bulking_composition requires Water, Ash and all six PaCTEA substrates')
        if bulking and (any(not isfinite(v) or v < 0 for v in b.values()) or abs(sum(b.values())-1) > 1e-9):
            raise ValueError('bulking_composition must be nonnegative wet fractions summing to one')
        keys = set(feed_composition)|set(b)
        mix = {k:(wet*feed_composition.get(k,0.)+bulking*b.get(k,0.))/mixture for k in keys}
        conversion = simulate_pactea(mix,c.active_days+c.curing_days,c.pactea)
        conversion.update(feedstock_carbohydrate_resolution=feed_resolution,
                          bulking_carbohydrate_resolution=bulking_resolution,
                          input_model_variant=('PaCTEA with lumped-carbohydrate inputs'
                              if 'lumped' in feed_resolution or (bulking and 'lumped' in bulking_resolution)
                              else 'PaCTEA with resolved carbohydrate pools'))
        dry = mixture*(1-mix['Water'])
        water = mixture*mix['Water']
        vs = dry-mixture*mix['Ash']
        inoculum = mixture*conversion['inoculum_kg_kg_mixture']
        methane = mixture*conversion['methane_kg_kg_mixture']
        n2o = mixture*conversion['n2o_kg_kg_mixture']
        co2 = mixture*conversion['co2_kg_kg_mixture']
        ammonia = mixture*conversion['ammonia_kg_kg_mixture']
        nitrogen_gas = mixture*conversion['nitrogen_gas_kg_kg_mixture']
        oxygen = mixture*conversion['oxygen_consumed_kg_kg_mixture']
        product_dry = mixture*conversion['dry_product_kg_kg_mixture']
        product_water = product_dry*c.product_moisture/(1-c.product_moisture)
        # Finish at the specified sale moisture. Water exchange closes the plant
        # boundary; the unadjusted PaCTEA water prediction is reported separately.
        balance_water = mixture+inoculum+oxygen-product_dry-methane-n2o-co2-ammonia-nitrogen_gas
        makeup = max(product_water-balance_water,0.)
        evaporated = max(balance_water-product_water,0.)
        carbon = mixture*conversion['initial_carbon_kg_kg_mixture']
        carbon_lost = carbon-mixture*conversion['carbon_kg_kg_mixture']
        if balance_water < 0:
            raise ValueError('PaCTEA product/emissions violate the plant mass balance')
        conversion['water_balance_adjustment_kg_kg_mixture'] = balance_water/mixture-conversion['model_water_kg_kg_mixture']
    product = product_dry+product_water
    sale = product*f.marketable_fraction
    h = e.operating_days*24
    inventory = dict(electricity_import_kwh=mixture/1000*c.electricity_kwh_t_mixture*h,
                     diesel_l=mixture/1000*c.diesel_l_t_mixture*h,
                     water_kg=makeup*h, bulking_kg=bulking*h,
                     methane_kg=methane*h, n2o_kg=n2o*h, biogenic_co2_kg=co2*h,
                     residue_kg=(product-sale)*h)
    # https://juneau.org/wp-content/uploads/2026/02/Phase-2-Technical-Report_Preliminary-Final-Draft_2026-02-18.pdf
    # Table 2, citing NREL/BR-6A20-81024: ASP 40,000 US t/yr, $8.9M (2020 USD).
    # Shares and exponents below are engineering assumptions; bundled site work
    # is retained. Land acquisition is modeled separately as a scenario addition.
    source = 'https://www.nrel.gov/docs/fy22osti/81024.pdf'
    shares = (
        ('Reception, scale and feedstock mixer', .12),
        ('Grinder and loader', .15),
        ('Aeration blowers, piping and controls', .18),
        ('Composting pad, drainage and leachate collection', .25),
        ('Biofilter and odor control', .10),
        ('Curing area and compost screen', .12),
        ('Electrical, engineering, permitting and startup', .08),
    )
    # Annual design throughput normalized at fixed 330 operating days, so uptime
    # changes utilization without changing equipment size.
    design_t_year = mixture*24/1000*330
    equipment = [Equipment(name, c.reference_capacity_t_year, 'wet mixture t/year at 330 days',
                           c.reference_installed_usd*share, c.reference_cepci,
                           exponent=c.scaling_exponent,
                           source=source+'; assumed equipment allocation').scale(design_t_year,e.target_cepci)
                 for name,share in shares]
    active_volume = mixture*24*c.active_days/c.bulk_density_kg_m3
    curing_volume = product*24*c.curing_days/c.bulk_density_kg_m3
    acres = max(c.minimum_acres, (active_volume+curing_volume)/c.pile_height_m
                /c.usable_site_fraction/4046.8564224)
    streams = dict(feedstock={'Water':wet*(1-f.total_solids),
                             'VS':wet*f.total_solids*f.vs_fraction,
                             'Ash':wet*f.total_solids*(1-f.vs_fraction)},
                   bulking_agent={'Water':bulking*(1-c.bulking_solids),
                                  'VS':bulking*c.bulking_solids*c.bulking_vs_fraction,
                                  'Ash':bulking*c.bulking_solids*(1-c.bulking_vs_fraction)},
                   makeup_water={'Water':makeup}, oxygen={'O2':oxygen},
                   compost={'Water':product_water,'DryMatter':product_dry},
                   emissions={'CH4':methane,'N2O':n2o,'CO2':co2,'Water':evaporated})
    if conversion:
        streams['bulking_agent'] = {'Water':bulking*b['Water'],
            'VS':bulking*(1-b['Water']-b['Ash']),'Ash':bulking*b['Ash']}
        streams['inoculum'] = {'DryMatter':inoculum}
        streams['emissions'].update(NH3=ammonia,N2=nitrogen_gas)
        inventory.update(ammonia_kg=ammonia*h,nitrogen_gas_kg=nitrogen_gas*h,
                         inoculum_kg=inoculum*h)
    return Result('Composting-ASP',f,e,capacity_tpd,equipment,streams,inventory,
                  dict(mixture_tpd=mixture*24/1000, active_volume_m3=active_volume,
                       curing_volume_m3=curing_volume, retained_carbon_kg_h=carbon-carbon_lost,
                       product_solids_fraction=1-c.product_moisture,
                       product_price_usd_t=f.product_price_usd_t,
                       product_price_basis='USD per wet metric tonne of saleable finished compost',
                       conversion_model=c.conversion_model,conversion=conversion,
                       product_quality=(dict(
                           carbon_dry_fraction=conversion['carbon_kg_kg_mixture']/conversion['dry_product_kg_kg_mixture'],
                           nitrogen_dry_fraction=conversion['nitrogen_kg_kg_mixture']/conversion['dry_product_kg_kg_mixture'],
                           organic_nitrogen_kg_h=mixture*conversion['organic_nitrogen_kg_kg_mixture'],
                           ammonium_n_kg_h=mixture*conversion['ammonium_n_kg_kg_mixture'],
                           nitrate_n_kg_h=mixture*conversion['nitrate_n_kg_kg_mixture'],
                           carbon_nitrogen_ratio=(conversion['carbon_kg_kg_mixture']/conversion['nitrogen_kg_kg_mixture']
                               if conversion['nitrogen_kg_kg_mixture'] else None),
                           maturity_certified=False) if conversion else None),
                       marketable_product_kg_year=sale*h),acres,
                  revenues=dict(compost=sale*h/1000*f.product_price_usd_t),
                  extra_costs=dict(bulking_agent=bulking*h/1000*c.bulking_usd_t),
                  notes=['ASP plant-level cost; equipment allocations and land layout are assumptions.',
                         'Lumped carbohydrate inputs do not resolve cellulose, hemicellulose or lignin kinetics; inspect input_model_variant before interpreting compost quality.',
                         'PaCTEA retains published kinetic/stoichiometric defaults; carbon/nitrogen residuals and final-water adjustment are reported in conversion diagnostics.',
                         'Carbon intensity includes direct treatment emissions; soil carbon, fertilizer displacement, transport and construction need explicit additions.',
                         'Marketable fraction represents outlet quality; no pathogen, metals or nutrient certification is predicted.'])


if __name__ == '__main__':
    print(create_system(50).to_json())
