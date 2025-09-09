#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Jianan Feng <jiananf2@illinois.edu>

    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

References:

(1) Jones, S. B.; Zhu, Y.; Anderson, D. B.; Hallen, R. T.; Elliott, D. C.; 
    Schmidt, A. J.; Albrecht, K. O.; Hart, T. R.; Butcher, M. G.; Drennan, C.; 
    Snowden-Swan, L. J.; Davis, R.; Kinchin, C. 
    Process Design and Economics for the Conversion of Algal Biomass to
    Hydrocarbons: Whole Algae Hydrothermal Liquefaction and Upgrading;
    PNNL--23227, 1126336; 2014; https://doi.org/10.2172/1126336.

(2) Davis, R. E.; Grundl, N. J.; Tao, L.; Biddy, M. J.; Tan, E. C.;
    Beckham, G. T.; Humbird, D.; Thompson, D. N.; Roni, M. S. Process Design
    and Economics for the Conversion of Lignocellulosic Biomass to Hydrocarbon
    Fuels and Coproducts: 2018 Biochemical Design Case Update; Biochemical
    Deconstruction and Conversion of Biomass to Fuels and Products via
    Integrated Biorefinery Pathways; NREL/TP--5100-71949, 1483234;
    2018; p NREL/TP--5100-71949, 1483234. https://doi.org/10.2172/1483234.

(3) Knorr, D.; Lukas, J.; Schoen, P. Production of Advanced Biofuels via
    Liquefaction - Hydrothermal Liquefaction Reactor Design: April 5, 2013;
    NREL/SR-5100-60462, 1111191; 2013; p NREL/SR-5100-60462, 1111191.
    https://doi.org/10.2172/1111191.
'''

from biosteam import TEA
import numpy as np, pandas as pd, thermosteam as tmo, biosteam as bst

__all__ = ('HTL_TEA', 'create_tea',)

class CAPEXTableBuilder:
    __slots__ = ('index', 'data')
    
    def __init__(self):
        self.index = []
        self.data =[]
    
    def entry(self, index: str, cost: list, notes: str = '-'):
        self.index.append(index)
        self.data.append([notes, *cost])

    @property
    def total_costs(self):
        N = len(self.data[0])
        return [sum([i[index] for i in self.data]) for index in range(1, N)]
    
    def table(self, names):
        return pd.DataFrame(self.data, 
                            index=self.index,
                            columns=('Notes', *[i + ' [MM$]' for i in names])
        )


class HTL_TEA(TEA):
    
    __slots__ = ('OSBL_units', 'warehouse', 'site_development',
                 'additional_piping', 'proratable_costs', 'field_expenses',
                 'construction', 'contingency', 'other_indirect_costs', 
                 '_labor_cost', 'labor_burden', 'property_insurance',
                 'maintenance', '_ISBL_DPI_cached', '_FCI_cached',
                 '_utility_cost_cached', '_steam_power_depreciation',
                 '_steam_power_depreciation_array',
                 'boiler_turbogenerator', 'land')
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI,  finance_interest,
                 finance_years, finance_fraction, OSBL_units, warehouse,
                 site_development, additional_piping, proratable_costs,
                 field_expenses, construction, contingency,
                 other_indirect_costs, labor_cost, labor_burden,
                 property_insurance, maintenance, steam_power_depreciation,
                 boiler_turbogenerator, land=0.):
        super().__init__(system, IRR, duration, depreciation, income_tax,
                         operating_days, lang_factor, construction_schedule,
                         startup_months, startup_FOCfrac, startup_VOCfrac,
                         startup_salesfrac, WC_over_FCI,  finance_interest,
                         finance_years, finance_fraction)
        self.OSBL_units = OSBL_units
        self.warehouse = warehouse
        self.site_development = site_development
        self.additional_piping = additional_piping
        self.proratable_costs = proratable_costs
        self.field_expenses = field_expenses
        self.construction = construction
        self.contingency = contingency
        self.other_indirect_costs = other_indirect_costs
        self.labor_cost = labor_cost
        self.labor_burden = labor_burden
        self.property_insurance = property_insurance
        self.maintenance = maintenance
        self.steam_power_depreciation = steam_power_depreciation
        self.boiler_turbogenerator = boiler_turbogenerator
        self.land = land
        
    @property
    def working_capital(self) -> float:
        '''Working capital calculated as the sum of WC_over_FCI*FCI and land.'''
        return self.WC_over_FCI * self.FCI+self.land
    
    @property
    def labor_cost(self):
        if hasattr(self, '_labor_cost'):
            if callable(self._labor_cost): return self._labor_cost()
            return self._labor_cost
        return 0.
    @labor_cost.setter
    def labor_cost(self, i):
        self._labor_cost = i

    @property
    def steam_power_depreciation(self):
        """[str] 'MACRS' + number of years (e.g. 'MACRS7')."""
        return self._steam_power_depreciation
    @steam_power_depreciation.setter
    def steam_power_depreciation(self, depreciation):
        self._steam_power_depreciation_array = self._depreciation_array_from_key(
            self._depreciation_key_from_name(depreciation)
        )
        self._steam_power_depreciation = depreciation
    
    @property
    def ISBL_installed_equipment_cost(self):
        return self.installed_equipment_cost - self.OSBL_installed_equipment_cost
    
    @property
    def OSBL_installed_equipment_cost(self):
        if self.lang_factor:
            raise NotImplementedError('lang factor cannot yet be used')
        elif isinstance(self.system, bst.AgileSystem):
            unit_capital_costs = self.system.unit_capital_costs
            OSBL_units = self.OSBL_units
            return sum([unit_capital_costs[i].installed_cost for i in OSBL_units])
        else:
            return sum([i.installed_cost for i in self.OSBL_units])
    
    def _fill_depreciation_array(self, D, start, years, TDC):
        depreciation_array = self._get_depreciation_array()
        N_depreciation_years = depreciation_array.size
        if N_depreciation_years > years:
            raise RuntimeError('depreciation schedule is longer than plant lifetime')
        system = self.system
        BT = self.boiler_turbogenerator
        if BT is None:
            D[start:start + N_depreciation_years] = TDC * depreciation_array
        else:
            if isinstance(system, bst.AgileSystem): BT = system.unit_capital_costs[BT]
            BT_TDC = BT.installed_cost 
            D[start:start + N_depreciation_years] = (TDC - BT_TDC) * depreciation_array
            
            depreciation_array = self._steam_power_depreciation_array
            N_depreciation_years = depreciation_array.size
            if N_depreciation_years > years:
                raise RuntimeError('steam power depreciation schedule is longer than plant lifetime')
            D[start:start + N_depreciation_years] += BT_TDC * depreciation_array
    
    def _ISBL_DPI(self, installed_equipment_cost):
        """Direct permanent investment of units inside battery limits."""
        if self.lang_factor:
            raise NotImplementedError('lang factor cannot yet be used')
        else:
            factors = self.warehouse + self.site_development + self.additional_piping
            self._ISBL_DPI_cached = self.ISBL_installed_equipment_cost * (1+factors)
        return self._ISBL_DPI_cached
        
    def _DPI(self, installed_equipment_cost):
        return self.OSBL_installed_equipment_cost + self._ISBL_DPI(installed_equipment_cost)
    
    def _indirect_costs(self, TDC):
        return TDC*(self.proratable_costs + self.field_expenses
                    + self.construction + self.contingency
                    + self.other_indirect_costs)
    
    def _FCI(self, TDC):
        self._FCI_cached = FCI = TDC + self._indirect_costs(TDC)
        return FCI
    
    def _FOC(self, FCI):
        return (FCI * self.property_insurance
                + self.ISBL_installed_equipment_cost * self.maintenance
                + self.labor_cost * (1 + self.labor_burden))


def create_tea(sys, OSBL_units=None, cls=None, **kwargs):
    OSBL_units = bst.get_OSBL(sys.cost_units)
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
    except:
        BT = None
    
    if cls is None: cls = HTL_TEA
    
    kwargs_keys = list(kwargs.keys())
    for i in ('IRR_value', 'income_tax_value', 'finance_interest_value', 'labor_cost_value'):
        if i in kwargs_keys: kwargs[i.rstrip('_value')] = kwargs.pop(i)
    
    default_kwargs = {
        'IRR': 0.03, # use 0%-3%-5% triangular distribution for waste management, and 5%-10%-15% triangular distribution for biofuel production
        'duration': (2022, 2052),
        'depreciation': 'MACRS7', # Jones et al. 2014
        'income_tax': 0.275, # Davis et al. 2018
        'operating_days': sys.operating_hours/24, # Jones et al. 2014
        'lang_factor': None, # related to expansion, not needed here
        'construction_schedule': (0.08, 0.60, 0.32), # Jones et al. 2014
        'startup_months': 6, # Jones et al. 2014
        'startup_FOCfrac': 1, # Davis et al. 2018
        'startup_salesfrac': 0.5, # Davis et al. 2018
        'startup_VOCfrac': 0.75, # Davis et al. 2018
        'WC_over_FCI': 0.05, # Jones et al. 2014
        'finance_interest': 0.03, # use 3% for waste management, use 8% for biofuel
        'finance_years': 10, # Jones et al. 2014
        'finance_fraction': 0.6, # debt: Jones et al. 2014
        'OSBL_units': OSBL_units,
        'warehouse': 0.04, # Knorr et al. 2013
        'site_development': 0.09, # Knorr et al. 2013
        'additional_piping': 0.045, # Knorr et al. 2013
        'proratable_costs': 0.10, # Knorr et al. 2013
        'field_expenses': 0.10, # Knorr et al. 2013
        'construction': 0.20, # Knorr et al. 2013
        'contingency': 0.10, # Knorr et al. 2013
        'other_indirect_costs': 0.10, # Knorr et al. 2013
        'labor_cost': 1e6, # use default value
        'labor_burden': 0.90, # Jones et al. 2014 & Davis et al. 2018
        'property_insurance': 0.007, # Jones et al. 2014 & Knorr et al. 2013
        'maintenance': 0.03, # Jones et al. 2014 & Knorr et al. 2013
        'steam_power_depreciation':'MACRS20',
        'boiler_turbogenerator': BT,
        'land':0
        }
    default_kwargs.update(kwargs)
    
    tea = cls(system=sys, **default_kwargs)
    return tea


def capex_table(teas, names=None):
    if isinstance(teas, bst.TEA): teas = [teas]
    capex = CAPEXTableBuilder()
    tea, *_ = teas
    ISBL_installed_equipment_costs = np.array([i.ISBL_installed_equipment_cost / 1e6 for i in teas])
    OSBL_installed_equipment_costs = np.array([i.OSBL_installed_equipment_cost / 1e6 for i in teas])
    capex.entry('ISBL installed equipment cost', ISBL_installed_equipment_costs)
    capex.entry('OSBL installed equipment cost', OSBL_installed_equipment_costs)
    ISBL_factor_entry = lambda name, value: capex.entry(name, ISBL_installed_equipment_costs * value, f"{value:.1%} of ISBL")
    ISBL_factor_entry('Warehouse', tea.warehouse)
    ISBL_factor_entry('Site development', tea.site_development)
    ISBL_factor_entry('Additional piping', tea.additional_piping)
    TDC = np.array(capex.total_costs)
    capex.entry('Total direct cost (TDC)', TDC)
    TDC_factor_entry = lambda name, value: capex.entry(name, TDC * value, f"{value:.1%} of TDC")
    TDC_factor_entry('Proratable costs', tea.proratable_costs)
    TDC_factor_entry('Field expenses', tea.field_expenses)
    TDC_factor_entry('Construction', tea.construction)
    TDC_factor_entry('Contingency', tea.contingency)
    TDC_factor_entry('Other indirect costs (start-up, permits, etc.)', tea.other_indirect_costs)
    TIC = np.array(capex.total_costs) - 2 * TDC
    capex.entry('Total indirect cost', TIC)
    FCI = TDC + TIC
    capex.entry('Fixed capital investment (FCI)', FCI)
    working_capital = FCI * tea.WC_over_FCI
    capex.entry('Working capital', working_capital, f"{tea.WC_over_FCI:.1%} of FCI")
    TCI = FCI + working_capital
    capex.entry('Total capital investment (TCI)', TCI)
    if names is None: names = [i.system.ID for i in teas]
    names = [i for i in names]
    return capex.table(names)

voc_table = bst.report.voc_table

def foc_table(teas, names=None):
    if isinstance(teas, bst.TEA): teas = [teas]
    tea, *_ = teas
    foc = bst.report.FOCTableBuilder()
    ISBL = np.array([i.ISBL_installed_equipment_cost / 1e6 for i in teas])
    FCI = np.array([i.FCI / 1e6 for i in teas])
    labor_cost = np.array([i.labor_cost / 1e6 for i in teas])
    foc.entry('Labor salary', labor_cost)
    foc.entry('Labor burden', tea.labor_burden * labor_cost, '90% of labor salary')
    foc.entry('Maintenance', tea.maintenance * ISBL, f'{tea.maintenance:.1%} of ISBL')
    foc.entry('Property insurance', tea.property_insurance * FCI, f'{tea.property_insurance:.1%} of FCI')
    if names is None: names = [i.system.ID for i in teas]
    names = [i + ' MM$/yr' for i in names]
    return foc.table(names)