#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Yalin Li <mailto.yalin.li@gmail.com>    

    Jianan Feng <jiananf2@illinois.edu>
    
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

import biosteam as bst, thermosteam as tmo

__all__ = ('SAF_TEA', 'create_tea',)

class SAF_TEA(bst.TEA):
    
    __slots__ = ('OSBL_units', 'warehouse', 'site_development',
                 'additional_piping', 'proratable_costs', 'field_expenses',
                 'construction', 'contingency', 'other_indirect_costs', 
                 'labor_cost', 'labor_burden', 'property_insurance',
                 'maintenance', '_ISBL_DPI_cached', '_FCI_cached',
                 '_utility_cost_cached', '_steam_power_depreciation',
                 '_steam_power_depreciation_array',
                 'boiler_turbogenerator')
    
    def __init__(self, system, IRR, duration, depreciation, income_tax,
                 operating_days, lang_factor, construction_schedule,
                 startup_months, startup_FOCfrac, startup_VOCfrac,
                 startup_salesfrac, WC_over_FCI,  finance_interest,
                 finance_years, finance_fraction, OSBL_units, warehouse,
                 site_development, additional_piping, proratable_costs,
                 field_expenses, construction, contingency,
                 other_indirect_costs, labor_cost, labor_burden,
                 property_insurance, maintenance, steam_power_depreciation,
                 boiler_turbogenerator):
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
    
        
    def _DPI(self, installed_equipment_cost):
        '''Direct permanent investment (total installed cost) considering additional factors (e.g., buildings).'''
        factors = self.warehouse + self.site_development + self.additional_piping
        return installed_equipment_cost + self.ISBL_installed_equipment_cost*factors
    
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


def create_tea(sys, OSBL_units=None, cls=None, IRR_value=0.03, income_tax_value=0.21, finance_interest_value=0.03, labor_cost_value=1e6):
    if OSBL_units is None: OSBL_units = bst.get_OSBL(sys.cost_units)
    try:
        BT = tmo.utils.get_instance(OSBL_units, (bst.BoilerTurbogenerator, bst.Boiler))
    except:
        BT = None
    if cls is None: cls = SAF_TEA
    tea = cls(
        system=sys, 
        IRR=IRR_value, # use 0%-3%-5% triangular distribution for waste management, and 5%-10%-15% triangular distribution for biofuel production
        duration=(2022, 2052), # Jones et al. 2014
        depreciation='MACRS7', # Jones et al. 2014
        income_tax=income_tax_value, # Davis et al. 2018
        operating_days=sys.operating_hours/24, # Jones et al. 2014
        lang_factor=None, # related to expansion, not needed here
        construction_schedule=(0.08, 0.60, 0.32), # Jones et al. 2014
        startup_months=6, # Jones et al. 2014
        startup_FOCfrac=1, # Davis et al. 2018
        startup_salesfrac=0.5, # Davis et al. 2018
        startup_VOCfrac=0.75, # Davis et al. 2018
        WC_over_FCI=0.05, # Jones et al. 2014
        finance_interest=finance_interest_value, # use 3% for waste management, use 8% for biofuel
        finance_years=10, # Jones et al. 2014
        finance_fraction=0.6, # debt: Jones et al. 2014
        OSBL_units=OSBL_units,
        warehouse=0.04, # Knorr et al. 2013
        site_development=0.09, # Knorr et al. 2013
        additional_piping=0.045, # Knorr et al. 2013
        proratable_costs=0.10, # Knorr et al. 2013
        field_expenses=0.10, # Knorr et al. 2013
        construction=0.20, # Knorr et al. 2013
        contingency=0.10, # Knorr et al. 2013
        other_indirect_costs=0.10, # Knorr et al. 2013
        labor_cost=labor_cost_value, # use default value
        labor_burden=0.90, # Jones et al. 2014 & Davis et al. 2018
        property_insurance=0.007, # Jones et al. 2014 & Knorr et al. 2013
        maintenance=0.03, # Jones et al. 2014 & Knorr et al. 2013
        steam_power_depreciation='MACRS20',
        boiler_turbogenerator=BT)
    return tea