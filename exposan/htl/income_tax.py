#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Jianan Feng <jiananf2@illinois.edu>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
"""

__all__ = ('state_income_tax_rate_2022',
           'state_income_tax_rate_2023',
           'state_income_tax_rate_2024',)

# https://taxfoundation.org/data/all/state/state-corporate-income-tax-rates-brackets-2022/
def state_income_tax_rate_2022(state='IL', sales=20000, net_income=10000):
    if state == 'AL': rate = 0.065
    
    if state == 'AK':
        if net_income <= 25000: rate = 0
        elif net_income <= 49000: rate = (25000*0 + (net_income - 25000)*0.02)/net_income
        elif net_income <= 74000: rate = (25000*0 + (49000 - 25000)*0.02 + (net_income - 49000)*0.03)/net_income
        elif net_income <= 99000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (net_income - 74000)*0.04)/net_income
        elif net_income <= 124000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (net_income - 99000)*0.05)/net_income 
        elif net_income <= 148000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (net_income - 124000)*0.06)/net_income
        elif net_income <= 173000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (net_income - 148000)*0.07)/net_income
        elif net_income <= 198000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (173000 - 148000)*0.07 + (net_income - 173000)*0.08)/net_income
        elif net_income <= 222000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (173000 - 148000)*0.07 + (198000 - 173000)*0.08 + (net_income - 198000)*0.09)/net_income
        else: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (173000 - 148000)*0.07 + (198000 - 173000)*0.08 + (222000 - 198000)*0.09 + (net_income - 222000)*0.094)/net_income
    
    if state == 'AR':
        if net_income <= 3000: rate = 0.01
        elif net_income <= 6000: rate = (3000*0.01 + (net_income - 3000)*0.02)/net_income
        elif net_income <= 11000: rate = (3000*0.01 + (6000 - 3000)*0.02 + (net_income - 6000)*0.03)/net_income
        elif net_income <= 25000: rate = (3000*0.01 + (6000 - 3000)*0.02 + (11000 - 6000)*0.03 + (net_income - 11000)*0.05)/net_income
        else: rate = (3000*0.01 + (6000 - 3000)*0.02 + (11000 - 6000)*0.03 + (25000 - 11000)*0.05 + (net_income - 25000)*0.059)/net_income
    
    if state == 'AZ': rate = 0.049
    
    if state == 'CA': rate = 0.0884
    
    if state == 'CO': rate = 0.0455
    
    if state == 'CT': rate = 0.075
    
    if state == 'DC': rate = 0.0825
    
    # both income tax and gross receipts tax
    # gross receipts tax: use Petroleum wholesaler (https://revenuefiles.delaware.gov/TaxTips/tt-petroleum_dealers2018.pdf)
    if state == 'DE': rate = (net_income*0.087 + sales*0.003983)/net_income
    
    if state == 'FL': rate = 0.055
    
    if state == 'GA': rate = 0.0575
    
    if state == 'HI':
        if net_income <= 25000: rate = 0.044
        elif net_income <= 100000: rate = (25000*0.044 + (net_income - 25000)*0.054)/net_income
        else: rate = (25000*0.044 + (100000 - 25000)*0.054 + (net_income - 100000)*0.064)/net_income
    
    if state == 'IA':
        if net_income <= 100000: rate = 0.055
        elif net_income <=250000: rate = (100000*0.055 + (net_income - 100000)*0.09)/net_income
        else: rate = (100000*0.055 + (250000 - 100000)*0.09 + (net_income - 250000)*0.098)/net_income
    
    if state == 'ID': rate = 0.06
    
    if state == 'IL': rate = 0.095
    
    if state == 'IN': rate = 0.049
    
    if state == 'KS':
        if net_income <= 50000: rate = 0.04
        else: rate = (50000*0.04 + (net_income - 50000)*0.07)/net_income
    
    if state == 'KY': rate = 0.05
    
    if state == 'LA':
        if net_income <= 50000: rate = 0.035
        elif net_income <= 150000: rate = (50000*0.035 + (net_income - 50000)*0.055)/net_income
        else: rate = (50000*0.035 + (150000 - 50000)*0.055 + (net_income - 150000)*0.075)/net_income
    
    if state == 'MA': rate = 0.08
    
    if state == 'MD': rate = 0.0825
    
    if state == 'ME':
        if net_income <= 350000: rate = 0.035
        elif net_income <= 1050000: rate = (350000*0.035 + (net_income - 350000)*0.0793)/net_income
        elif net_income <= 3500000: rate = (350000*0.035 + (1050000 - 350000)*0.0793 + (net_income - 1050000)*0.0833)/net_income
        else: rate = (350000*0.035 + (1050000 - 350000)*0.0793 + (3500000 - 1050000)*0.0833 + (net_income - 3500000)*0.0893)/net_income
    
    if state == 'MI': rate = 0.06
    
    if state == 'MN': rate = 0.098
    
    if state == 'MO': rate = 0.04
    
    if state == 'MS':
        if net_income <= 5000: rate = 0
        elif net_income <= 10000: rate = (5000*0 + (net_income - 5000)*0.04)/net_income
        else: rate = (5000*0 + (10000 - 5000)*0.04 + (net_income - 10000)*0.05)/net_income
    
    if state == 'MT': rate = 0.0675
    
    if state == 'NC': rate = 0.025
    
    if state == 'ND':
        if net_income <= 25000: rate = 0.0141
        elif net_income <= 50000: rate = (25000*0.0141 + (net_income - 25000)*0.0355)/net_income
        else: rate = (25000*0.0141 + (50000 - 25000)*0.0355 + (net_income - 50000)*0.0431)/net_income
    
    if state == 'NE':
        if net_income <= 100000: rate = 0.0558
        else: rate = (100000*0.0558 + (net_income - 100000)*0.075)/net_income
    
    if state == 'NH': rate = 0.076
    
    # note in NJ, the rates apply to the entire net income rather than just income over the threshold
    if state == 'NJ':
        if net_income <= 50000: rate = 0.065
        elif net_income <= 100000: rate = 0.075
        elif net_income <= 1000000: rate = 0.09
        else: rate = 0.115
    
    if state == 'NM':
        if net_income <= 500000: rate = 0.048
        else: rate = (500000*0.048 + (net_income - 500000)*0.059)/net_income
    
    # gross receipts tax: use value for waste management and remediation services (https://tax.nv.gov/uploadedFiles/taxnvgov/Content/Commerce/Commerce_Tax_Instructions.pdf, link invalid)
    if state == 'NV':
        rate = 0 if sales <= 4000000 else (sales - 4000000)*0.00261/net_income
    
    if state == 'NY':
        if net_income <= 5000000: rate = 0.065
        else: rate = (5000000*0.065 + (net_income - 5000000)*0.0725)/net_income
    
    # gross receipts tax (https://taxfoundation.org/data/all/state/state-gross-receipts-taxes-2022/)
    if state == 'OH':
        rate = sales*0.0026/net_income
    
    if state == 'OK': rate = 0.04
    
    # both income tax and gross receipts tax (https://taxfoundation.org/data/all/state/state-gross-receipts-taxes-2022/)
    if state == 'OR':
        if net_income <= 1000000: rate = 0.066
        else: rate = (1000000*0.066 + (net_income - 1000000)*0.076)/net_income
        rate = (net_income*rate + sales*0.0057)/net_income
    
    if state == 'PA': rate = 0.0999
    
    if state == 'RI': rate = 0.07
    
    if state == 'SC': rate = 0.05
    
    # no income tax and no gross receipts tax
    if state == 'SD': rate = 0
    
    # both income tax and business tax, no gross receipts tax
    # gross receipts tax is 0 (https://www.tn.gov/revenue/taxes/gross-receipts-taxes.html)
    # business tax (Class 1E: https://www.tn.gov/revenue/taxes/business-tax/due-dates-and-tax-rates.html)
    if state == 'TN': rate = (net_income*0.065 + sales*0.0003125)/net_income
    
    # gross receipts tax: use retail or wholesaler (https://comptroller.texas.gov/taxes/franchise/)
    if state == 'TX':
        rate = sales*0.00375/net_income
    
    # retroactively reduced from 0.0495 to 0.0485
    if state == 'UT': rate = 0.0485
    
    if state == 'VA': rate = 0.06
    
    if state == 'VT':
        if net_income <= 10000: rate = 0.06
        elif net_income <= 25000: rate = (10000*0.06 + (net_income - 10000)*0.07)/net_income
        else: rate = (10000*0.06 + (25000 - 10000)*0.07 + (net_income - 25000)*0.085)/net_income
    
    # gross receipts tax: use value for Manufacturing of Wood Biomass Fuel (https://dor.wa.gov/taxes-rates/business-occupation-tax/business-occupation-tax-classifications)
    if state == 'WA':
        rate = sales*0.00138/net_income
    
    if state == 'WI': rate = 0.079
    
    if state == 'WV': rate = 0.065
    
    # no income tax and no gross receipts tax
    if state == 'WY': rate = 0
    
    return rate

# https://taxfoundation.org/data/all/state/state-corporate-income-tax-rates-brackets-2023/
def state_income_tax_rate_2023(state='IL', sales=20000, net_income=10000):
    if state == 'AL': rate = 0.065
    
    if state == 'AK':
        if net_income <= 25000: rate = 0
        elif net_income <= 49000: rate = (25000*0 + (net_income - 25000)*0.02)/net_income
        elif net_income <= 74000: rate = (25000*0 + (49000 - 25000)*0.02 + (net_income - 49000)*0.03)/net_income
        elif net_income <= 99000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (net_income - 74000)*0.04)/net_income
        elif net_income <= 124000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (net_income - 99000)*0.05)/net_income 
        elif net_income <= 148000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (net_income - 124000)*0.06)/net_income
        elif net_income <= 173000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (net_income - 148000)*0.07)/net_income
        elif net_income <= 198000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (173000 - 148000)*0.07 + (net_income - 173000)*0.08)/net_income
        elif net_income <= 222000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (173000 - 148000)*0.07 + (198000 - 173000)*0.08 + (net_income - 198000)*0.09)/net_income
        else: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (173000 - 148000)*0.07 + (198000 - 173000)*0.08 + (222000 - 198000)*0.09 + (net_income - 222000)*0.094)/net_income
    
    if state == 'AR':
        if net_income <= 3000: rate = 0.01
        elif net_income <= 6000: rate = (3000*0.01 + (net_income - 3000)*0.02)/net_income
        elif net_income <= 11000: rate = (3000*0.01 + (6000 - 3000)*0.02 + (net_income - 6000)*0.03)/net_income
        elif net_income <= 25000: rate = (3000*0.01 + (6000 - 3000)*0.02 + (11000 - 6000)*0.03 + (net_income - 11000)*0.05)/net_income
        else: rate = (3000*0.01 + (6000 - 3000)*0.02 + (11000 - 6000)*0.03 + (25000 - 11000)*0.05 + (net_income - 25000)*0.053)/net_income
    
    if state == 'AZ': rate = 0.049
    
    if state == 'CA': rate = 0.0884
    
    if state == 'CO': rate = 0.044
    
    if state == 'CT': rate = 0.075
    
    if state == 'DC': rate = 0.0825
    
    # both income tax and gross receipts tax
    # gross receipts tax: use Petroleum wholesaler (https://revenuefiles.delaware.gov/TaxTips/tt-petroleum_dealers2018.pdf)
    if state == 'DE': rate = (net_income*0.087 + sales*0.003983)/net_income
    
    if state == 'FL': rate = 0.055
    
    if state == 'GA': rate = 0.0575
    
    if state == 'HI':
        if net_income <= 25000: rate = 0.044
        elif net_income <= 100000: rate = (25000*0.044 + (net_income - 25000)*0.054)/net_income
        else: rate = (25000*0.044 + (100000 - 25000)*0.054 + (net_income - 100000)*0.064)/net_income
    
    if state == 'IA':
        if net_income <= 100000: rate = 0.055
        else: rate = (100000*0.055 + (net_income - 100000)*0.084)/net_income
    
    if state == 'ID': rate = 0.058
    
    if state == 'IL': rate = 0.095
    
    if state == 'IN': rate = 0.049
    
    if state == 'KS':
        if net_income <= 50000: rate = 0.04
        else: rate = (50000*0.04 + (net_income - 50000)*0.07)/net_income
    
    if state == 'KY': rate = 0.05
    
    if state == 'LA':
        if net_income <= 50000: rate = 0.035
        elif net_income <= 150000: rate = (50000*0.035 + (net_income - 50000)*0.055)/net_income
        else: rate = (50000*0.035 + (150000 - 50000)*0.055 + (net_income - 150000)*0.075)/net_income
    
    if state == 'MA': rate = 0.08
    
    if state == 'MD': rate = 0.0825
    
    if state == 'ME':
        if net_income <= 350000: rate = 0.035
        elif net_income <= 1050000: rate = (350000*0.035 + (net_income - 350000)*0.0793)/net_income
        elif net_income <= 3500000: rate = (350000*0.035 + (1050000 - 350000)*0.0793 + (net_income - 1050000)*0.0833)/net_income
        else: rate = (350000*0.035 + (1050000 - 350000)*0.0793 + (3500000 - 1050000)*0.0833 + (net_income - 3500000)*0.0893)/net_income
    
    if state == 'MI': rate = 0.06
    
    if state == 'MN': rate = 0.098
    
    if state == 'MO': rate = 0.04
    
    if state == 'MS':
        if net_income <= 5000: rate = 0
        elif net_income <= 10000: rate = (5000*0 + (net_income - 5000)*0.04)/net_income
        else: rate = (5000*0 + (10000 - 5000)*0.04 + (net_income - 10000)*0.05)/net_income
    
    if state == 'MT': rate = 0.0675
    
    if state == 'NC': rate = 0.025
    
    if state == 'ND':
        if net_income <= 25000: rate = 0.0141
        elif net_income <= 50000: rate = (25000*0.0141 + (net_income - 25000)*0.0355)/net_income
        else: rate = (25000*0.0141 + (50000 - 25000)*0.0355 + (net_income - 50000)*0.0431)/net_income
    
    if state == 'NE':
        if net_income <= 100000: rate = 0.0558
        else: rate = (100000*0.0558 + (net_income - 100000)*0.0725)/net_income
    
    if state == 'NH': rate = 0.075
    
    # note in NJ, the rates apply to the entire net income rather than just income over the threshold
    if state == 'NJ':
        if net_income <= 50000: rate = 0.065
        elif net_income <= 100000: rate = 0.075
        elif net_income <= 1000000: rate = 0.09
        else: rate = 0.115
    
    if state == 'NM':
        if net_income <= 500000: rate = 0.048
        else: rate = (500000*0.048 + (net_income - 500000)*0.059)/net_income
    
    # gross receipts tax: use value for waste management and remediation services (https://tax.nv.gov/uploadedFiles/taxnvgov/Content/Commerce/Commerce_Tax_Instructions.pdf, link invalid)
    if state == 'NV':
        rate = 0 if sales <= 4000000 else (sales - 4000000)*0.00261/net_income
    
    if state == 'NY':
        if net_income <= 5000000: rate = 0.065
        else: rate = (5000000*0.065 + (net_income - 5000000)*0.0725)/net_income
    
    # gross receipts tax (https://taxfoundation.org/data/all/state/state-gross-receipts-taxes-2023/)
    if state == 'OH':
        rate = sales*0.0026/net_income
    
    if state == 'OK': rate = 0.04
    
    # both income tax and gross receipts tax (https://taxfoundation.org/data/all/state/state-gross-receipts-taxes-2023/)
    if state == 'OR':
        if net_income <= 1000000: rate = 0.066
        else: rate = (1000000*0.066 + (net_income - 1000000)*0.076)/net_income
        rate = (net_income*rate + sales*0.0057)/net_income
    
    if state == 'PA': rate = 0.0899
    
    if state == 'RI': rate = 0.07
    
    if state == 'SC': rate = 0.05
    
    # no income tax and no gross receipts tax
    if state == 'SD': rate = 0
    
    # both income tax and business tax, no gross receipts tax
    # gross receipts tax is 0 (https://www.tn.gov/revenue/taxes/gross-receipts-taxes.html)
    # business tax (Class 1E: https://www.tn.gov/revenue/taxes/business-tax/due-dates-and-tax-rates.html)
    if state == 'TN': rate = (net_income*0.065 + sales*0.0003125)/net_income
    
    # gross receipts tax: use retail or wholesaler (https://comptroller.texas.gov/taxes/franchise/)
    if state == 'TX':
        rate = sales*0.00375/net_income
    
    # retroactively reduced from 0.0485 to 0.0465
    if state == 'UT': rate = 0.0465
    
    if state == 'VA': rate = 0.06
    
    if state == 'VT':
        if net_income <= 10000: rate = 0.06
        elif net_income <= 25000: rate = (10000*0.06 + (net_income - 10000)*0.07)/net_income
        else: rate = (10000*0.06 + (25000 - 10000)*0.07 + (net_income - 25000)*0.085)/net_income
    
    # gross receipts tax: use value for Manufacturing of Wood Biomass Fuel (https://dor.wa.gov/taxes-rates/business-occupation-tax/business-occupation-tax-classifications)
    if state == 'WA':
        rate = sales*0.00138/net_income
    
    if state == 'WI': rate = 0.079
    
    if state == 'WV': rate = 0.065
    
    # no income tax and no gross receipts tax
    if state == 'WY': rate = 0
    
    return rate

# https://taxfoundation.org/data/all/state/state-corporate-income-tax-rates-brackets-2024/
def state_income_tax_rate_2024(state='IL', sales=20000, net_income=10000):
    if state == 'AL': rate = 0.065
    
    if state == 'AK':
        if net_income <= 25000: rate = 0
        elif net_income <= 49000: rate = (25000*0 + (net_income - 25000)*0.02)/net_income
        elif net_income <= 74000: rate = (25000*0 + (49000 - 25000)*0.02 + (net_income - 49000)*0.03)/net_income
        elif net_income <= 99000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (net_income - 74000)*0.04)/net_income
        elif net_income <= 124000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (net_income - 99000)*0.05)/net_income 
        elif net_income <= 148000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (net_income - 124000)*0.06)/net_income
        elif net_income <= 173000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (net_income - 148000)*0.07)/net_income
        elif net_income <= 198000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (173000 - 148000)*0.07 + (net_income - 173000)*0.08)/net_income
        elif net_income <= 222000: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (173000 - 148000)*0.07 + (198000 - 173000)*0.08 + (net_income - 198000)*0.09)/net_income
        else: rate = (25000*0 + (49000 - 25000)*0.02 + (74000 - 49000)*0.03 + (99000 - 74000)*0.04 + (124000 - 99000)*0.05 + (148000 - 124000)*0.06 + (173000 - 148000)*0.07 + (198000 - 173000)*0.08 + (222000 - 198000)*0.09 + (net_income - 222000)*0.094)/net_income
    
    if state == 'AR':
        if net_income <= 3000: rate = 0.01
        elif net_income <= 6000: rate = (3000*0.01 + (net_income - 3000)*0.02)/net_income
        elif net_income <= 11000: rate = (3000*0.01 + (6000 - 3000)*0.02 + (net_income - 6000)*0.03)/net_income
        else: rate = (3000*0.01 + (6000 - 3000)*0.02 + (11000 - 6000)*0.03 + (net_income - 11000)*0.048)/net_income
    
    if state == 'AZ': rate = 0.049
    
    if state == 'CA': rate = 0.0884
    
    if state == 'CO': rate = 0.044
    
    if state == 'CT': rate = 0.075
    
    if state == 'DC': rate = 0.0825
    
    # both income tax and gross receipts tax
    # gross receipts tax: use Petroleum wholesaler (https://revenuefiles.delaware.gov/TaxTips/tt-petroleum_dealers2018.pdf)
    if state == 'DE': rate = (net_income*0.087 + sales*0.003983)/net_income
    
    if state == 'FL': rate = 0.055
    
    if state == 'GA': rate = 0.0575
    
    if state == 'HI':
        if net_income <= 25000: rate = 0.044
        elif net_income <= 100000: rate = (25000*0.044 + (net_income - 25000)*0.054)/net_income
        else: rate = (25000*0.044 + (100000 - 25000)*0.054 + (net_income - 100000)*0.064)/net_income
    
    if state == 'IA':
        if net_income <= 100000: rate = 0.055
        else: rate = (100000*0.055 + (net_income - 100000)*0.071)/net_income
    
    if state == 'ID': rate = 0.058
    
    if state == 'IL': rate = 0.095
    
    if state == 'IN': rate = 0.049
    
    if state == 'KS':
        if net_income <= 50000: rate = 0.035
        else: rate = (50000*0.035 + (net_income - 50000)*0.065)/net_income
    
    if state == 'KY': rate = 0.05
    
    if state == 'LA':
        if net_income <= 50000: rate = 0.035
        elif net_income <= 150000: rate = (50000*0.035 + (net_income - 50000)*0.055)/net_income
        else: rate = (50000*0.035 + (150000 - 50000)*0.055 + (net_income - 150000)*0.075)/net_income
    
    if state == 'MA': rate = 0.08
    
    if state == 'MD': rate = 0.0825
    
    if state == 'ME':
        if net_income <= 350000: rate = 0.035
        elif net_income <= 1050000: rate = (350000*0.035 + (net_income - 350000)*0.0793)/net_income
        elif net_income <= 3500000: rate = (350000*0.035 + (1050000 - 350000)*0.0793 + (net_income - 1050000)*0.0833)/net_income
        else: rate = (350000*0.035 + (1050000 - 350000)*0.0793 + (3500000 - 1050000)*0.0833 + (net_income - 3500000)*0.0893)/net_income
    
    if state == 'MI': rate = 0.06
    
    if state == 'MN': rate = 0.098
    
    if state == 'MO': rate = 0.04
    
    if state == 'MS':
        if net_income <= 5000: rate = 0
        elif net_income <= 10000: rate = (5000*0 + (net_income - 5000)*0.04)/net_income
        else: rate = (5000*0 + (10000 - 5000)*0.04 + (net_income - 10000)*0.05)/net_income
    
    if state == 'MT': rate = 0.0675
    
    if state == 'NC': rate = 0.025
    
    if state == 'ND':
        if net_income <= 25000: rate = 0.0141
        elif net_income <= 50000: rate = (25000*0.0141 + (net_income - 25000)*0.0355)/net_income
        else: rate = (25000*0.0141 + (50000 - 25000)*0.0355 + (net_income - 50000)*0.0431)/net_income
    
    if state == 'NE':
        if net_income <= 100000: rate = 0.0558
        else: rate = (100000*0.0558 + (net_income - 100000)*0.0584)/net_income
    
    if state == 'NH': rate = 0.075
    
    # note in NJ, the rates apply to the entire net income rather than just income over the threshold
    if state == 'NJ':
        if net_income <= 50000: rate = 0.065
        elif net_income <= 100000: rate = 0.075
        else: rate = 0.09
    
    if state == 'NM':
        if net_income <= 500000: rate = 0.048
        else: rate = (500000*0.048 + (net_income - 500000)*0.059)/net_income
    
    # gross receipts tax: use value for waste management and remediation services (https://tax.nv.gov/uploadedFiles/taxnvgov/Content/Commerce/Commerce_Tax_Instructions.pdf, link invalid)
    if state == 'NV':
        rate = 0 if sales <= 4000000 else (sales - 4000000)*0.00261/net_income
    
    if state == 'NY':
        if net_income <= 5000000: rate = 0.065
        else: rate = (5000000*0.065 + (net_income - 5000000)*0.0725)/net_income
    
    # gross receipts tax (https://taxfoundation.org/data/all/state/state-gross-receipts-taxes-2024/)
    if state == 'OH':
        rate = sales*0.0026/net_income
    
    if state == 'OK': rate = 0.04
    
    # both income tax and gross receipts tax (https://taxfoundation.org/data/all/state/state-gross-receipts-taxes-2024/)
    if state == 'OR':
        if net_income <= 1000000: rate = 0.066
        else: rate = (1000000*0.066 + (net_income - 1000000)*0.076)/net_income
        rate = (net_income*rate + sales*0.0057)/net_income
    
    if state == 'PA': rate = 0.0849
    
    if state == 'RI': rate = 0.07
    
    if state == 'SC': rate = 0.05
    
    # no income tax and no gross receipts tax
    if state == 'SD': rate = 0
    
    # both income tax and business tax, no gross receipts tax
    # gross receipts tax is 0 (https://www.tn.gov/revenue/taxes/gross-receipts-taxes.html)
    # business tax (Class 1E: https://www.tn.gov/revenue/taxes/business-tax/due-dates-and-tax-rates.html)
    if state == 'TN': rate = (net_income*0.065 + sales*0.0003125)/net_income
    
    # gross receipts tax: use retail or wholesaler (https://comptroller.texas.gov/taxes/franchise/)
    if state == 'TX':
        rate = sales*0.00375/net_income
    
    if state == 'UT': rate = 0.0465
    
    if state == 'VA': rate = 0.06
    
    if state == 'VT':
        if net_income <= 10000: rate = 0.06
        elif net_income <= 25000: rate = (10000*0.06 + (net_income - 10000)*0.07)/net_income
        else: rate = (10000*0.06 + (25000 - 10000)*0.07 + (net_income - 25000)*0.085)/net_income
    
    # gross receipts tax: use value for Manufacturing of Wood Biomass Fuel (https://dor.wa.gov/taxes-rates/business-occupation-tax/business-occupation-tax-classifications)
    if state == 'WA':
        rate = sales*0.00138/net_income
    
    if state == 'WI': rate = 0.079
    
    if state == 'WV': rate = 0.065
    
    # no income tax and no gross receipts tax
    if state == 'WY': rate = 0
    
    return rate