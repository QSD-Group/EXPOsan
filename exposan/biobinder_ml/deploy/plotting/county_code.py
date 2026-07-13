# -*- coding: utf-8 -*-
'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Ali Ahmad <aliahmad1331@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import pandas as pd

# Load generated summary file
csv_path = r"C:\Work\Rutgers\QSDsan\EXPOsan\exposan\biobinder_ml\deploy\Data\county_lane_miles_summary.csv"
df = pd.read_csv(csv_path)

# Standardized FIPS Mapping for State Prefixes
state_fips_map = {
    'AL': '01', 'AK': '02', 'AZ': '04', 'AR': '05', 'CA': '06', 'CO': '08', 'CT': '09', 'DE': '10',
    'DC': '11', 'FL': '12', 'GA': '13', 'HI': '15', 'ID': '16', 'IL': '17', 'IN': '18', 'IA': '19',
    'KS': '20', 'KY': '21', 'LA': '22', 'ME': '23', 'MD': '24', 'MA': '25', 'MI': '26', 'MN': '27',
    'MS': '28', 'MO': '29', 'MT': '30', 'NE': '31', 'NV': '32', 'NH': '33', 'NJ': '34', 'NM': '35',
    'NY': '36', 'NC': '37', 'ND': '38', 'OH': '39', 'OK': '40', 'OR': '41', 'PA': '42', 'RI': '44',
    'SC': '45', 'SD': '46', 'TN': '47', 'TX': '48', 'UT': '49', 'VT': '50', 'VA': '51', 'WA': '53',
    'WV': '54', 'WI': '55', 'WY': '56', 'PR': '72'
}

# 1. Map the state code string prefix
df['state_prefix'] = df['state_abbr'].map(state_fips_map).fillna('00')

# 2. Convert county fips to string, pad with leading zeros to 3 digits (e.g., 3 -> 003)
df['county_clean'] = df['county_fips'].astype(int).astype(str).str.zfill(3)

# 3. Combine into the unified national 5-digit GEOID string
df['county_fips'] = df['state_prefix'] + df['county_clean']

# Drop the temporary calculation helper columns
df = df.drop(columns=['state_prefix', 'county_clean'])

# Overwrite with the corrected version
df.to_csv(csv_path, index=False)
print("FIPS Codes successfully repaired to 5-digit standard strings!")