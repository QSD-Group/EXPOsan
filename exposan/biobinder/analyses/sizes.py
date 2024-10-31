#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:

    Ali Ahmad <aa3056@scarletmail.rutgers.edu>

    Yalin Li <mailto.yalin.li@gmail.com>
    
This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.
'''

import numpy as np

# __all__ = ('',)

# =============================================================================
# Feedstock & Biocrude Transportation
# =============================================================================

biocrude_radius = 100 * 1.61  # 100 miles to km, PNNL 29882

def biocrude_distances(N_decentralized_HTL, biocrude_radius):
    """
    Generate a list of distances for biocrude transport from decentralized HTL facilities.
    
    Parameters:
    N_decentralized_HTL (int): Number of decentralized HTL facilities.
    biocrude_radius (float): Maximum distance for transportation.
    
    Returns:
    list: Distances for each facility.
    """
    distances = []
    scale = 45  # scale parameter for the exponential distribution

    for _ in range(N_decentralized_HTL):
        r = np.random.exponential(scale)
        r = min(r, biocrude_radius)  # cap distance at biocrude_radius
        distances.append(r)

    return distances

def total_biocrude_distance(N_decentralized_HTL, biocrude_radius):
    """
    Calculate the total biocrude transportation distance.
    
    Parameters:
    N_decentralized_HTL (int): Number of decentralized HTL facilities.
    biocrude_radius (float): Maximum distance for transportation.
    
    Returns:
    float: Total transportation distance.
    """
    distances = biocrude_distances(N_decentralized_HTL, biocrude_radius)
    total_distance = np.sum(distances)  # Sum of individual distances
    return total_distance

# def simulate_biobinder_and_gwp(N_decentralized_HTL):
#     """
#     Simulates the biobinder's price and calculates its Global Warming Potential (GWP) 
#     along with various financial metrics.

#     Parameters:
#     N_decentralized_HTL (int): The number of decentralized HTL units to simulate.

#         """
#     FeedstockScaler = u.Scaler(
#     'FeedstockScaler', ins=scaled_feedstock, outs='feedstock',
#     scaling_factor=N_decentralized_HTL, reverse=True,
# )

#     FeedstockScaler.simulate()
#     sys.simulate()
  
#     biobinder.price = biobinder_price = tea.solve_price(biobinder)
#     print(f"Number of Reactors: {N_decentralized_HTL}, Biobinder Price: {biobinder_price}")
#     c = qs.currency
#     metrics = {}
#     for attr in ('NPV', 'AOC', 'sales', 'net_earnings'):
#         uom = c if attr in ('NPV', 'CAPEX') else (c + '/yr')
#         metrics[attr] = getattr(tea, attr)  # Use getattr to access attributes dynamically

#     # Calculate allocated impacts for GWP
#     all_impacts = lca.get_allocated_impacts(streams=(biobinder,), operation_only=True, annual=True)
#     GWP = all_impacts['GlobalWarming'] / (biobinder.F_mass * lca.system.operating_hours)
    
#     return biobinder_price, GWP, metrics


if __name__ == '__main__':
    N_range = np.arange(100, 2001, 100)  # Range of HTL reactors
    
    N_decentralized_HTL = 1300
    biocrude_transportation_distance = total_biocrude_distance(N_decentralized_HTL, biocrude_radius)
    print("Total biocrude transportation distance:", biocrude_transportation_distance)
    
#     biobinder_prices = []
#     gwps = []
#     npv_list = []
#     aoc_list = []
#     sales_list = []
#     net_earnings_list = []
    
#     for N in N_range:
#         price, gwp, metrics = simulate_biobinder_and_gwp(N)
#         print("Reactor Count and Corresponding Biobinder Prices:")
#     for N, price in zip(N_range, biobinder_prices):
#         print(f"Reactors: {N}, Price: {price}")
        
#         # Store the results
#         biobinder_prices.append(price)
#         gwps.append(gwp)
#         npv_list.append(metrics['NPV'])
#         aoc_list.append(metrics['AOC'])
#         sales_list.append(metrics['sales'])
#         net_earnings_list.append(metrics['net_earnings'])

#     plt.figure(figsize=(10, 5))
#     plt.plot(N_range, biobinder_prices, marker='o', color='b')
#     plt.title('Biobinder Price vs. Number of Decentralized HTL Reactors')
#     plt.xlabel('Number of HTL Reactors')
#     plt.ylabel('Biobinder Price ($/kg)')
#     plt.grid()
#     plt.tight_layout()
#     plt.show()

#     plt.figure(figsize=(10, 5))
#     plt.plot(N_range, gwps, marker='o', color='g')
#     plt.title('GWP vs. Number of Decentralized HTL Reactors')
#     plt.xlabel('Number of HTL Reactors')
#     plt.ylabel('GWP (kg CO2e/kg)')
#     plt.grid()
#     plt.tight_layout()
#     plt.show()

#     bar_width = 0.2  # Width of the bars
#     index = np.arange(len(N_range))  # X locations for the groups

#     plt.figure(figsize=(10, 5))
#     plt.bar(index - bar_width * 1.5, np.array(npv_list) / 1_000_000, bar_width, label='NPV (millions)', color='blue')
#     plt.bar(index - bar_width / 2, np.array(aoc_list) / 1_000_000, bar_width, label='AOC (millions)', color='orange')
#     plt.bar(index + bar_width / 2, np.array(sales_list) / 1_000_000, bar_width, label='Sales (millions)', color='green')
#     plt.bar(index + bar_width * 1.5, np.array(net_earnings_list) / 1_000_000, bar_width, label='Net Earnings (millions)', color='red')

#     plt.title('Metrics vs. Number of Decentralized HTL Reactors')
#     plt.xlabel('Number of HTL Reactors')
#     plt.ylabel('Value (in millions of dollars)')
#     plt.xticks(index, N_range)
#     plt.legend()
#     plt.grid()
#     plt.tight_layout()
#     plt.show()

