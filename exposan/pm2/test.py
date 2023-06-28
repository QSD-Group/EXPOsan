# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 06:53:33 2023

@author: Ga-Yeong Kim
"""

import time
import random

start_time = time.time()
print(start_time)

our_list = list(range(10000000))
element = 7000000

def time_track(t, y):
    global start_time

    elapsed_time = time.time() - start_time

    return elapsed_time-600

#%%
def objective_function():
    global start_time
    start_time = time.time()    # reset start_time to be here

    random_choice = random.choice(our_list)
    while random_choice != element:
        random_choice = random.choice(our_list)


    return obj

