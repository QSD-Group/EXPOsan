# -*- coding: utf-8 -*-
"""
EXPOsan: Exposition of sanitation and resource recovery systems

This module is developed by:
    
    Joy Zhang <joycheung1994@gmail.com>

This module is under the University of Illinois/NCSA Open Source License.
Please refer to https://github.com/QSD-Group/EXPOsan/blob/main/LICENSE.txt
for license details.

"""

from math import sqrt
from itertools import permutations
import numpy as np, scipy.spatial.distance as spa
# from qsdsan.utils import time_printer

#%%
class Location:
    def __init__(self, name, x, y):
        self.name = name
        self.x = x
        self.y = y
    
    def distance(self, other):
        return sqrt((self.x - other.x)**2 + (self.y - other.y)**2)

    def __repr__(self):
        return f'<Location {self.name}: x={self.x:.2f}, y={self.y:.2f}>\n'

def _distance_matrix(xs, ys, metric='euclidean'):
    coords = np.vstack((xs, ys)).T
    arr = spa.pdist(coords)
    return spa.squareform(arr)

class Locations:
    @classmethod
    def random_within_area(cls, n=5, x_range=(0,1), y_range=(0,1), seed=None, 
                           distance_metric='euclidean'):
        rng = np.random.default_rng(seed)
        xs = rng.uniform(*x_range, size=n)
        ys = rng.uniform(*y_range, size=n)
        return cls.from_coordinates(xs, ys, distance_metric)
    
    @classmethod
    def from_coordinates(cls, xs, ys, distance_metric='euclidean'):
        new = cls(distance_metric=distance_metric)
        new.xs = np.asarray(xs)
        new.ys = np.asarray(ys)
        names = []
        for i in range(len(xs)):
            ID = f'l{i}'
            names.append(ID)
            setattr(new, ID, Location(ID, xs[i], ys[i]))
        new.size = len(xs)
        new.names = names
        new.index = dict(zip(range(new.size), new.names))
        return new
    
    def __init__(self, locations=(), distance_metric='euclidean'):
        xs = []
        ys = []
        self._M_dist = None
        self._metric=distance_metric
        for loc in locations:
            setattr(self, loc.name, loc)
            xs.append(loc.x)
            ys.append(loc.y)
        self.xs = np.array(xs)
        self.ys = np.array(ys)
        self.size = len(locations)
        self.names = [loc.name for loc in locations]
        self.index = dict(zip(range(self.size), self.names))
    
    @property
    def distance_metric(self):
        return self._metric
    @distance_metric.setter
    def distance_metric(self, m):
        if m != self._metric:
            self._metric = m
            self._M_dist = None
    @property
    def distance_matrix(self):
        if self._M_dist is None:
            self._M_dist = _distance_matrix(self.xs, self.ys, self.distance_metric)
        return self._M_dist
    
    def __repr__(self):
        return f'<Locations: {self.names}>\n'

#%%
def solve_tsp(locations):
    M_dist = locations.distance_matrix()

    # Solve the TSP using a brute force algorithm
    best_tour = None
    best_distance = float('inf')
    tours = permutations(range(1, locations.size))
    for tour in tours:
        row = [0, *tour]
        col = [*tour, 0]
        distance = sum(M_dist[i][j] for i,j in zip(row, col))
        # print(f'tour: {row}, distance: {distance}')
        if distance < best_distance:
            best_tour = row
            best_distance = distance
            
    return best_tour, best_distance


#%%
n = 100
locs = Locations.random_within_area()
# bt, bd = solve_tsp(locs)
