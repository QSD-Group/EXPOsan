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
import numpy as np, scipy.spatial.distance as spa, matplotlib.pyplot as plt
from tsp_solver.greedy import solve_tsp  # pypi: tsp-solver2; https://github.com/dmishin/tsp-solver/tree/master
# from itertools import permutations
# from qsdsan.utils import time_printer

__all__ = ('Location', 'Locations',)

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


class Locations:
    
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
    
    @classmethod
    def random_within_area(cls, n=5, x_range=(0,1), y_range=(0,1), seed=None, 
                           distance_metric='euclidean'):
        rng = np.random.default_rng(seed)
        xs = rng.uniform(*x_range, size=n)
        ys = rng.uniform(*y_range, size=n)
        return cls.from_coordinates(xs, ys, distance_metric)

    
    @property
    def distance_metric(self):
        return self._metric
    @distance_metric.setter
    def distance_metric(self, m):
        if m != self._metric:
            self._metric = m
            self._M_dist = None
    
    def _distance_matrix(self):
        coords = np.vstack((self.xs, self.ys)).T
        arr = spa.pdist(coords, self.distance_metric)
        return spa.squareform(arr)
    
    @property
    def distance_matrix(self):
        if self._M_dist is None:
            self._M_dist = self._distance_matrix()
        return self._M_dist
    
    def __repr__(self):
        return f'<Locations: {self.names}>\n'
    
    def solve_best_roundtrip(self):
        return solve_tsp(self.distance_matrix, endpoints=(0,0))
    
    def path_cost(self, path):
        M = self.distance_matrix
        return sum(M[path[i]][path[i+1]] for i in range(self.size))
    
    def plot_path(self, path):
        fig, ax = plt.subplots()
        xs, ys = self.xs, self.ys
        for i in range(len(path)-1):
            locs = path[i: i+2]
            ax.plot(xs[locs], ys[locs], c='r', 
                    marker='o', mfc='b', mec='b')
        ax.set_aspect('equal')
        return fig, ax

#%%
# @time_printer
# def bf_tsp(locations):
#     M_dist = locations.distance_matrix

#     # Solve the TSP using a brute force algorithm
#     best_tour = None
#     best_distance = float('inf')
#     tours = permutations(range(1, locations.size))
#     for tour in tours:
#         row = [0, *tour]
#         col = [*tour, 0]
#         distance = sum([M_dist[i][j] for i,j in zip(row, col)])
#         # print(f'tour: {row}, distance: {distance}')
#         if distance < best_distance:
#             best_tour = [0, *tour, 0]
#             best_distance = distance
            
#     return best_tour, best_distance

#%%
# bt, bd = bf_tsp(locs)
# path = solve_tsp(locs.distance_matrix, endpoints=(0,0))
# dist = sum(locs.distance_matrix[path[i]][path[i+1]] for i in range(len(path)-1))

#%%
n = 100
seed = 52
locs = Locations.random_within_area(n, seed=seed, distance_metric='cityblock',
                                    x_range=(0, 8.30),
                                    y_range=(0, 7.03))
path = locs.solve_best_roundtrip()
dist = locs.path_cost(path)
locs.plot_path(path)


