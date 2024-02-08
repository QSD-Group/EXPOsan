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
from ortools.constraint_solver import routing_enums_pb2
from ortools.constraint_solver import pywrapcp
from math import ceil

__all__ = ('Location', 'Locations', 'SimpleCVRP')

#%%
class Location:
    def __init__(self, name, x, y, demand=None):
        self.name = name
        self.x = x
        self.y = y
        self.demand = demand
    
    def distance(self, other):
        return sqrt((self.x - other.x)**2 + (self.y - other.y)**2)

    def __repr__(self):
        return f'<Location {self.name}: x={self.x:.2f}, y={self.y:.2f}>\n'


class Locations:
    
    def __init__(self, locations=(), distance_metric='euclidean'):
        xs = []
        ys = []
        ds = []
        self._M_dist = None
        self._metric=distance_metric
        for loc in locations:
            setattr(self, loc.name, loc)
            xs.append(loc.x)
            ys.append(loc.y)
            d = loc.demand if loc.demand else 0
            ds.append(d)
        self.xs = np.array(xs)
        self.ys = np.array(ys)
        self.demands = np.array(ds)
        self.size = len(locations)
        self.names = [loc.name for loc in locations]
        self.index = dict(zip(range(self.size), self.names))

    @classmethod
    def from_coordinates(cls, xs, ys, demands=None, distance_metric='euclidean'):
        new = cls(distance_metric=distance_metric)
        new.xs = np.asarray(xs)
        new.ys = np.asarray(ys)
        if demands is not None:
            assert len(demands) == len(xs)
            new.demands = np.asarray(demands)
        names = []
        for i in range(len(xs)):
            ID = f'l{i}'
            names.append(ID)
            setattr(new, ID, Location(ID, xs[i], ys[i], demands[i]))
        new.size = len(xs)
        new.names = names
        new.index = dict(zip(range(new.size), new.names))
        return new
    
    @classmethod
    def random_within_area(cls, n=5, x_range=(0,1), y_range=(0,1), demands=None,
                           seed=None, distance_metric='euclidean'):
        rng = np.random.default_rng(seed)
        xs = rng.uniform(*x_range, size=n)
        ys = rng.uniform(*y_range, size=n)
        if demands is not None:
            if isinstance(demands, (int, float)):
                demands = [0] + [demands] * (n-1)
        return cls.from_coordinates(xs, ys, demands, distance_metric)

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
        return sum(M[path[i]][path[i+1]] for i in range(len(path)-1))
    
    def plot_path(self, path):
        fig, ax = plt.subplots()
        xs, ys = self.xs, self.ys
        for i in range(len(path)-1):
            locs = path[i: i+2]
            ax.plot(xs[locs], ys[locs], c='r', 
                    marker='o', mfc='b', mec='b')
        ax.plot(xs[path[0]], ys[path[0]], marker='*', markersize=12,
                mfc='yellow', mec='black')
        ax.set_aspect('equal')
        return fig, ax
    
#%%
class SimpleCVRP:
    
    def __init__(self, locations, vehicle_capcity):
        self._locations = locations
        self._dm = dm = np.int64(locations.demands)
        self.vehicle_capcity = vc = int(vehicle_capcity)
        self._n_vehicle = nv = ceil(sum(dm) / vc)
        self._M_dist = np.int64(locations.distance_matrix)
        self._capacities = [vehicle_capcity] * nv
        self.manager = manager = pywrapcp.RoutingIndexManager(locations.size, nv, 0)
        self.routing = pywrapcp.RoutingModel(manager)
        self.register()
    
    def register(self):
        M_dist = self._M_dist
        manager = self.manager
        routing = self.routing
        dm = self._dm
        capacities = self._capacities

        def distance_callback(from_index, to_index):
            from_node = manager.IndexToNode(from_index)
            to_node = manager.IndexToNode(to_index)
            return M_dist[from_node][to_node]

        transit_callback_index = routing.RegisterTransitCallback(distance_callback)
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

        def demand_callback(from_index):
            from_node = manager.IndexToNode(from_index)
            return dm[from_node]

        demand_callback_index = routing.RegisterUnaryTransitCallback(demand_callback)
        routing.AddDimensionWithVehicleCapacity(
            demand_callback_index,
            0,              # null capacity slack
            capacities,     # vehicle maximum capacities
            True,           # start cumul to zero
            "Capacity",
        )
        
        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.first_solution_strategy = (
            routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
        )
        # search_parameters.first_solution_strategy = (
        #     routing_enums_pb2.FirstSolutionStrategy.CHRISTOFIDES
        # )
        search_parameters.local_search_metaheuristic = (
            routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
        )
        self.search_parameters = search_parameters 
        self.solution = None
        self._routes = None
    
    def solve(self, time_limit=300, print_solution=True): 
        self._routes = None
        params = self.search_parameters
        routing = self.routing
        params.time_limit.FromSeconds(time_limit)
        self.solution = routing.SolveWithParameters(params)
        if print_solution: self._print_solution()
        
    @property
    def routes(self):
        if not self._routes: 
            self._get_routes()
        return self._routes
    
    def _get_routes(self):
        manager = self.manager
        routing = self.routing
        solution = self.solution
        if solution is None: return
        routes = []
        for vehicle in range(self._n_vehicle):
            index = routing.Start(vehicle)
            route = [manager.IndexToNode(index),]
            while not routing.IsEnd(index):
                index = solution.Value(routing.NextVar(index))
                route.append(manager.IndexToNode(index))
            routes.append(route)
        self._routes = routes

    def _print_solution(self):
        """Prints solution on console."""
        routing = self.routing
        locs = self._locations
        print(f"Search status: {routing.status()}\n")
        total_distance = 0
        total_load = 0
        routes = self.routes
        for i, route in enumerate(routes):
            plan_output = f"Route for vehicle {i}:\n"
            plan_output += str(route).replace(',', ' ->')[1:-1]
            distance = locs.path_cost(route)
            load = sum(self._dm[route])
            plan_output += f"\nDistance of the route: {distance:.0f} m\n"
            plan_output += f"Load of the route: {load}\n"
            print(plan_output)
            total_distance += distance
            total_load += load
        print(f"Total distance of all routes: {total_distance:.0f} m")
        print(f"Total load of all routes: {total_load}")
    
    # def plot_routes(self):

#%%
n = 100
seed = 52
locs = Locations.random_within_area(n, seed=seed, distance_metric='cityblock',
                                    x_range=(0, 13.36e3), # in meter
                                    y_range=(0, 11.31e3),
                                    demands=55,)
# path = locs.solve_best_roundtrip()
# dist = locs.path_cost(path)
# locs.plot_path(path)
cvr = SimpleCVRP(locs, vehicle_capcity=2000)
cvr.solve(10)
