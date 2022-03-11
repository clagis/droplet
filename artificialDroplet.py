#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 14:50:39 2021

@author: clement
"""

# =============================================================================
## Library dependancies
# =============================================================================

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib as mpl

from utilitaries import simpleToEdgeList
from clusteringCoefficient import weightedCC, unweightedCC, multigraphClustering

# =============================================================================
## Basic functions
# =============================================================================

def euclideanDistance(p0,p1):
    """
    Return the Euclidean distance between two points.

    Parameters:
    -----------
    p0: pandas dataframe
        p0 must coantain three columns labeled 'x', 'y' and 'z'.
    p1: pandas dataframe
        p1 must coantain three columns labeled 'x', 'y' and 'z'.

    Returns:
    --------
    np.sqrt(x**2 + y**2 + z**2): float
        The Euclidean distance between p0 and p1.
    """
    
    x, y, z = (p0-p1)[['x','y','z']]
    
    return np.sqrt(x**2 + y**2 + z**2)


def artificialDroplet(size, max_dist):
    """
    Return a pandas dataframe of lattice points linked to all lattice points
    in a sphere of radius max_dist.

    Parameters:
    -----------
    size: integer
        The approximative number of junctions in the dataframe. In facts, it 
        will use the nearest cube.
    max_dist: float
        The maximal distance between two neighbors. If max_dist < 1, then the 
        droplet will be totaly disconnected. Typical max_dist are 1.5 (> sqrt(2))
        and 1.8 (> sqrt(3)).

    Returns:
    --------
    df: pandas dataframe
        A dataframe containing the position of the junctions and their neighbors
        in the droplet.
    """
    
    Coordinates = []
    n = round(np.cbrt(size))
    
    count = 0
    for i in range(n):
        for j in range(n):
            for k in range(n):
                Coordinates += [[count, i, j, k]]
                count += 1
    
    df = pd.DataFrame(Coordinates, columns=['junction', 'x', 'y', 'z'])
    
    Neighbors = []
    
    for i in range(n**3):
        Neighbors_i = []
        for j in range(n**3):
            if 0 < euclideanDistance(df.iloc[i],df.iloc[j]) <= max_dist:
                Neighbors_i += [j]
        Neighbors += [Neighbors_i]
    
    df['neighbors'] = Neighbors
    return df


# df = artificialDroplet(125, 1.5)
# df_edges = simpleToEdgeList(df, 2)
# G, avgWCC = weightedCC(df_edges)
# color_lookup = {k:v for k, v in enumerate([value for value in nx.clustering(G).values()])}

# low, *_, high = sorted(color_lookup.values())
# norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
# mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)
# nx.draw(G,
#         node_color=[mapper.to_rgba(i) 
#                     for i in color_lookup.values()]color
#         )


df = artificialDroplet(125, 1.5)

dct = {}
# for k in range(1,2):
df_edges = simpleToEdgeList(df, 1)
G, avgWCC = unweightedCC(df_edges)
dct[1] = avgWCC
color_lookup = {k:v for k, v in enumerate([value for value in multigraphClustering(G)])}

low, *_, high = sorted(color_lookup.values())
norm = mpl.colors.Normalize(vmin=low, vmax=high, clip=True)
mapper = mpl.cm.ScalarMappable(norm=norm, cmap=mpl.cm.coolwarm)
nx.draw(G,
        node_color=[mapper.to_rgba(i) 
                    for i in color_lookup.values()]
        )

def squareCC(size):
    """
    Return the average clustering coefficient of a cubic droplet of a certain 
    size, when the square diagonals link the junctions.

    Parameters:
    -----------
    size: integer
        The size of the side of the cube.

    Returns:
    --------
    CC/size**3: float
        The average clustering coefficient of the dataframe.
    """
    vertice = 4/5
    edge = 7/12
    facet = 6/13
    core = 20/51
    
    CC = 8*vertice + 12*(size-2)*edge + 6*(size-2)**2*facet + (size-2)**3*core
    return CC/size**3
    
def cubeCC(size):
    """
    Return the average clustering coefficient of a cubic droplet of a certain 
    size, when the cube diagonals link the junctions.

    Parameters:
    -----------
    size: integer
        The size of the side of the cube.

    Returns:
    --------
    CC/size**3: float
        The average clustering coefficient of the dataframe.
    """
    vertice = 1
    edge = 39/55
    facet = 9/17
    core = 132/325
    
    CC = 8*vertice + 12*(size-2)*edge + 6*(size-2)**2*facet + (size-2)**3*core
    return CC/size**3