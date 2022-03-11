#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 14:15:17 2021

@author: clement
"""

# =============================================================================
## Library dependancies
# =============================================================================

import numpy as np
import networkx as nx

# =============================================================================
## Basic functions
# =============================================================================


def multigraphClustering(G):
    """
    Return the list of multigraph clustering coefficient of each node of the
    graph G. The multigraph clustering coefficient is defined as in the same 
    manner as the clustering coefficient but with each edge of the triangles 
    multiplied by the number of actual links between the two extremal junctions.

    Parameters:
    -----------
    G: networkx multigraph

    Returns:
    --------
    C: list
        A list containing the multigraph clustering coefficient of each node of 
        the graph G.
    """
    
    A = nx.adjacency_matrix(G)
    A = A.todense()
    n = A.shape[0]
    C = []
    
    p = 0
    for i in range(n):
        for j in range(n):
            p += A[i,j]
    
    for i in range(n):
        A[i,i] = 0
    
    T = np.matmul(np.matmul(A,A),A)
    for i in range(n):
        C_i = T[i,i]
        k_i = sum([min(A[i,j],1) for j in range(n)])
          
        C_i /= max((k_i-1)*(k_i),1)
        C += [C_i]
    
    return C


def unweightedCC(data_edges):
    """
    Return a graph and its average clustering coefficient.

    Parameters:
    -----------
    data_edges: pandas dataframe
        data must contain at least three columns named source, target and links.   

    Returns:
    --------
    G: networkx graph
        The graph corresponfing to the input data.
    avgUWCC: float
        Average clustering coefficient of G.
    """
    
    G = nx.from_pandas_edgelist(data_edges, 'source', 'target', create_using=nx.Graph())
    avgUWCC = nx.average_clustering(G)
    
    return G, avgUWCC


def weightedCC(data_edges):
    """
    Return a multigraph and its average multigraph clustering coefficient.

    Parameters:
    -----------
    data_edges: pandas dataframe
        data must contain at least three columns named source, target and links.   

    Returns:
    --------
    G: networkx multigraph
        The graph corresponfing to the input data.
    avgWCC: float
        Average multigrah clustering coefficient of G.
    """
    
    G = nx.from_pandas_edgelist(data_edges, 'source', 'target', create_using=nx.MultiGraph())
    M = multigraphClustering(G)
    avgWCC = sum(M)/max(len(M),1)

    return G, avgWCC