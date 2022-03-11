#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 10:14:10 2022

@author: clement
"""

import pandas as pd
import numpy as np
import seaborn as sns
import networkx as nx
import os
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3

from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch
from pathlib import Path
from utilitaries import toEdgeList, dfListConnex, includeSJ
from clusteringCoefficient import unweightedCC, weightedCC, multigraphClustering




def graphLinks(dir_in, lst, unweighted = True, SJ = False):
    """
    Return a seaborn plot showing the mass of the junction with or without the
    loop counted and with or without the small junctions.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the directory named name.
    lst: list of strings
        List of subdirectories in the input directory in which the files lie.
    unweighted: boolean, default True
        Indicate if you want the weighted or unweighted degree.
    SJ: boolean, default False
        Indicate if you want the small junctions counted or not.

    Returns:
    --------
    figure: seaborn violinplot
        The plot showing the distribution of the average clustering coefficient
        of the network.
    """
    
    
    
    affinityList = []
    typeList = []


    for name in lst:
        df_lst = dfListConnex(dir_in, name)
        n = len(df_lst)
        m = int(np.ceil(n/2))
        fig = plt.figure()
        axes = fig.subplots(nrows=2, ncols=m)
        typeList += [str(df_lst[0].at[0,'type'])]
        for k in range(len(df_lst)):
            df_edges = toEdgeList(includeSJ(df_lst[k], SJ))
            G, avgUWCC = weightedCC(df_edges)
            
            N = list(G.nodes())
            E = list(G.edges())
            p = len(N)
            
            
            links = []
            for i in range(p):
                for j in range(i+1,p):
                    links += [E.count((i,j))]
            
            lightLinks = [elem for elem in links if elem != 0]
            q = len(lightLinks)
            L = []
            for elem in set(lightLinks):
                L += [[elem, lightLinks.count(elem)/q*100]]
    
            df = pd.DataFrame(L,columns=['links', 'percentage'])
            if k < m:
                a = 0
                b = k
            else:
                a = 1
                b = k-m
            sns.barplot(x='links', y = 'percentage', data=df, ax = axes[a, b]).set(title='{}'.format((k+1)*50000))
        
        fig.suptitle(name, fontsize=20)
        plt.show()
    
    return E





R = ['2B/6091',
     '4C/4025', 
     '5D/11003',
     '6E/14003']

def artificialLinks(dir_in, lst, unweighted = True, SJ = False):
    L = []
    for name in lst:
        df_lst = dfListConnex(dir_in, name)
        for k in range(len(df_lst)):
            data_edges = toEdgeList(includeSJ(df_lst[k], SJ))
            if unweighted:
                G = nx.from_pandas_edgelist(data_edges, 'source', 'target', create_using=nx.Graph())
            else:
                G = nx.from_pandas_edgelist(data_edges, 'source', 'target', create_using=nx.MultiGraph())
            e = len(G.edges())
            n = len(G.nodes())
            q = e//n
            r = e - n*q
            
            H = nx.MultiGraph()
            H.add_nodes_from([i for i in range(n)])
            for i in range(n):
                for j in range(1,q+1):
                    if not ((i,(i+j) % n) in H.edges()):
                        H.add_edge(i,(i+j) % n)
            
            for j in range(r):
                H.add_edge(j,(j+q+1) % n)

            M = multigraphClustering(H)
            avgWCC = sum(M)/max(len(M),1)
            L += [[name, avgWCC]]
    
    if unweighted:
        df = pd.DataFrame(L,columns=['run','avgCC'])
        sns.violinplot(x='run',y='avgCC', data=df)
    else:
        df = pd.DataFrame(L,columns=['run','avgWCC'])
        sns.violinplot(x='run',y='avgWCC', data=df)
        

for k in range(4):
    graphLinks('neo_junctions',R[k:k+1])
    