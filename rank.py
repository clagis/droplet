#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 18:06:40 2022

@author: clement
"""

# =============================================================================
## Library dependancies
# =============================================================================

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
from clusteringCoefficient import unweightedCC, weightedCC




def zigzag(dir_in, lst, unweighted = True, SJ = False):
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
    
    M = []
    for name in lst:
        df_lst = dfListConnex(dir_in, name)
        for k in range(len(df_lst)):
            rankList = [0 for i in range(6)]
            for j in range(len(df_lst[k])):
                L = sorted(df_lst[k].iloc[j]['members'])
                count = 1
                if len(L) == 1:
                    rankList[count] += 1
                else:
                    while len(L) > 1:
                        testPolymer = L[0][0]
                        testBead = L[0][1]
                        for elem in L[1:]:
                            if (testPolymer == elem[0]) and (elem[1] > (testBead+1)):
                                count += 1
                                testBead = elem[1]
                                idx = L.index(elem)
                                L.pop(idx)
                        L.pop(0)
                        rankList[count] += 1
                        count = 1
            rankList2 = [x/sum(rankList)*100 for x in rankList]
            M += [[name, 50000*(k+1), rankList2[1:]]]
    
    df = pd.DataFrame(M, columns=['run', 'time', 'rankList'])
    return df
                            
                            
R = ['2B/6091',
     '4C/4025', 
     '5D/11003',
     '6E/14003']

df2 = zigzag('neo_junctions',R)

def meanList(L):
    s = 0
    for i in range(len(L)):
        s += (i+1)*L[i]/100
    return s

def stdList(L):
    s = 0
    for i in range(len(L)):
        s += (i+1)*L[i]/100
    e = 0
    for j in range(len(L)):
        e += ((i+1 - s)**2)*L[i]/100
    return np.sqrt(e)

def expandList(L,k):
    return L[k]

df2['rank1'] = df2['rankList'].apply(lambda x : expandList(x,0))
df2['rank2'] = df2['rankList'].apply(lambda x : expandList(x,1))
df2['rank3'] = df2['rankList'].apply(lambda x : expandList(x,2))
df2['rank4'] = df2['rankList'].apply(lambda x : expandList(x,3))
df2['rank5'] = df2['rankList'].apply(lambda x : expandList(x,4))