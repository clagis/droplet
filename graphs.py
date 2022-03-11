#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 09:04:55 2021

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


# =============================================================================
## Graphical functions
# =============================================================================

###############################################################################

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

###############################################################################


def network3D(dir_in,name, step):
    """
    A retester
    """
    
    df_lst = dfListConnex(dir_in,name)
    k = step

    fig = plt.figure(k)
    ax = p3.Axes3D(fig)
    ax.set_axis_off()
    sc = ax.scatter(df_lst[k]['x'],df_lst[k]['z'],df_lst[k]['y'],
                    s = df_lst[k]['weight']*10,
                    c = df_lst[k]['small_junction'],
                    cmap = 'bwr_r',
                    alpha=0.8)
    
    for cluster in df_lst[k]['junction']:
        # check if the cluster is sufficiently heavy
        # if not df_lst[k].loc[df_lst[k]['cluster'] == cluster, 'small_cluster'].item():
            L = df_lst[k].loc[df_lst[k]['junction'] == cluster, 'neighbors_list'].to_list()[0]
            for neighbour in L:
                # check if the neighbour is sufficiently heavy
                if not df_lst[k].loc[df_lst[k]['junction'] == neighbour,'small_junction'].item():
                    df_cluster = df_lst[k].loc[df_lst[k]['junction'] == cluster].copy()
                    df_neighbour = df_lst[k].loc[df_lst[k]['junction'] == neighbour].copy()
                    # return df_cluster, df_neighbour
                    arw = Arrow3D([df_cluster['x'].item(),df_neighbour['x'].item()],
                              [df_cluster['z'].item(),df_neighbour['z'].item()],
                              [df_cluster['y'].item(),df_neighbour['y'].item()],
                              color="grey",
                              lw = 1, mutation_scale=0.5,alpha=0.95)
                    ax.add_artist(arw)
                else:
                    pass
        # else:
        #     pass
    # L = [[5,0,0,'black'],[0,5,0,'green'],[0,0,5,'orange']]
    # df = pd.DataFrame(L,columns=['x','y','z','c'])
    # sc = ax.scatter(df['x'],df['y'],df['z'],
    #                 c = df['c'],
    #                 alpha=1)
    ax.view_init(0,0)
    
    return sc

def graphMass(dir_in, lst, loop = True, SJ = False):
    """
    Return a seaborn plot showing the mass of the junction with or without the
    loop counted and with or without the small junctions.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the directory named name.
    lst: list of strings
        List of subdirectories in the input directory in which the files lie.
    loop: boolean, default True
        Indicate if you want the loop counted or not.
    SJ: boolean, default False
        Indicate if you want the small junctions counted or not.

    Returns:
    --------
    figure: seaborn violinplot
        The plot showing the distribution of the mass of the junctions in the 
        network.
    """
    
    L = []
    for elem in lst:
        data_list = dfListConnex(dir_in, elem)
        for j in range(len(data_list)):
            aff = data_list[j].iloc[0]['affinity']
            M = []
            df = includeSJ(data_list[j], SJ)
            for k in range(len(df)):
                if loop:
                    M += [[aff, df.iloc[k]['mass']]]
                else:
                    M += [[aff, df.iloc[k]['degreeMG']]]
            L += M
    
    df = pd.DataFrame(L,columns=['aff','mass'])
    T = list(set(df['aff'].to_list()))
    if (5.5 in T) or (6.5 in T):
        for elem in [0.5*k for k in range(2,14)]:
            if elem not in T:
                df =  df.append({'aff': elem,
                                 'mass': -1},
                                ignore_index=True)
    df['affinity'] = df['aff'].apply(lambda x: round(1 - (x/25),2))
    figure = sns.violinplot(x="affinity",y='mass',data=df)
    
    for tick in figure.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)

    for tick in figure.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    
    plt.xlabel('Affinity', size=28, family='helvetica')
    if loop:
        plt.ylabel('Junction mass with loops included', size=28, family='helvetica')
    else:
        plt.ylabel('Junction mass without loops', size=28, family='helvetica')
    figure.set(ylim=(0, None))
    return figure


def graphDegree(dir_in, lst, multi = True, SJ = False):
    """
    Return a seaborn plot showing the degree or multoidegree of the junction 
    with or without the loop counted and with or without the small junctions.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the directory named name.
    lst: list of strings
        List of subdirectories in the input directory in which the files lie.
    multi: boolean, default True
        Indicate if you want the degree or the multidegree.
    SJ: boolean, default False
        Indicate if you want the small junctions counted or not. The degree 
        (or multidegree) of a small junction is always 0.

    Returns:
    --------
    figure: seaborn violinplot
        The plot showing the distribution of the degree or multidegree of the 
        junctions in the network.
    """
    
    
    L = []
    for elem in lst:
        data_list = dfListConnex(dir_in, elem)
        for j in range(len(data_list)):
            aff = data_list[j].iloc[0]['affinity']
            M = []
            df = includeSJ(data_list[j], SJ)
            for k in range(len(df)):
                if multi:
                    M += [[aff, df.iloc[k]['degreeMG']]]
                else:
                    M += [[aff, df.iloc[k]['degreeG']]]
            L += M
    
    df = pd.DataFrame(L,columns=['aff','degree'])
    T = list(set(df['aff'].to_list()))
    if (5.5 in T) or (6.5 in T):
        for elem in [0.5*k for k in range(2,14)]:
            if elem not in T:
                df =  df.append({'aff': elem,
                                 'degree': -1},
                                ignore_index=True)
    df['affinity'] = df['aff'].apply(lambda x: round(1 - (x/25),2))
    
    figure = sns.violinplot(x="affinity",y='degree',data=df)
    
    for tick in figure.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)

    for tick in figure.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    
    plt.xlabel('Affinity', size=28, family='helvetica')
    if multi:
        plt.ylabel('Multidegree', size=28, family='helvetica')
    else:
        plt.ylabel('Degree', size=28, family='helvetica')
    figure.set(ylim=(0, None))
    return figure


def graphClustering(dir_in, lst, unweighted = True, SJ = False):
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
    
    df_clustering = pd.DataFrame()
    df_random = pd.DataFrame()
    
    affinityList = []
    typeList = []
    
    affinityRandomList = []
    CCList = []
    CCRandomList = []
    runList = []
    runRandomList = []
    if unweighted:
        for name in lst:
            df_lst = dfListConnex(dir_in, name)
            n = len(df_lst)
            try:
                affinityList += [df_lst[0].at[0,'affinity']]*n
            except KeyError:
                print(name)
                break
            typeList += [str(df_lst[0].at[0,'type'])]*n
            affinityRandomList += ['small world']*n
            #label = df_lst[0].at[0,'label']
            for k in range(len(df_lst)):
                df_edges = toEdgeList(includeSJ(df_lst[k], SJ))
                G, avgUWCC = unweightedCC(df_edges)
                
                # Get the parameters
                n = len(df_lst[k])
                L = list(G.degree())
                M = [x[1] for x in L]
                m = sum(M)/2
               
                # Create a small-world graph and exploit
                k = int(np.ceil((m//n) + 1))
                K = nx.connected_watts_strogatz_graph(n, k, 0.95, tries=100, seed=None)
                
                avgSW = nx.average_clustering(K)
                
                CCList += [avgUWCC]
                CCRandomList += [avgSW]
                runList += [df_lst[k].at[0,'label']]
                runRandomList += ['small world']
    else:
        for name in lst:
            df_lst = dfListConnex(dir_in, name)
            n = len(df_lst)
            try:
                affinityList += [df_lst[0].at[0,'affinity']]*n
            except KeyError:
                print(name)
                break
            typeList += [str(df_lst[0].at[0,'type'])]*n
            # affinityRandomList += ['small world']*n
            for k in range(len(df_lst)):
                df_edges = toEdgeList(includeSJ(df_lst[k], SJ))
                G, avgWCC = weightedCC(df_edges)
                
                # # Get the parameters
                # n = len(df_lst[k])
                # L = list(G.degree())
                # M = [x[1] for x in L]
                # m = sum(M)/2
               
                # # Create a small-world graph and exploit
                # k = int(np.ceil((m//n) + 1))
                # K = nx.connected_watts_strogatz_graph(n, k, 0.95, tries=100, seed=None)
                
                # degSW = nx.degree(K)
                # avgSW = 0
                # for key in K.nodes():
                #     avgSW += degSW[key]
                # avgSW /= len(degSW)
                
                CCList += [avgWCC]
                # CCRandomList += [avgSW]
                runList += [name]
                # runRandomList += ['small world']
            

    df_random['type'] = affinityRandomList
    df_random['affinity'] = affinityRandomList
    df_random['run'] = runRandomList
    df_random['CC'] = CCRandomList
    
    # newRunList = []
    # for elem in runList:
    #     if elem == '6E48/14001':
    #         newRunList += ['A']
    #     elif elem == '6E48/14002':
    #         newRunList += ['B']
    #     elif elem == '6E48/14003':
    #         newRunList += ['C']
    #     elif elem == '6E48/14004':
    #         newRunList += ['D']
    #     elif elem == '7E48/14001':
    #         newRunList += ['A two off']
    #     elif elem == '7E48/14002':
    #         newRunList += ['B two off']
    #     elif elem == '7E48/14003':
    #         newRunList += ['C two off']
    #     elif elem == '7E48/14004':
    #         newRunList += ['D two off']
    
    
    df_clustering['type'] = typeList
    df_clustering['aff'] = affinityList
    df_clustering['run'] = runList
    df_clustering['CC'] = CCList

    T = list(set(df_clustering['aff'].to_list()))
    if (5.5 in T) or (6.5 in T):
        for elem in [0.5*k for k in range(2,14)]:
            if elem not in T:
                df_clustering =  df_clustering.append({'aff': elem,
                                                        'CC': -1},
                                                      ignore_index=True)
    df_clustering['affinity'] = df_clustering['aff'].apply(lambda x: round(1 - (float(x)/25),2))
   
    #df_clustering = df_clustering.sort_values(by=['run'])
    
    df_random = df_random.append(df_clustering,ignore_index=True)
    sns.set_style({'axes.grid': True})
    figure = sns.violinplot(x = 'run', y = 'CC', 
                     scale = "width", 
                     data = df_clustering
                     )
    
    # Set the graphical properties
    for tick in figure.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in figure.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    
    # Labels for the polymer type and concentration graphs
    # ax.set_xticks(range(len(df_random['run'].unique()))) # <--- set the ticks first
    # ax.set_xticklabels(['small world',
    #                      '2B48 - 0.0008', '0.001', '0.002', '0.004', '0.006', '0.008',
    #                      '4B16 - 0.0005', '0.001', '0.002', '0.004',
    #                      '5B8 - 0.0004', '0.0006', '0.0008', '0.001',
    #                      '6B6 - 0.0004', '0.0006', '0.0008', '0.001'
    #                      ])
    plt.xticks(rotation=19, horizontalalignment='right')
        
    plt.xlabel('Polymer Type and concentration', size=28, family='helvetica')

    if unweighted:
        plt.ylabel('Average clustering coefficient', size=28, family='helvetica')
        figure.set(ylim=(0,1))
    else:
        plt.ylabel('Average weighted clustering coefficient', size=28, family='helvetica')
        figure.set(ylim=(0,None))
        

def graphClusteringComparison(dir_in, lst, SJ = False):
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
    
    df_clusteringCC = pd.DataFrame()
    df_random = pd.DataFrame()
    
    CCaffinityList = []
    CCtypeList = []
    
    affinityRandomList = []
    CCList = []
    CCRandomList = []
    CCrunList = []
    runRandomList = []
 
    for name in lst:
        df_lst = dfListConnex(dir_in, name)
        n = len(df_lst)
        try:
            CCaffinityList += [df_lst[0].at[0,'affinity']]*n
        except KeyError:
            print(name)
            break
        CCtypeList += [str(df_lst[0].at[0,'type'])]*n
        affinityRandomList += ['small world']*n
        for k in range(len(df_lst)):
            df_edges = toEdgeList(includeSJ(df_lst[k], SJ))
            G, avgUWCC = unweightedCC(df_edges)
            
            # Get the parameters
            n = len(df_lst[k])
            L = list(G.degree())
            M = [x[1] for x in L]
            m = sum(M)/2
           
            # Create a small-world graph and exploit
            k = int(np.ceil((m//n) + 1))
            K = nx.connected_watts_strogatz_graph(n, k, 0.95, tries=100, seed=None)
            
            avgSW = nx.average_clustering(K)
            
            CCList += [avgUWCC]
            CCRandomList += [avgSW]
            CCrunList += [name]
            runRandomList += ['small world']
  

    df_random['type'] = affinityRandomList
    df_random['affinity'] = affinityRandomList
    df_random['run'] = runRandomList
    df_random['CC'] = CCRandomList
    
    df_clusteringCC['type'] = CCtypeList
    df_clusteringCC['aff'] = CCaffinityList
    df_clusteringCC['run'] = CCrunList
    df_clusteringCC['CC'] = CCList

    T = list(set(df_clusteringCC['aff'].to_list()))
    if (5.5 in T) or (6.5 in T):
        for elem in [0.5*k for k in range(2,14)]:
            if elem not in T:
                df_clusteringCC =  df_clusteringCC.append({'aff': elem,
                                                        'CC': -1},
                                                      ignore_index=True)
    df_clusteringCC['affinity'] = df_clusteringCC['aff'].apply(lambda x: round(1 - (float(x)/25),2))
   
    df_clusteringCC = df_clusteringCC.sort_values(by=['type'])
    
    df_random = df_random.append(df_clusteringCC,ignore_index=True)
    
    f, (subfigure1, subfigure2) = plt.subplots(1, 2, gridspec_kw={'width_ratios': [5, 4]})
    sns.set_style({'axes.grid': True})
    plt.rcParams["axes.titlesize"] = 25
    sns.violinplot(x = 'type', y = 'CC', 
                     scale = "width", 
                     data = df_random,
                     ax = subfigure1
                     )

    
    # Set the graphical properties
    for tick in subfigure1.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in subfigure1.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    
    # WCC
    df_clusteringWCC = pd.DataFrame()
    
    WCCaffinityList = []
    WCCtypeList = []
    
    affinityRandomList = []
    WCCList = []
    WCCrunList = []

    
    for name in lst:
        df_lst = dfListConnex(dir_in, name)
        n = len(df_lst)
        try:
            WCCaffinityList += [df_lst[0].at[0,'affinity']]*n
        except KeyError:
            print(name)
            break
        WCCtypeList += [str(df_lst[0].at[0,'type'])]*n
        for k in range(len(df_lst)):
            df_edges = toEdgeList(includeSJ(df_lst[k], SJ))
            G, avgWCC = weightedCC(df_edges)
            
            
            WCCList += [avgWCC]
            WCCrunList += [name]
    
    df_clusteringWCC['type'] = WCCtypeList
    df_clusteringWCC['aff'] = WCCaffinityList
    df_clusteringWCC['run'] = WCCrunList
    df_clusteringWCC['CC'] = WCCList

    T = list(set(df_clusteringWCC['aff'].to_list()))
    if (5.5 in T) or (6.5 in T):
        for elem in [0.5*k for k in range(2,14)]:
            if elem not in T:
                df_clusteringWCC =  df_clusteringWCC.append({'aff': elem,
                                                        'CC': -1},
                                                      ignore_index=True)
    df_clusteringWCC['affinity'] = df_clusteringWCC['aff'].apply(lambda x: round(1 - (float(x)/25),2))
   
    df_clusteringWCC = df_clusteringWCC.sort_values(by=['type'])

    sns.set_style({'axes.grid': True})
    sns.violinplot(x = 'type', y = 'CC', 
                     scale = "width", 
                     data = df_clusteringWCC,
                     ax = subfigure2
                     )

    
    # Set the graphical properties
    for tick in subfigure2.xaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    for tick in subfigure2.yaxis.get_major_ticks():
        tick.label.set_fontsize(20)
    
    # Labels for the polymer type and concentration graphs
    # ax.set_xticks(range(len(df_random['run'].unique()))) # <--- set the ticks first
    # ax.set_xticklabels(['small world',
    #                      '2B48 - 0.0008', '0.001', '0.002', '0.004', '0.006', '0.008',
    #                      '4B16 - 0.0005', '0.001', '0.002', '0.004',
    #                      '5B8 - 0.0004', '0.0006', '0.0008', '0.001',
    #                      '6B6 - 0.0004', '0.0006', '0.0008', '0.001'
    #                      ])
    # plt.xticks(rotation=19, horizontalalignment='right')
    
    subfigure1.set_xlabel('Polymer Type', size=28, family='helvetica')
    subfigure1.set_ylabel('Average clustering coefficient', size=28, family='helvetica')
    subfigure1.set(ylim=(0,1))
        
    subfigure2.set_xlabel('Polymer Type', size=28, family='helvetica')
    subfigure2.set_ylabel('Average weighted clustering coefficient', size=28, family='helvetica')
    subfigure2.set(ylim=(0,None))
    
   
    
    return f.tight_layout()


def graphFluidity(dir_in, filename):
    """
    Return a seaborn plot showing the mass of the junction with or without the
    loop counted and with or without the small junctions.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the file filename.
    filename: string
        Name of the file that will de converted in a dataframe and where we 
        already have saved the fluidity (it's a time consuming 'function thats
        why it's already done).

    Returns:
    --------
    figure: seaborn lineplot
        The plot showing the evolution of the average fluidity of the network.
    """
    
    dir_path = Path(os.getcwd()) 
    file_path = dir_path.parent.joinpath(dir_in, filename)
    df = pd.read_pickle(file_path)
    
    figure = sns.lineplot(x = 'step', y = 'fluidity', 
                     hue = "affinity", 
                     data = df
                     )
    
    return figure


# =============================================================================
## List to test
# =============================================================================


P = ['1011', '2414', '2450', '622', '2214', '2806']
Q = ['2B/6091', '2B/6092', '2B/6093', '2B/6094', '2B/6095', '2B/6096',
     '4C/4025', '4C/4026', '4C/4027', '4C/4028',
     '5D/11001', '5D/11002', '5D/11003', '5D/11004',
     '6E/14001', '6E/14002', '6E/14003', '6E/14004']
R = ['2B/6091',
     '4C/4025', 
     '5D/11003',
     '6E/14003']
S = ['6E48/14001','7E48/14001',
     '6E48/14002','7E48/14002',
     '6E48/14003','7E48/14003',
     '6E48/14004','7E48/14004']
T = ['4/6e48/14001', '4/6e48/14002', '4/6e48/14003', '4/6e48/14004',
     '4/6f48/15001', '4/6f48/15002', '4/6f48/15003', '4/6f48/15004',
     #'65/6d48/13001r6', '65/6d48/13002r6', '65/6d48/13003r6', '65/6d48/13004r6',
     '65/6e48/14001', '65/6e48/14002', '65/6e48/14003', '65/6e48/14004',
     '65/6d48/13001', '65/6d48/13002', '65/6d48/13003', '65/6d48/13004',
     '65/6f48/15001', '65/6f48/15002', '65/6f48/15003', '65/6f48/15004'
     ]

# =============================================================================
## Graphs generations
# =============================================================================

graphClustering('neo_junctions',T[:8])

