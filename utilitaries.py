#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  3 10:04:19 2021

@author: clement
"""

# =============================================================================
## Library dependancies
# =============================================================================

import pandas as pd
import numpy as np
import os
from pathlib import Path
from os.path import isfile, join

# =============================================================================
## Utilitary functions
# =============================================================================


def simpleToEdgeList(data, links):
    """
    Return a new dataframe containing the edges spanned by the junctions and
    their neighbors, multiple edges are returned if there is multiple links 
    between them.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least two columns named junction and neighbors.   

    Returns:
    --------
    df: pandas dataframe
        A new dataframe containing the edges spanned by the junctions and
        their neighbors, multiple edges are returned if there is multiple links 
        between them.
    """
    
    L = []
    
    for i in range(len(data)):
        for j in range(len(data.at[i,'neighbors'])):
            if ([data.at[i,'junction'], data.at[i,'neighbors'][j]] in L) or ([data.at[i,'neighbors'][j],data.at[i,'junction']] in L):
                pass
            else:
                L += [[data.at[i,'junction'], data.at[i,'neighbors'][j]]]*links

    df = pd.DataFrame(L,columns=['source', 'target'])
    
    return df

def toEdgeList(data):
    """
    Return a new dataframe containing the edges spanned by the junctions and
    their neighbors, multiple edges are returned if there is multiple links 
    between them.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least two columns named junction and neighbors.   

    Returns:
    --------
    df: pandas dataframe
        A new dataframe containing the edges spanned by the junctions and
        their neighbors, multiple edges are returned if there is multiple links 
        between them.
    """
    
    L = []
    
    for i in range(len(data)):
        for j in range(len(data.at[i,'neighbors'])):
            if ([data.at[i,'junction'], data.at[i,'neighbors'][j][0], data.at[i,'neighbors'][j][2]] in L) or ([data.at[i,'neighbors'][j][0],data.at[i,'junction'], data.at[i,'neighbors'][j][2]] in L):
                pass
            else:
                L += [[data.at[i,'junction'], data.at[i,'neighbors'][j][0], data.at[i,'neighbors'][j][2]]]*data.at[i,'neighbors'][j][2]
            #print(L)

    df = pd.DataFrame(L,columns=['source', 'target', 'links'])
    
    return df


def dfListConnex(dir_in, name):
    """
    Return a list of dataframes containing only the main connex component 
    (containing the most junctions).

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the directory named name.
    name: string
        Subdirectory where the script will take the files.

    Returns:
    --------
    df_list: list of pandas dataframes
        A list of dataframes containing only the main connex component
        (containing the most junctions).
    """
    
    dir_path = Path(os.getcwd()) 
    path = dir_path.parent.joinpath(dir_in,name)
    
    # work on each files
    onlyfiles = [f for f in os.listdir(path) if isfile(join(path, f))]
    
    file_list = []
    for filename in onlyfiles:
        if filename.endswith('.pkl'):
            file_list.append(filename)
    file_list.sort()
    
    df_list = []
    for i in range(len(file_list)):
        file_path = os.path.join(path,file_list[i])
        df = pd.read_pickle(file_path)
        
        n = len(df)
        L = df['component'].to_list()
        L_max = max(L) + 1
        component = -1
        component_members = 0
        for j in range(L_max):
            if L.count(j) > np.ceil(n/2):
                component = j
                break
            else:
                if L.count(j) > component_members:
                    component_members = L.count(j)
                    component = j
                else:
                    pass
        
        df_connex = df.loc[df['component'] == int(component)].reset_index()
        df_list += [df_connex]
    
    return df_list


def includeSJ(data, y = True):
    """
    Return the same dataframe with or without the small junctions (mass < 3),
    but keeping the bridges (junctions with mass 2 connecting two junctions
    of mass > 2).

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least one column named neighbors.
    y: bool, default True
        Indicate if you include the small junctions or not.

    Returns:
    --------
    data: pandas dataframe
        The input dataframe with or without the small junctions.
    """
    
    if y:
        return data
    else:
        data2 = data.loc[data['small_junction'] == 0]
        data3 = data.loc[(data['small_junction'] == 1) & (data['bridge'] == 1)]
        
        data4 = pd.concat([data2, data3], ignore_index=True)
        data4.sort_values(by=['junction'],inplace=True)
        data4.reset_index(inplace=True)
        data4.drop('index',axis=1,inplace=True)
        
        return data4