#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  6 11:26:20 2021

@author: clement
"""

# =============================================================================
## Library dependancies
# =============================================================================

import pandas as pd
import os
import seaborn as sns
from pathlib import Path
from os.path import isfile, isdir, join

from utilitaries import dfListConnex

# =============================================================================
## Basic functions
# =============================================================================


def children(data_list, max_step = 6):
    """
    Return a dataframe containing the children of each node for all steps 
    in max_step.

    Parameters:
    -----------
    data_list: list of pandas dataframe
        Each element in data_list must contain at least three columns named 
        members, weight and small_junction.
    max_steps: interger, default 6.
        Number of time steps the children search will de done.

    Returns:
    --------
    df: pandas dataframe
        A new dataframe with the weight, the step, the children and children 
        number for each cluster and each time step.
    """
    
    maxStep = min(len(data_list), max_step)
    if maxStep < max_step:
        print("Not enough data to reach {} steps. Aborting operation.".format(max_step))
    n = len(data_list[0])
    L_len = []
    A = []
    for i in range(n):
        # Create a list of all the members in the junction
        L = data_list[0].iloc[i]['members']
        L_len += [data_list[0].iloc[i]['mass']]*maxStep
        for step in range(1,1+maxStep):
            M = []
            # Check if the members of the cluster in the dataframe corresponding to step are in L
            try:
                #print('length: {}'.format(len(data_list[step])))
                for k in range(len(data_list[step])):
                    for elem in data_list[step].at[k,'members']:
                        if elem in L:
                            M += [(k, data_list[step].at[k,'small_junction'])]
            except IndexError:
                print(k)
            # Eliminate the repetition in M
            N = sorted(list(set(M)))
            P = []
            # Count the occurences of each children
            for elem in N:
                P += [(elem[0], elem[1], M.count(elem))]
            A += [[step, P, len(P)]]
    
    df = pd.DataFrame(A, columns=['step', 'children', 'children_number'])
    df['mass'] = L_len
    
    return df


# def childrenAlternative(data_list, max_step = 6):
#     """
#     Return a dataframe containing the children of each node for all steps 
#     in max_step.

#     Parameters:
#     -----------
#     data_list: list of pandas dataframe
#         Each element in data_list must contain at least three columns named 
#         members, weight and small_junction.
#     max_steps: interger, default 6.
#         Number of time steps the children search will de done.

#     Returns:
#     --------
#     df: pandas dataframe
#         A new dataframe with the weight, the step, the children and children 
#         number for each cluster and each time step.
#     """
    
#     maxStep = min(len(data_list), max_step)
#     if maxStep < max_step:
#         print("Not enough data to reach {} steps. Aborting operation.".format(max_step))
#     n = len(data_list[0])
#     L_len = []
#     A = []
#     for i in range(n):
#         # Create a list of all the members in the junction
#         L = data_list[0].iloc[i]['members']
#         L_len += [data_list[0].iloc[i]['mass']]*maxStep
#         for step in range(1,1+maxStep):
#             M = []
#             # Check if the members of the cluster in the dataframe corresponding to step are in L
#             for elem in L:
#                 k = data_list[step].loc[elem in data_list[step]['members']]
#             # for k in range(len(data_list[step])):
#             #     for elem in data_list[step].at[k,'members']:
#             #         if elem in L:
#             #             M += [(k, data_list[step].at[k,'small_junction'])]
#             # Eliminate the repetition in M
#             N = sorted(list(set(M)))
#             P = []
#             # Count the occurences of each children
#             for elem in N:
#                 P += [(elem[0], elem[1], M.count(elem))]
#             A += [[step, P, len(P)]]
    
#     df = pd.DataFrame(A, columns=['step', 'children', 'children_number'])
#     df['mass'] = L_len
    
#     return df


def fluidity(data_list, max_offset, max_step):
    """
    Return a dataframe containing the children of each node for all steps 
    in max_step and for all starting data in max_offset.

    Parameters:
    -----------
    data_list: list of pandas dataframe
        Each element in data_list must contain at least three columns named 
        members, weight and small_junction.
    max_offset: integer
        Until which file in the data_list does the function make the children 
        search?
    max_steps: interger
        Number of time steps the children search will de done.

    Returns:
    --------
    df: pandas dataframe
        A new dataframe with the weight, the step, the children, children 
        number and fluidity for each cluster and each time step.
    """
    
    dfluidity = []
    
    for offset in range(max_offset):
        data = children(data_list[offset:], max_step)
        # Add the starting time to teh dataframe
        data['time'] = [int(50000*(offset+1)) for j in range(len(data))]
        K = []
        # Add the fluidity measure to each junction
        for i in range(len(data)):
            K += [(data.at[i,'children_number'] - 1)/ data.at[i,'mass']]
        data['fluidity'] = K
        dfluidity += [data]
    
    df = pd.concat(dfluidity, ignore_index=True)
    return df


# =============================================================================
## Main function
# =============================================================================


def mainFluidity(dir_in, dir_out, save_name, max_offset, max_step):
    """
    Save a dataframe containing the fluidity of each node for a range of time 
    steps and strating points.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search take all the subdirectories.
    dir_out: string
        Output directory, that the script will eventually create and where it 
        will save the outputed dataframe.
    max_offset: integer
        Until which file in the data_list does the function make the children 
        search?
    max_steps: interger
        Number of time steps the children search will de done.

    Returns:
    --------
    df: pandas dataframe
        A new dataframe with the weight, the step, the children, children 
        number and fluidity for each cluster and each time step for all the 
        runs in the input directory.
    """
    
    dir_path = Path(os.getcwd()) 
    path = dir_path.parent.joinpath(dir_in)
    save_path = dir_path.parent.joinpath(dir_out)
    save_path.mkdir(parents=True, exist_ok=True)
    
    dfluidity = []
    
    onlydirs = [d for d in os.listdir(path) if isdir(join(path, d))]
    n = len(onlydirs)
    count = 1
    
    for name in onlydirs:
        print("file {:02d} / {}".format(count,n))
        subpath = path.joinpath(name)
        # work on each files
        onlyfiles = [f for f in os.listdir(subpath) if isfile(join(subpath, f))]
        
        file_list = []
        for filename in onlyfiles:
            if filename.endswith('.pkl'):
                file_list.append(filename)
        file_list.sort()

        df_list = []
        for i in range(len(file_list)):
            
            file_path = os.path.join(subpath,file_list[i])
            df_list += [pd.read_pickle(file_path)]
        
        df1 = fluidity(df_list, max_offset, max_step)
        df1['run'] = [name] * len(df1)
        #df1['affinity'] = [df_list.iloc[0]['affinity']] * len(df1)
        dfluidity += [df1]
        
        count += 1
    
    df = pd.concat(dfluidity, ignore_index=True)
    saveName = save_path.joinpath(save_name + '.pkl')
    df.to_pickle(saveName)