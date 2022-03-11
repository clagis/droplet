#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 15:58:58 2021

@author: clement
"""

# =============================================================================
## Library dependancies
# =============================================================================

import os
import pandas as pd
from pathlib import Path
from os.path import isfile, isdir, join

# =============================================================================
## Main functions
# =============================================================================


def affinity(dir_in):
    """
    Add affinity to the dataframes in the given directory and its subdiretories.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the directory named name.

    Returns:
    --------
    Nothing directly. The script will save the files in the same directories 
    it found them but with the affinity column added.
    """
    
    L = []
    
    dir_path = Path(os.getcwd()) 
    path = dir_path.parent.joinpath(dir_in)
    
    onlydirs = [d for d in os.listdir(path) if isdir(join(path, d))]
    
    for name in onlydirs:
        subpath = path.joinpath(name)
        # work on each files
        onlyfiles = [f for f in os.listdir(subpath) if isfile(join(subpath, f))]
        
        file_list = []
        for filename in onlyfiles:
            if filename.endswith('.pkl'):
                file_list.append(filename)
        file_list.sort()
        
        for file in file_list:
            df = pd.read_pickle(subpath.joinpath(file))
            n = len(df)
            if name in ['13001', '13002', '13003', '13004', '13001r6', '13002r6', '13003r6', '13004r6','14001', '14002', '14003', '14004','15001', '15002', '15003', '15004']:
                df['affinityParameter'] = [4 for i in range(n)]
                df['affinity'] = df['affinityParameter'].apply(lambda x: 1 - (x/25))
                df.to_pickle(subpath.joinpath(file))
            # if name in ['4025','4026','4027','4028', '11001', '11002', '11003', '11004', '14001', '14002', '14003', '14004','6091', '6092', '6093', '6094', '6095', '6096']:
            #     df['affinityParameter'] = [4 for i in range(n)]
            #     ['1008','1009','1010','1011']:
            #     df['affinity'] = [1 for i in range(n)]
            # # elif name in ['8022']:
            #     df['affinity'] = [2 for i in range(n)]
            # elif name in ['2412','2414','8016']:
            #     df['affinity'] = [3 for i in range(n)]
            # elif name in ['2410', '2411', '2412', '2413', '2414',
            #               '2415', '2416', '2417', '2418',
            #               '2448', '2450', '8010', '4026', '4027',
            #               'd4_13001', 'd4_13002', 'd4_13003', 'd4_13004',
            #               'e4_14001', 'e4_14002', 'e4_14003', 'e4_14004',
            #               'f4_15001', 'f4_15002', 'f4_15003', 'f4_15004',
            #               '7e_14001', '7e_14002', '7e_14003', '7e_14004',
            #               '7f_15001', '7f_15002', '7f_15003', '7f_15004']:
            #     df['affinity'] = [4 for i in range(n)]
            # elif name in ['5D', '6E', '2B_b', '4C'
            #               '2446', '2447', '2448', '2449',
            #               '2451', '2452', '2453', '2454', 
            #               '620','621','622','623','624','625','626','8004','8005',
            #               '4030','4031','6001','6002','6003','6004','704','705',
            #               '706','707',
            #               'd5_13001', 'd5_13002', 'd5_13003', 'd5_13004',
            #               'e5_14001', 'e5_14002', 'e5_14003', 'e5_14004',
            #               'f5_15001', 'f5_15002', 'f5_15003', 'f5_15004']:
            #     df['affinity'] = [5 for i in range(n)]
            # elif name in ['4034','4035']:
            #     df['affinity'] = [5.5 for i in range(n)]
            # elif name in ['2214','2216','3002','4002','4003',
            #               'd6_13001', 'd6_13002', 'd6_13003', 'd6_13004',
            #               'e6_14001', 'e6_14002', 'e6_14003', 'e6_14004',
            #               'f6_15001', 'f6_15002', 'f6_15003', 'f6_15004']:
            #     df['affinity'] = [6 for i in range(n)]
            # elif name in ['2806','2808','3010','4006','4007','6028']:
            #     df['affinity'] = [6.5 for i in range(n)]
            # elif name in ['4010','3018','4011']:
            #     df['affinity'] = [7 for i in range(n)]
            else:
                # L += [name]
                pass
            
            
    
    return list(set(L))


def polymerType(dir_in):
    """
    Add the type of the polymer to the dataframes in the given directory and 
    its subdiretories.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the directory named name.

    Returns:
    --------
    Nothing directly. The script will save the files in the same directories 
    it found them but with the (polymer) type column added.
    """
    
    L = []
    
    dir_path = Path(os.getcwd()) 
    path = dir_path.parent.joinpath(dir_in)
    
    onlydirs = [d for d in os.listdir(path) if isdir(join(path, d))]
    
    for name in onlydirs:
        subpath = path.joinpath(name)
        # work on each files
        onlyfiles = [f for f in os.listdir(subpath) if isfile(join(subpath, f))]
        
        file_list = []
        for filename in onlyfiles:
            if filename.endswith('.pkl'):
                file_list.append(filename)
        file_list.sort()
        
        for file in file_list:
            df = pd.read_pickle(subpath.joinpath(file))
            n = len(df)
            if name in ['1011','622','2450','2214','2806','2414']:
                df['type'] = ['2B16' for i in range(n)]
            elif name in ['6091', '6092', '6093', '6094', '6095', '6096']:
                df['type'] = ['2B48' for i in range(n)]
            elif name in ['4025','4026','4027','4028']:
                df['type'] = ['4B16' for i in range(n)]
            elif name in ['11001', '11002', '11003', '11004']:
                df['type'] = ['5B8' for i in range(n)]
            elif name in ['13001', '13002', '13003', '13004', '13001r6', '13002r6', '13003r6', '13004r6','14001', '14002', '14003', '14004','15001', '15002', '15003', '15004']:
                df['type'] = ['6B6' for i in range(n)]
        
            df.to_pickle(subpath.joinpath(file))
    
    return list(set(L))



def label(dir_in):
    """
    Add the type of the polymer to the dataframes in the given directory and 
    its subdiretories.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the directory named name.

    Returns:
    --------
    Nothing directly. The script will save the files in the same directories 
    it found them but with the (polymer) type column added.
    """
    
    L = []
    
    dir_path = Path(os.getcwd()) 
    path = dir_path.parent.joinpath(dir_in)
    
    onlydirs1 = [d for d in os.listdir(path) if isdir(join(path, d))]
    
    for name1 in onlydirs1:
        subpath = path.joinpath(name1)
        # work on each files
        
        onlydirs2 = [d for d in os.listdir(subpath) if isdir(join(subpath, d))]
        
        for name2 in onlydirs2:
            subpath2 = subpath.joinpath(name2)
            onlyfiles = [f for f in os.listdir(subpath2) if isfile(join(subpath2, f))]
            
            file_list = []
            for filename in onlyfiles:
                if filename.endswith('.pkl'):
                    file_list.append(filename)
            file_list.sort()
            
            for file in file_list:
                df = pd.read_pickle(subpath2.joinpath(file))
                n = len(df)
                if name1 == '6e48':
                    if name2 == '14001':
                        df['label'] = ['6B6 A' for i in range(n)]
                    if name2 == '14002':
                        df['label'] = ['6B6 B' for i in range(n)]
                    if name2 == '14003':
                        df['label'] = ['6B6 C' for i in range(n)]
                    if name2 == '14004':
                        df['label'] = ['6B6 D' for i in range(n)]
                elif name1 == '6d48':
                    if name2 == '13001':
                        df['label'] = ['6B8 A' for i in range(n)]
                    if name2 == '13002':
                        df['label'] = ['6B8 B' for i in range(n)]
                    if name2 == '13003':
                        df['label'] = ['6B8 C' for i in range(n)]
                    if name2 == '13004':
                        df['label'] = ['6B8 D' for i in range(n)]
                elif name1 == '6f48':
                     if name2 == '15001':
                         df['label'] = ['6B10 A' for i in range(n)]
                     if name2 == '15002':
                         df['label'] = ['6B10 B' for i in range(n)]
                     if name2 == '15003':
                         df['label'] = ['6B10 C' for i in range(n)]
                     if name2 == '15004':
                         df['label'] = ['6B10 D' for i in range(n)]
            
                df.to_pickle(subpath2.joinpath(file))
        
    return list(set(L))