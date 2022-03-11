#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 10:48:45 2022

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


def graphLabels(dir_in, kind):
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
        print(name)
        subpath = path.joinpath(name)
        # work on each files
        onlyfiles = [f for f in os.listdir(subpath) if isfile(join(subpath, f))]
        
        file_list = []
        for filename in onlyfiles:
            if filename.endswith('.pkl'):
                file_list.append(filename)
        file_list.sort()
        if kind == '6E48':
            for file in file_list:
                df = pd.read_pickle(subpath.joinpath(file))
                n = len(df)
                if name in ['14001']:
                    df['label'] = ['A' for i in range(n)]
                    df.to_pickle(subpath.joinpath(file))
                
                elif name in ['14002']:
                    df['label'] = ['B' for i in range(n)]
                    df.to_pickle(subpath.joinpath(file))
               
                elif name in ['14003']:
                    df['label'] = ['C' for i in range(n)]
                    df.to_pickle(subpath.joinpath(file))
               
                elif name in ['14004']:
                    df['label'] = ['D' for i in range(n)]
                    df.to_pickle(subpath.joinpath(file))
                else:
                    # L += [name]
                    pass
        elif kind == '7E48':
            for file in file_list:
                df = pd.read_pickle(subpath.joinpath(file))
                n = len(df)
                if name in ['14001']:
                    df['label'] = ['A two off' for i in range(n)]
                    df.to_pickle(subpath.joinpath(file))
                elif name in ['14002']:
                    df['label'] = ['B two off' for i in range(n)]
                    df.to_pickle(subpath.joinpath(file))
                elif name in ['14003']:
                    df['label'] = ['C two off' for i in range(n)]
                    df.to_pickle(subpath.joinpath(file))
                elif name in ['14004']:
                    df['label'] = ['D two off' for i in range(n)]
                    df.to_pickle(subpath.joinpath(file))
                else:
                    # L += [name]
                    pass
        else:
            # L += [name]
            pass
    return list(set(L))