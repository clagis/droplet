#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 13:31:37 2021

@author: clement
"""

# =============================================================================
## Library dependancies
# =============================================================================

import pandas as pd
import os
import numpy as np
from pathlib import Path
from os.path import isfile, join

# =============================================================================
## Basic functions
# =============================================================================


def dist_matrix(data, box_size = 48, boundary = 'pbc'):
    """
    Computes the distance between every points in a 3D box with periodic or non perciodic boundary.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least three columns named x, y and z.
    box_size: int
        Size of the simulation box.
    boundary: 'pbc' or 'npbc', default 'pbc'
        Type of boundary.

    Returns:
    --------
    D: numpy array
        Distance matrix between every point of the data, in truth only the upper
        triangular part of the matrix is returned. Uncomment N if you wish the
        full distance matrix.
    """
    
    D = []
    n = len(data)
    
    if boundary == 'pbc':
        for i in range(n):
            pi = data.iloc[i][['x','y','z']]
            L = [0 for i in range(i+1)]
            for j in range(i+1,n):
                pj = data.iloc[j][['x','y','z']]
                x_dist, y_dist, z_dist = abs(pi-pj)
                
                euclidean_distance = np.sqrt(min(box_size - x_dist,x_dist)**2
                                               + min(box_size - y_dist,y_dist)**2
                                               + min(box_size - z_dist,z_dist)**2)
                L += [euclidean_distance]
            D += [L]
    elif boundary == 'npbc':
        for i in range(n):
            a = data.iloc[i][['x','y','z']]
            L = [0 for i in range(i+1)]
            for j in range(i+1,n):
                b = data.iloc[j][['x','y','z']]
                x_dist = (a-b)['x']
                y_dist = (a-b)['y']
                z_dist = (a-b)['z']
                
                euclidean_distance = np.sqrt(x_dist**2
                                               + y_dist**2
                                               + z_dist**2)
                L += [euclidean_distance]
            D += [L]
    else:
        print('Boundaries of type {} are not understood, please use pbc for'.format(boundary) +
              ' periodic boundaries or npbc for non periodic ones.')
    M = np.array(D)
    #N = M.transpose() # Create the lower triangular matrix of distances.
    D = M #+N # Create the full matrix of distances.
    return D


def neighbours(data, searchRadius, box_size, boundary):
    """
    Return the data with a new column containing the points that are least in a
    sphere of radius searchRadius around each point of the data.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least three columns named x, y and z.
    searchRadius: float
        Radius of the sphere.

    Returns:
    --------
    data2: pandas dataframe
        A new dataframe with the neighbours added for each point.
    """
    
    M = dist_matrix(data, box_size, boundary)
    n = M.shape[0]
    D = []
    
    for i in range(n):
        #print("row {:02d} / {}".format(i+1,n))
        L = M[i].tolist()
        N = [[L[j],j,data.at[j,'polyIndex'],data.at[j,'beadPosition']] for j in range(n)]
        N = sorted(N)
        D_prime = []
        #print(N)
        j = i+1
        test = True
        while test:
            try:
                p = N[j][0]
                if p < searchRadius:
                    D_prime += [N[j]]
                    j += 1
                else:
                    test = False
            except IndexError:
                test = False      
        D += [D_prime]
    
    data2 = data.copy(deep=True)
    data2['nearest'] = D
    
    return data2


# =============================================================================
## Main function
# =============================================================================


def mainNeighbors(dir_in, dir_out, name, searchRadius = 1.5, box_size = 48, boundary = 'pbc'):
    """
    Do the neighbours search for a bunch of file and save the result in a (new)
    directory.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the directory named name.
    dir_out: string
        Output directory, that the script will eventually create and where it 
        will save the {name} subdirectory.
    name: string
        Subdirectory where the script will take the files.
    searchRadius: float, default 1.5
        Radius of the sphere.

    Returns:
    --------
    Nothing directly. The script will save the files in a new subdirectory 
    {name} under the output directory.
    """
    
    dir_path = Path(os.getcwd()) 
    path = dir_path.parent.joinpath(dir_in,name)
    save_path = dir_path.parent.joinpath(dir_out,name)
    save_path.mkdir(parents=True, exist_ok=True)
    
    # work on each files
    onlyfiles = [f for f in os.listdir(path) if isfile(join(path, f))]
    
    file_list = []
    for filename in onlyfiles:
        if filename.endswith('.pkl'):
            file_list.append(filename)
    file_list.sort()
    n = len(file_list)
    
    for i in range(n):
        print("file {:02d} / {}".format(i+1,n))
        file_path = os.path.join(path,file_list[i])
        df0 = pd.read_pickle(file_path)
        df1 = neighbours(df0, searchRadius, box_size, boundary)
        
        save_name = save_path.joinpath(file_list[i])
        df1.to_pickle(save_name)
        

L = ['6d48']
for elem in L:
    #mainNeighbors('neo_clean_pbc/4', 'neo_KNN/4', elem)
    mainNeighbors('neo_clean_pbc/65', 'neo_KNN/65', elem)