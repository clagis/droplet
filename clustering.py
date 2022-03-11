#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 13:48:26 2021

@author: clement
"""

# =============================================================================
## Library dependancies
# =============================================================================

import os
import pandas as pd
import numpy as np
from pathlib import Path
from os.path import isfile, join

# =============================================================================
## Basic functions
# =============================================================================


def clustering(data, sort = True):
    """
    Return a new dataframe containing the junctions, the coordinates of theirs
    centers of masses, theirs members and weight.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least three columns named x, y and z and one column
        named nearest, and two columns named polyIndex and beadPosition
    sort: bool, default True
        Indicate if you want the junction sorted or not.

    Returns:
    --------
    data3: pandas dataframe
        A new dataframe with the junctions, the coordinates of theirs
        centers of masses, theirs members and weight.
    """
    
    data2 = data.copy(deep=True)
    
    clusters_counter = -1
    L = [0 for i in range(len(data2))]
    full_members = {}
    for k in range(len(data2)):
        # if row has no cluster assigned, assign one and update the counter
        if L[k] == 0:
            L[k] = clusters_counter
            clusters_counter -= 1
        
        members = [data.iloc[k][['polyIndex','beadPosition']].to_list()]
        if data.at[k,'nearest'] != []:
            M = data.at[k,'nearest']
            for elem in M:
                neighbour = elem[1]
                members += [elem[2:]]
                if L[neighbour] == 0:
                    L[neighbour] = L[k]
                else:
                    I = [index for index, value in enumerate(L) if value == L[k]]
                    for idx in I:
                        L[idx] = L[neighbour]
        try:
            for member in members:
                if member not in full_members[L[k]]:
                    full_members[L[k]] += [member]
        except KeyError:
            full_members[L[k]] = members
    
    data2['junction'] = L
    
    data2.drop(['polyIndex','beadPosition','beadType','nearest'],axis=1,inplace=True)
    
    data3 = data2.groupby(['junction']).mean().reset_index()
    
    M = []
    for key in data3['junction'].to_list():
        M += [full_members[key]]
    
    data3['members'] = M
    data3['mass'] = data3['members'].apply(len)
    
    # reindexing junction
    index = 0
    for junction_number in data3['junction'].unique().tolist():
        data3.loc[data3['junction'] == junction_number,
                              ['junction']] = index
        index += 1
    if sort:
        data3.sort_values(by = 'junction', inplace=True, kind='quicksort')
        data3 = data3.astype({'junction': 'int64'})
        return data3
    else:
        data3 = data3.astype({'junction': 'int64'})
        return data3
    

def links(data_junctions, data_polymers, boundary ,box_size):
    """
    Return a new dataframe containing the junctions, the coordinates of theirs
    centers of masses, theirs members and weight.

    Parameters:
    ----------- 
    data_junctions: pandas dataframe
        data_junctions must contain at least one column named members.
    data_polymers: pandas dataframe
        data_polymers must contain at least three columns named x, y and z and 
        two columns named polyIndex and beadPosition.
    boundary: 'pbc' or 'npbc', default 'pbc'
        Type of boundary.
    box_size: int
        Size of the simulation box.     

    Returns:
    --------
    data_junctions: pandas dataframe
        The input data_junctions but with new columns giving the linked junctions,
        the minimal distance and the number of links between them.
    """
    
    #data_junctions['neighbours'] = [[] for i in range(len(data_junctions))]
    L = []
    length = data_polymers.groupby('polyIndex').count()['beadPosition'].head(1).item()
    if boundary == 'pbc':
        for i in range(len(data_junctions)):
            M = []
            for member in data_junctions.iloc[i]['members']:
                # L is a list with all the binding sites on the polymer save the one we look at.
                L = [(member[0],i) for i in range(length) if i != member[1]]
                N = []
                K = []
                for j in range(len(data_junctions)):
                    if j != i:
                        if not set(data_junctions.iloc[j]['members']).isdisjoint(L):
                            N += [j]
                            K += [[value for value in data_junctions.iloc[j]['members'] if value in L]]      

                for elem in N:
                    B = K[N.index(elem)]
                    D = []
                    for elem2 in B:
                        BS_number = elem2[1]
                        p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == member[1])][['x','y','z']]
                        p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == BS_number)][['x','y','z']]
 
                        x_dist, y_dist, z_dist = abs(p0-p1)

                        D += [np.sqrt(min(x_dist,box_size-x_dist)**2
                                      + min(y_dist,box_size-y_dist)**2
                                      + min(z_dist,box_size-z_dist)**2)]
                    
                    #data_junctions.at[i,'neighbours'].append((data_junctions.at[elem,'cluster'],min(D),len(B)))
                    M += [[data_junctions.at[elem,'cluster'],min(D),len(B)]]
            #data_junctions.at[i,'neighbours'].sort()
            M.sort()
            L += [M]
        data_junctions['neighbors'] = L
    elif boundary == 'npbc':
        pass
    else:
        print('Boundaries of type {} are not understood, please use pbc for'.format(boundary) +
              ' periodic boundaries or npbc for non periodic ones.')
        
    return data_junctions


def linksLinear(data_junctions, data_polymers, boundary, box_size):
    """
    Return a new dataframe containing the junctions, the coordinates of theirs
    centers of masses, theirs members and weight.

    Parameters:
    -----------
    data_junctions: pandas dataframe
        data_junctions must contain at least one column named members.
    data_polymers: pandas dataframe
        data_polymers must contain at least three columns named x, y and z and 
        two columns named polyIndex and beadPosition.
    boundary: 'pbc' or 'npbc'
        Type of boundary.
    box_size: int
        Size of the simulation box.     

    Returns:
    --------
    data_junctions: pandas dataframe
        The input data_junctions but with new columns giving the linked junctions,
        the minimal distance and the number of links between them.
    """
    
    data_junctions['neighbors'] = [[] for i in range(len(data_junctions))]
    length = data_polymers.groupby('polyIndex').count()['beadPosition'].head(1).item()
    if boundary == 'pbc':    
        for i in range(len(data_junctions)):
            N = []
            for member in data_junctions.iloc[i]['members']:
                K = [[member[0], i] for i in range(length) if i != member[1]]
                #print('K: {}'.format(K))
                if all(elem in data_junctions.iloc[i]['members'] for elem in K):
                    pass
                else:
                    if (member[1] == 0):
                        j = 0
                        while not ([member[0],1] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == 0)].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == 1)].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])

                            distance = np.sqrt(min(x_dist,box_size-x_dist)**2
                                          + min(y_dist,box_size-y_dist)**2
                                          + min(z_dist,box_size-z_dist)**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
                            #print('Neighbors of {}: {}'.format(i,N))
                        
                    elif (member[1] == (length-1)):
                        j = 0
                        while not ([member[0],length-2] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == (length-1))].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == (length-2))].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                            
                            distance = np.sqrt(min(x_dist,box_size-x_dist)**2
                                          + min(y_dist,box_size-y_dist)**2
                                          + min(z_dist,box_size-z_dist)**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
                    else:
                        b = member[1]
                        j = 0
                        while not ([member[0],b-1] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b)].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b-1)].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                            
                            distance = np.sqrt(min(x_dist,box_size-x_dist)**2
                                          + min(y_dist,box_size-y_dist)**2
                                          + min(z_dist,box_size-z_dist)**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
                            
                        j = 0
                        while not ([member[0],b+1] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b)].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b+1)].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                            
                            distance = np.sqrt(min(x_dist,box_size-x_dist)**2
                                          + min(y_dist,box_size-y_dist)**2
                                          + min(z_dist,box_size-z_dist)**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
            N.sort()  
            N1 = list(set([x[0] for x in N]))
            M = []
            for elem in N1:
                N2 = [elm[1] for elm in N if elm[0] == elem]
                dist = min(N2)
                strength = len(N2)
                M += [[elem, dist, strength]]
            data_junctions.at[i,'neighbors'] = M
    elif boundary == 'npbc':    
        for i in range(len(data_junctions)):
            N = []
            for member in data_junctions.iloc[i]['members']:
                K = [[member[0], i] for i in range(length) if i != member[1]]
                #print('K: {}'.format(K))
                if all(elem in data_junctions.iloc[i]['members'] for elem in K):
                    pass
                else:
                    if (member[1] == 0):
                        j = 0
                        while not ([member[0],1] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == 0)].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == 1)].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])

                            distance = np.sqrt(x_dist**2
                                          + y_dist**2
                                          + z_dist**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
                            #print('Neighbors of {}: {}'.format(i,N))
                        
                    elif (member[1] == (length-1)):
                        j = 0
                        while not ([member[0],length-2] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == (length-1))].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == (length-2))].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                            
                            distance = np.sqrt(x_dist**2
                                          + y_dist**2
                                          + z_dist**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
                    else:
                        b = member[1]
                        j = 0
                        while not ([member[0],b-1] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b)].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b-1)].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                            
                            distance = np.sqrt(x_dist**2
                                          + y_dist**2
                                          + z_dist**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
                            
                        j = 0
                        while not ([member[0],b+1] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b)].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b+1)].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                            
                            distance = np.sqrt(x_dist**2
                                          + y_dist**2
                                          + z_dist**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
            N.sort()  
            N1 = list(set([x[0] for x in N]))
            M = []
            for elem in N1:
                N2 = [elm[1] for elm in N if elm[0] == elem]
                dist = min(N2)
                strength = len(N2)
                M += [[elem, dist, strength]]
            data_junctions.at[i,'neighbors'] = M
    else:
        print('Boundaries of type {} are not understood, please use pbc for'.format(boundary) +
              ' periodic boundaries or npbc for non periodic ones.')
    return data_junctions


def linksFarLinear(data_junctions, data_polymers, boundary, box_size):
    """
    Return a new dataframe containing the junctions, the coordinates of theirs
    centers of masses, theirs members and weight.

    Parameters:
    -----------
    data_junctions: pandas dataframe
        data_junctions must contain at least one column named members.
    data_polymers: pandas dataframe
        data_polymers must contain at least three columns named x, y and z and 
        two columns named polyIndex and beadPosition.
    boundary: 'pbc' or 'npbc'
        Type of boundary.
    box_size: int
        Size of the simulation box.     

    Returns:
    --------
    data_junctions: pandas dataframe
        The input data_junctions but with new columns giving the linked junctions,
        the minimal distance and the number of links between them.
    """
    
    data_junctions['neighbors'] = [[] for i in range(len(data_junctions))]
    length = data_polymers.groupby('polyIndex').count()['beadPosition'].head(1).item()
    if boundary == 'pbc':    
        for i in range(len(data_junctions)):
            N = []
            for member in data_junctions.iloc[i]['members']:
                K = [[member[0], i] for i in range(length) if i != member[1]]
                #print('K: {}'.format(K))
                if all(elem in data_junctions.iloc[i]['members'] for elem in K):
                    pass
                else:
                    if (member[1] == 0):
                        j = 0
                        t = True
                        while t:
                            for k in range(1,length):
                                while not ([member[0],k] in data_junctions.iloc[j]['members']):
                                    j += 1
                                    if j >= len(data_junctions):
                                        break
                                if j < len(data_junctions):
                                    t = False
                                    p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == 0)].reset_index()[['x','y','z']]
                                    p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == 1)].reset_index()[['x','y','z']]
             
                                    x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
        
                                    distance = np.sqrt(min(x_dist,box_size-x_dist)**2
                                                  + min(y_dist,box_size-y_dist)**2
                                                  + min(z_dist,box_size-z_dist)**2)
                                    N += [[data_junctions.at[j,'junction'],distance]]
                                    #print('Neighbors of {}: {}'.format(i,N))
                                    break
                        
                    elif (member[1] == (length-1)):
                        j = 0
                        t = True
                        while t:
                            for k in range(length-2, 0, -1):
                                while not ([member[0],k] in data_junctions.iloc[j]['members']):
                                    j += 1
                                    if j >= len(data_junctions):
                                        break
                                if j < len(data_junctions):
                                    t = False
                                    p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == (length-1))].reset_index()[['x','y','z']]
                                    p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == (length-2))].reset_index()[['x','y','z']]
             
                                    x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                                    
                                    distance = np.sqrt(min(x_dist,box_size-x_dist)**2
                                                  + min(y_dist,box_size-y_dist)**2
                                                  + min(z_dist,box_size-z_dist)**2)
                                    N += [[data_junctions.at[j,'junction'],distance]]
                                    break
                    else:
                        b = member[1]
                        j = 0
                        t = True
                        while t:
                            for k in range(b-1, 0, -1):
                                while not ([member[0],k] in data_junctions.iloc[j]['members']):
                                    j += 1
                                    if j >= len(data_junctions):
                                        break
                                if j < len(data_junctions):
                                    t = False
                                    p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b)].reset_index()[['x','y','z']]
                                    p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b-1)].reset_index()[['x','y','z']]
             
                                    x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                                    
                                    distance = np.sqrt(min(x_dist,box_size-x_dist)**2
                                                  + min(y_dist,box_size-y_dist)**2
                                                  + min(z_dist,box_size-z_dist)**2)
                                    N += [[data_junctions.at[j,'junction'],distance]]
                                    break
                        j = 0
                        t = True
                        while t:
                            for k in range(b+1, length):
                                while not ([member[0],k] in data_junctions.iloc[j]['members']):
                                    j += 1
                                    if j >= len(data_junctions):
                                        break
                                if j < len(data_junctions):
                                    t = False
                                    p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b)].reset_index()[['x','y','z']]
                                    p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b+1)].reset_index()[['x','y','z']]
             
                                    x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                                    
                                    distance = np.sqrt(min(x_dist,box_size-x_dist)**2
                                                  + min(y_dist,box_size-y_dist)**2
                                                  + min(z_dist,box_size-z_dist)**2)
                                    N += [[data_junctions.at[j,'junction'],distance]]
                                    break
            N.sort()  
            N1 = list(set([x[0] for x in N]))
            M = []
            for elem in N1:
                N2 = [elm[1] for elm in N if elm[0] == elem]
                dist = min(N2)
                strength = len(N2)
                M += [[elem, dist, strength]]
            data_junctions.at[i,'neighbors'] = M
    elif boundary == 'npbc':    
        for i in range(len(data_junctions)):
            N = []
            for member in data_junctions.iloc[i]['members']:
                K = [[member[0], i] for i in range(length) if i != member[1]]
                #print('K: {}'.format(K))
                if all(elem in data_junctions.iloc[i]['members'] for elem in K):
                    pass
                else:
                    if (member[1] == 0):
                        j = 0
                        while not ([member[0],1] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == 0)].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == 1)].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])

                            distance = np.sqrt(x_dist**2
                                          + y_dist**2
                                          + z_dist**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
                            #print('Neighbors of {}: {}'.format(i,N))
                        
                    elif (member[1] == (length-1)):
                        j = 0
                        while not ([member[0],length-2] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == (length-1))].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == (length-2))].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                            
                            distance = np.sqrt(x_dist**2
                                          + y_dist**2
                                          + z_dist**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
                    else:
                        b = member[1]
                        j = 0
                        while not ([member[0],b-1] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b)].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b-1)].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                            
                            distance = np.sqrt(x_dist**2
                                          + y_dist**2
                                          + z_dist**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
                            
                        j = 0
                        while not ([member[0],b+1] in data_junctions.iloc[j]['members']):
                            j += 1
                            if j >= len(data_junctions):
                                break
                        if j < len(data_junctions):
                            p0 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b)].reset_index()[['x','y','z']]
                            p1 = data_polymers[(data_polymers['polyIndex'] == member[0]) & (data_polymers['beadPosition'] == b+1)].reset_index()[['x','y','z']]
     
                            x_dist, y_dist, z_dist = abs((p0-p1).iloc[0])
                            
                            distance = np.sqrt(x_dist**2
                                          + y_dist**2
                                          + z_dist**2)
                            N += [[data_junctions.at[j,'junction'],distance]]
            N.sort()  
            N1 = list(set([x[0] for x in N]))
            M = []
            for elem in N1:
                N2 = [elm[1] for elm in N if elm[0] == elem]
                dist = min(N2)
                strength = len(N2)
                M += [[elem, dist, strength]]
            data_junctions.at[i,'neighbors'] = M
    else:
        print('Boundaries of type {} are not understood, please use pbc for'.format(boundary) +
              ' periodic boundaries or npbc for non periodic ones.')
    return data_junctions

def neighboursList(data):
    """
    Return the same dataframe containing the list of neighbors for each junction.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least one column named neighbors.   

    Returns:
    --------
    data: pandas dataframe
        The input data but with a new column giving the list of neighbors
        for each junction.
    """
    
    data['neighbors_list'] = [[] for i in range(len(data))]
    
    for j in range(len(data)):
        L = []
        for elem in data.at[j,'neighbors']:
            L += [elem[0]]
        data.at[j,'neighbors_list'] = L
    
    return data


def connex(data,K,L):
    """
    Return the list of members of the connex component containing the elements
    originaly in L.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least one column named neighbors.
    K: list
        The list of members of the connex component so far.
    L: list
        The starting element, then the list of members of the connex component 
        reached so far minus the ones that are in K.

    Returns:
    --------
    K: list
        The junctions in the connex component containing the elements
        originaly in L.
    """
    
    if L == []:
        K.sort()
        return K
    
    M = []
    for elem in L:
        M += data.at[elem,'neighbors_list']
    
    N = []
    for elem2 in M:
        if not (elem2 in K):
            N += [elem2]
            
    N = list(set(N))
    L = N
    K += L
    return connex(data,K,L)


def connexComponents(data):
    """
    Return the same dataframe with a new column labelling the connex components.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least one column named neighbors_list.   

    Returns:
    --------
    data: pandas dataframe
        The input data but with a new column giving the connex components.
    """
    
    data['component'] = [-1 for i in range(len(data))]
    
    start = 0
    component = 0
    while (-1 in data['component'].to_list()):
        L = connex(data,[start],[start])
        for elem in L:
            data.at[elem,'component'] = component
        try:
            start = data['component'].to_list().index(-1)
        except ValueError:
            pass
        component += 1
        
    return data


def smallJunctions(data):
    """
    Return the same dataframe containing a boolean for each junction labelling 
    it as a small junction (mass < 3) or not.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least one column named mass.   

    Returns:
    --------
    data: pandas dataframe
        The input data but with a new column giving the list of neighbors
        for each junction.
    """
    
    data['small_junction'] = data['mass'].apply(lambda x: 1 if x < 3 else 0)
    return data


def bridges(data):
    """
    Return the same dataframe containing a boolean for each junction labelling 
    it as a bridge (mass = 2) connecting two not small junctions (mass > 2).

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least two columns named mass and neighbors_list.   

    Returns:
    --------
    data: pandas dataframe
        The input data but with a new column labelling each junction as a 
        bridge or not.
    """
    
    B = []
    for i in range(len(data)):
        if data.iloc[i]['mass'] == 2:
            test = [True,True]
            L = data.iloc[i]['neighbors_list']
            try:
                for i in range(2):
                    if data.iloc[L[i]]['mass'] < 2:
                        test[i] = False
                B += [test[0]*test[1]]
            except IndexError:
                B += [0]
        else:
            B += [0]
    data['bridge'] = B
    return data


def degree(data):
    """
    Return the same dataframe containing the degree and multidegree of each
    junction.

    Parameters:
    -----------
    data: pandas dataframe
        data must contain at least one column named neighbors.   

    Returns:
    --------
    data: pandas dataframe
        The input data but with new columns giving the degree and multidegree 
        of each junction.
    """
    
    degreeG = []
    degreeMG = []
    
    for i in range(len(data)):
        A = data.iloc[i][['small_junction','bridge']].to_list()
        if A == [1,0]:
            degreeG += [0]
            degreeMG += [0]
        elif A == [1,1]:
            degreeG += [2]
            degreeMG += [2]
        else:
            L = data.iloc[i]['neighbors_list']
            degG = 0
            degMG = 0
            for neighbor in L:
                if data.iloc[neighbor][['small_junction','bridge']].to_list() == [1,0]:
                    degG += 0
                    degMG += 0
                elif data.iloc[neighbor][['small_junction','bridge']].to_list() == [1,1]:
                    degG += 1
                    degMG += 1
                else:
                    degG += 1
                    M = data.iloc[i]['neighbors']
                    n = len(M)
                    j = 0
                    t = False
                    while ((j < n) and not t):
                        if M[j][0] == neighbor:
                            degMG += M[j][2]
                            t = True
                        else:
                            j += 1
            degreeG += [degG]
            degreeMG += [degMG]
    
    data['degreeG'] = degreeG
    data['degreeMG'] = degreeMG
    
    return data
        


# =============================================================================
## Main function
# =============================================================================


def mainClustering(dir_in, dir_out, name, box_size = 48, boundary = 'pbc'):
    """
    Do the clustering, linksLinear, neighbours_list, connex_components, 
    small_junctions, bridges and degree functions for a bunch of file and save
    the result in a (new) directory.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will search the directory named name.
    dir_out: string
        Output directory, that the script will eventually create and where it 
        will save the {name} subdirectory.
    name: string
        Subdirectory in the input directory where the files lie.

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
        print("{} file {:02d} / {}".format(name, i+1, n))
        file_path = os.path.join(path,file_list[i])
        df0 = pd.read_pickle(file_path)
        df1 = clustering(df0)
        linksLinear(df1, df0, boundary, box_size)
        neighboursList(df1)
        connexComponents(df1)
        smallJunctions(df1)
        bridges(df1)
        degree(df1)
        
        save_name = save_path.joinpath(file_list[i])
        df1.to_pickle(save_name)



L = ['6d48', '6e48', '6f48']
M = ['6e48', '6f48']

for elem in L:
    mainClustering('neo_KNN/65', 'neo_junctions/65', elem)
for elem in M:
    mainClustering('neo_KNN/4', 'neo_junctions/4', elem)



