#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 18:42:01 2020

@author: clement
"""

#=============================================================================
# Libraries
#=============================================================================

import pandas as pd
import os
import shutil
import re
from pathlib import Path
from os.path import isfile, join

#=============================================================================
# Cleaning
#=============================================================================

# Import data and convert them to dataframe
def clean(filename, boundary):
    """
    Clean the files, only keeping the relevant data.

    Parameters:
    -----------
    filename: string
        Name of the file that we want to use as input. It must be a .rst file
        with 14 columns.
    boundary: string, accept only 'pbc' or 'npbc'
        Type of boundary used. Determine which coordinates we will keep.

    Returns:
    --------
    df_pbc or df_npbc: pandas dataframe
        The cleaned dataframe.
    """
    
    col_names = ['polyIndex',
           'polyType',
           'beadIndex',
           'beadType',
           'a',
           'x', 'y', 'z',
           'x_npbc', 'y_npbc', 'z_npbc',
           'b', 'c', 'd'] 
    df = pd.read_csv(filename,
                 names=col_names,
                 header=None,
                 sep = ' ',
                 index_col=False)
    
    if boundary == 'pbc':
        # Clean data (delete columns 2, 5, 12-14) by removing non used columns
        df_pbc = df.drop(['polyType', 'a', 'b', 'c', 'd', 'x_npbc', 'y_npbc', 'z_npbc'], axis=1)
        header_list_pbc = ['polyIndex', 'beadIndex', 'beadType', 'beadPosition', 'x', 'y',
                'z']
        df_pbc = df_pbc.reindex(columns = header_list_pbc)
        # Drop useless labels
        df_pbc.drop(['beadIndex'], axis=1, inplace=True)
        
        return df_pbc
    
    elif boundary == 'npbc':
        # Clean data (delete columns 2, 5, 12-14) by removing non used columns
        df_npbc = df.drop(['polyType', 'a', 'b', 'c', 'd', 'x', 'y', 'z'], axis=1)
       
        header_list_npbc = ['polyIndex', 'beadIndex', 'beadType', 'beadPosition', 'x_npbc', 'y_npbc',
                'z_npbc']
        df_npbc = df_npbc.reindex(columns = header_list_npbc)
        df_npbc.rename({'x_npbc': 'x', 'y_npbc': 'y', 'z_npbc': 'z'},inplace=True)
        # Drop useless labels
        df_npbc.drop(['beadIndex'], axis=1, inplace=True)
             
        return df_npbc
    
    else:
        print('Boundaries of type {} are not understood, please use pbc for'.format(boundary) +
              ' periodic boundaries or npbc for non periodic ones.')

#=============================================================================
# Formating
#=============================================================================

# Regroup by polyIndex & beadPosition to calculte coordinates of binding sites
def BS_tag(data):
    """
    Remove the non-bidding beads and fuse the bidding ones in the same site.

    Parameters:
    -----------
    data: pandas dataframe
        must contain at least three columns named polyIndex, beadPosition and 
        beadType.

    Returns:
    --------
    data3: pandas dataframe
        The cdataframe with only the bidding sites for each polymer.
    """
    
    data2 = data.copy(deep=True)
    L = []
    
    """
    For each polyIndex we create a small dataframe and label the beads if they
    are not in the backbone by the index of the first one in the chain and by
    -1 otherwise. So the beads in the first bidding site will all be labeled 0, etc.
    """
    for elem in data2.groupby(['polyIndex']).mean().index.to_list():
        beadP = 0
        working_df = data2.loc[data2['polyIndex'] == elem, 'beadType']
        idx = working_df.index.to_list()
        for k in idx:
            if working_df.at[k] != 2:
                L += [beadP]
            else:
                L += [-1]
                beadP += 1
    
    # Drop the irrelevant beads (backbone)
    data2['beadPosition'] = L
    index_drops = data2[data2['beadType'] == 2].index
    data2.drop(index_drops,inplace=True)
    
    # Aggregate each bidding site
    data3 = data2.groupby(['polyIndex','beadPosition']).mean().reset_index()
    
    # Relabel by their position on the polymer (needed for linear links)
    positions = list(data3['beadPosition'].unique())
    data3['beadPosition'] = data3['beadPosition'].apply(lambda x: positions.index(x))
    
    data3 = data3.astype({'beadPosition': 'int64'})
    return data3


#=============================================================================
# Main function to use
#=============================================================================

def main_clean(dir_in, dir_out, boundary = 'pbc'):
    """
    Do the neighbours search for a bunch of file and save the result in a (new)
    directory.

    Parameters:
    -----------
    dir_in: string
        Input directory, where the script will take all the files.
    dir_out: string
        Output directory, that the script will eventually create and where it 
        will save the files.
    boundary: string, accept only 'pbc' or 'npbc', default 'pbc'
        Type of boundary used. Determine which coordinates we will keep.

    Returns:
    --------
    Nothing directly. The script will save the files in the output directory.
    """
    
    dir_path = Path(os.getcwd())
    path = dir_path.parent.joinpath(dir_in)

    # make the outputs directories
    clean_path = dir_path.parent.joinpath(dir_out)
    clean_path.mkdir(parents=True, exist_ok=True)

    # start working on the raw data
    onlyfiles = [f for f in os.listdir(path) if isfile(join(path, f))]
    onlyfiles = sorted(onlyfiles)
    
    # starting parameters for the counter
    k = 1
    n = len(onlyfiles)
    for file in onlyfiles:
        if file.endswith('.rst'):
            # counter to see progression
            print("file {:02d} / {}".format(k,n))
            k += 1
            file_name = file[:-4]
            # create clean dataframe from raw data
            file_path = os.path.join(path,file)
            df = BS_tag(clean(file_path,boundary))
            
            # save cleanded data file in clean_data subfolder
            if re.match('[\w.]*con.\d{3,5}$',file_name):
                if re.match('[\w.]*con.\d\d\d$',file_name):
                    file_name = file_name[:-3] + '000' + file_name[-3:]
                elif re.match('[\w.]*con.\d\d\d\d$',file_name):
                    file_name = file_name[:-4] + '00' + file_name[-4:]
                else:
                    file_name = file_name[:-5] + '0' + file_name[-5:]
            savedfile_name = file_name + '.pkl'
            df.to_pickle(savedfile_name)
            
            savedfile_path = clean_path.joinpath(savedfile_name)
            if os.path.exists(savedfile_path):
                os.remove(savedfile_path)
                shutil.move(savedfile_name, clean_path)
            else:
                shutil.move(savedfile_name, clean_path)
       
        else:
            pass
        
        

#=============================================================================
# Tests
#=============================================================================

# dir_path = Path(os.getcwd())
# path = dir_path.parent.joinpath('rst3')
# file_path = os.path.join(path,'dmpccs.rod4030r5.con.50000.rst')
# df = final_clean(file_path,'pbc')
