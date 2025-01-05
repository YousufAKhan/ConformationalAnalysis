#!/usr/bin/env python
# coding: utf-8

# This is a script that will serve as a generic input file that will run the conformation analysis

# In[ ]:


# Basic imports
import os
import argparse
import ast
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pylab
from collections import OrderedDict
import matplotlib
# Sklearn
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
# Import biopython and setup parser
from Bio.PDB import *
parser = PDBParser()
# Conformation tool specific files
import structure_processing
import angle_processing
import order_and_statistics
import protein_pca


# In[ ]:


def process_files(directory, comps=5, eval_comp=1, common_chain_list=[], residue_range_dic={}):
    # Read in the files from the directory
    print("Processing input files")
    parser = PDBParser()
    identifier_list = []
    input_files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
    for file in input_files:
        name = file.split('_')[0]
        run = file.split('_real')[0]
        run = run.split('run')[1]
        fn = name+run
        identifier_list.append(parser.get_structure(fn, os.path.join(directory,file)))

    # Identify the common chains between the structures
    print("Identifying common chains between structures if not specified earlier")
    if common_chain_list == []:
        common_chain_list = structure_processing.identify_common_chains(identifier_list)
    # Create an array of dihedrals
    print("Extracting dihedral information")
    array_dic, add_dic = structure_processing.create_dihedral_array(common_chain_list, identifier_list, residue_range_dic)
    # Set up an easily accessible pandas dataframe
    all_classes_df = angle_processing.create_dihedral_dataframe(array_dic)
    # Parametrize the angles into x and y values 
    parametrized_df = angle_processing.parametrize_angles(all_classes_df)
    data = parametrized_df.drop('Class',axis=1)
    label = parametrized_df[['Class']]
    # Scaling the phi and psi features prior to PCA
    # Currently scaling is turned on
    print("Normalizing data")
    scaled_data = StandardScaler().fit_transform(data)
    # scaled_data = data
    # Perform PCA
    print("Performing PCA")
    var_contribution, pca = protein_pca.protein_pca(comps,eval_comp,scaled_data,label,identifier_list)
    var_dic = protein_pca.variance_contribution(var_contribution, eval_comp, data)
    # Normalize values to 0-1
    print("Final statistics and outputting files")
    sort_orders = order_and_statistics.min_max(var_dic)
    final_res_dict = order_and_statistics.get_top_n(sort_orders, common_chain_list, add_dic)
    # Save the values as b factors
    for protein in identifier_list:
        for chain in protein[0]:
            if chain.id not in add_dic:
                for res in protein[0][chain.id]:
                    for atom in res:
                        atom.set_bfactor(0)
            else:
                for res in protein[0][chain.id]:
                    res_num = int(res.id[1])
                    if res_num not in final_res_dict[chain.id]:
                        for atom in res:
                            atom.set_bfactor(0)
                    else:
                        for atom in res:
                            atom.set_bfactor(final_res_dict[chain.id][res_num])
    for protein in identifier_list:
        io = PDBIO()
        io.set_structure(protein)
        filename = str(protein.id)+'_all_bfactor_colored_REV2.pdb'
        io.save(filename)


# In[ ]:


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process files in a directory with optional parameters.")
    parser.add_argument("directory", help="The directory path where the files are located")
    parser.add_argument("--comps", type=int, default=5, help="The number of components for analysis (default is 5)")
    parser.add_argument("--eval_comp", type=int, default=1, help="The number of components for evaluating differences (default is 1)")
    parser.add_argument("--common_chains", nargs="*", default=[], help="Chains to explicitly compare (space-separated list, default is none)")
    parser.add_argument("--residue_range_dic", type=str, default="{}", help="Residue range specifier (default is empty dictionary)")
    # Example of what this might look like
    # python script_name.py /path/to/directory --comps 10 --eval_comp 3 --common_chains A B C --residue_range_dic "{'A':[205,583], 'B':[205,583], 'C':[205,583], 'D':[205,583], 'E':[205,583], 'F':[205,583]}"
    args = parser.parse_args()
    directory_path = args.directory
    comps = args.comps
    eval_comp = args.eval_comp
    common_chain_list = args.common_chains
    # Safely evaluate the residue_range_dic from string to dictionary
    try:
        residue_range_dic = ast.literal_eval(args.residue_range_dic)
    except ValueError:
        print("Error: Invalid residue_range_dic format. Please provide a valid dictionary-like input.")
        exit(1)

    # Check if the provided path exists and is a directory
    if os.path.exists(directory_path) and os.path.isdir(directory_path):
        # Call the function to process files in the given directory
        process_files(directory_path, comps=comps, eval_comp=eval_comp, common_chain_list=common_chain_list, residue_range_dic=residue_range_dic)
    else:
        print("Invalid directory path. Please provide a valid path to the directory.")

