#!/usr/bin/env python
# coding: utf-8

# In[8]:


from Bio.PDB import *
import numpy as np
from collections import OrderedDict


# In[9]:


# We will only consider the traditional 20 amino acids in this three letter list.
aa_3let_list = ['PHE',
 'GLU',
 'ASP',
 'LEU',
 'GLY',
 'VAL',
 'LYS',
 'THR',
 'ILE',
 'ARG',
 'ALA',
 'SER',
 'PRO',
 'HIS',
 'TYR',
 'MET',
 'ASN',
 'GLN',
 'CYS',
 'TRP']


# In[10]:


def identify_common_chains(identifier_list):
    """
    Description:
        This function will identify the common chains between several protein structures. This will 
        tell us which chains we want to perform a meaningful comparison upon 
    Input: 
        -identifier_list: A list of biopython structures for which we wish to compare
    Returns:
        A python list of the common chains that exist between all of the structures
    """
    common_chain_list = []
    for protein in identifier_list:
        current_chains = []
        for chain in protein[0]:
            current_chains.append(chain.id)
        if len(common_chain_list) == 0:
            common_chain_list = current_chains
        else:
            common_chain_list = list(set(common_chain_list) & set(current_chains))
    common_chain_list = sorted(common_chain_list)
    return common_chain_list


# In[11]:


def create_dihedral_array(common_chain_list, identifier_list, residue_range_dic):
    """
    Description:
        This function will create an array that represents all of the dihedral angles that exist
        within all of the protein structures in the identifier_list that are in the protein
        chains that are common to them, as listed in common_chain_list
    Input:
        -common_chain_list: Generated from identify_common_chains function, a list chains
        that are common among all of our structures
        -identifier_list: A list of biopython structures for which we wish to compare
        -residue_range_dic: An optional dictionary specifying which residues in what chains to retain for analysis
    Returns:
        -array_dic: A dictionary object with keys (label) that is a string of the protein.id and 
        the corresponding entry is an array of length *max_residues* by width 2 (phi and psi)
        -add_dic: A dictionary that helps keep track of the indexing since we unravel the chains
        into one long chain
    """
    if residue_range_dic != {}:
        for chain in residue_range_dic:
            start = residue_range_dic[chain][0]
            end = residue_range_dic[chain][1]
            residue_range_dic[chain] = range(start,end)
    # Set up the array size and fill it with nan values
    add_dic = OrderedDict()
    for chain in common_chain_list:
        add_dic[chain] = 0
    max_residues = 0
    for protein in identifier_list:
        num_residues = 0
        chain_max = 0
        for chain in common_chain_list:
            if add_dic[chain] < chain_max:
                add_dic[chain] = chain_max
            res_list = list(protein[0][chain])
            start = int(res_list[1].id[1])
            num_residues += (int(res_list[-1].id[1]))
            chain_max += (int(res_list[-1].id[1]))
        if num_residues > max_residues:
            max_residues = num_residues
    array_dic = {}
    for protein in identifier_list:
        label = str(protein.id)
        array_dic[label] = np.empty(shape=(max_residues+1,2), dtype=object)
        array_dic[label][:,:] = np.nan

    # Extract the phi-psi angles from all the proteins
    for protein in identifier_list:
        for chain in protein[0]:
            if chain.id in common_chain_list:
                index_add = add_dic[chain.id]
                res_list = list(protein[0][chain.id])
                start = int(res_list[1].id[1])
                end = int(res_list[-2].id[1])
                poly = Polypeptide.Polypeptide(chain)
                angle_array = poly.get_phi_psi_list()
                j = 0
                for res in (protein[0][chain.id]):
                    if res.resname in aa_3let_list:
                        res_num = int(res.id[1])+index_add
                        angles = angle_array[j]
                        if None not in angles:
                            if residue_range_dic == {}:
                                array_dic[protein.id][res_num,:] = angles
                            # Added code for specifying specific chains/resnums
                            else:
                                raw_resnum = int(res.id[1])
                                if raw_resnum in residue_range_dic[chain.id]:
                                    array_dic[protein.id][res_num,:] = angles
                                else:
                                    array_dic[protein.id][res_num,:] = np.nan
                        j += 1
    return array_dic, add_dic

