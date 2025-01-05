#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np


# In[2]:


def create_dihedral_dataframe(array_dic):
    """
    Description:
        Take the array dic and convert it into an easily accessible pandas dataframe
    Input:
        -array_dic: Dictionary of all possible residue positions and their respective phi-psi angles
    Returns:
        -all_classes_df: A dataframe of all the phi-psi angles for all the structures
    """
    all_classes_df = pd.DataFrame()
    for protein in array_dic:
        prot_df = pd.DataFrame(array_dic[protein])
        prot_phi = list(prot_df[0])
        prot_psi = list(prot_df[1])
        prot_phi_index = []
        for i in range(len(prot_phi)):
            prot_phi_index.append('Phi:'+str(i))
        prot_psi_index = []
        for i in range(len(prot_psi)):
            prot_psi_index.append('Psi:'+str(i))
        prot_phi_df = pd.DataFrame(pd.DataFrame(prot_phi).T)
        prot_phi_df.columns = prot_phi_index
        prot_phi_df.reset_index
        prot_psi_df = pd.DataFrame(pd.DataFrame(prot_psi).T)
        prot_psi_df.columns = prot_psi_index
        prot_label_df = pd.DataFrame([protein])
        prot_label_df.columns = ['Class']
        prot_phi_psi = pd.concat([prot_label_df,prot_phi_df,prot_psi_df], axis=1)
        if len(all_classes_df) == 0:
            all_classes_df = prot_phi_psi
        else:
            all_classes_df = pd.concat([all_classes_df,prot_phi_psi], axis=0)
    all_classes_df = all_classes_df.dropna(axis=1)
    return all_classes_df


# In[3]:


def parametrize_angles(all_classes_df):
    """
    Description:
        This function takes the angles and encodes them into two numerical x and y values using a
        sine and cosine transformation
    Input:
        -all_classes_df: A dataframe with rows representing different structures, and the columns
        the phi-psi angles
    Returns:
        -parametrized_df: A dataframe with now double the amount of columns due to more features
    """
    input_dic = {'Class':list(all_classes_df['Class'])}
    parametrized_df = pd.DataFrame(input_dic)
    for column in all_classes_df:
        if column != 'Class':
            x = []
            y = []
            for value in all_classes_df[column]:
                x.append(np.cos(value))
                y.append(np.sin(value))
            name_x = column+'_x'
            name_y = column+'_y'
            parametrized_df[name_x] = x
            parametrized_df[name_y] = y
    return parametrized_df


# In[ ]:




