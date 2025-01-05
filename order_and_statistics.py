#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
from collections import OrderedDict


# In[2]:


def min_max(var_dic):
    """
    Description:
        Take an input dictionary of values and return a sorted, 0-1 normalized dictionary
    Input:
        -var_dic: A dictionary of variability for each feature
    Returns:
        -sort_orders: A 0-1 normalized, sorted dictionary of variances
    """
    min_v = np.min(list(var_dic.values()))
    max_v = np.max(list(var_dic.values()))
    norm_sig = {}
    factor = 1/(max_v-min_v)
    for i in var_dic:
        adj = var_dic[i]*factor - min_v*factor
        norm_sig[i] = adj
    sort_orders = sorted(norm_sig.items(), key=lambda x: x[1], reverse=True)
    return sort_orders


# In[5]:


def get_top_n(sort_orders, common_chain_list, add_dic, n = 5):
    """
    Description:
        Take an input dictionary of sorted, normalized variances and return the top n%
        features that contribute to the features with the final chain and res-number
        that maps 1:1 to the original structure's labeling scheme
    Input:
        -sort_orders: A 0-1 normalized, sorted dictionary of variances
        -common_chain_list: A python list of the common chains that exist between all of the structures
        -add_dic: A dictionary that helps keep track of the indexing since we unravel the chains
        into one long chain
    Returns:
        -final_res_dict: A final dictionary of chain, residue_numbers that we consider to be of interest
    """
    # top n% of residues by variance
    factor = n/100
    num_top = int(len(sort_orders)*factor)
    top_orders = OrderedDict()
    top_resis = sort_orders[0:num_top]
    for i in (top_resis):
        top_orders[i[0]] = i[1]
    # The function below only takes the top 5%
    res_keep = []
    res_dict = {}
    res_dict = top_orders
    final_res_dict = {}
    for chain in common_chain_list:
        final_res_dict[chain] = {}
    for res in res_dict:
        res = int(res)
        prev_chain = ''
        for chain in add_dic:
            if res > add_dic[chain]:
                prev_chain = chain
            else:
                break
        final_chain = prev_chain
        final_resnum = res-add_dic[final_chain]
        final_res_dict[final_chain][final_resnum] = res_dict[str(res)]
    return final_res_dict


# In[ ]:




