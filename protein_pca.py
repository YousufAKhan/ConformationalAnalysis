#!/usr/bin/env python
# coding: utf-8

# In[13]:


import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
import pandas as pd
import numpy as np
# To help with 3D plotting
# from adjustText import adjust_text
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.lines import Line2D


# In[2]:


def protein_pca(comps, eval_comp, scaled_data, label, identifier_list):
    """
    Description:
        Take an input dataframes of features for protein sequences and perform PCA
    Input:
        -comps: The number of components to perform for PCA
        -eval_comp: The number of components to consider of the principal components
        -scaled_data: Dataframe of features/values for your proteins of interest
        -label: The label for the scaled_data dataframe
        -identifier_list: The initial list of protein structure objects
    Returns:
        -var_contribution: The variance contributed by each feature for the selected PCA
        -pca: PCA results
    """
    if len(identifier_list) < comps:
        n_comp = len(identifier_list)
        eval_comp = len(identifier_list)
    else:
        n_comp = comps
    pca = PCA(n_components=n_comp)
    pc = pca.fit_transform(scaled_data)
    column_list = []
    for i in range(n_comp):
        name = 'pc'+str(i+1)
        column_list.append(name)
    final_df = pd.DataFrame(data = pc, columns = column_list)
    index_list = list(label['Class'])
    plot_pca(final_df, n_comp, index_list, pca)
    var_contribution = (abs(pca.components_))
    return var_contribution, pca


# In[11]:


def get_cmap(n):
    col_list = []
    for i in range(n):
        col = (np.random.random(), np.random.random(), np.random.random())
        col_list.append(col)
    return col_list


# In[4]:


def plot_pca(final_df, n_comp, index_list, pca):
    """
    Description:
        Display all the plots related to the PCA run.
    Input:
        -final_df: The dataframe with the output of the PCA results
        -n_comp: The number
        -index_list: label of PC components
        -pca: The data from the pca
    """
    # Generate a 3D PCA plot with the first three principal components
    colors = get_cmap(len(index_list))
    color_dic = {}
    color_list = []
    for i,item in enumerate(index_list):
        prefix = item.split('_')[0]
        if (len(color_dic) == 0) or (prefix not in color_dic.keys()):
            color_dic[prefix] = colors[i]
            color_list.append(color_dic[prefix])
        else:
            color_list.append(color_dic[prefix])
    if n_comp > 2:
        fig = plt.figure(figsize = (12,12))
        ax_p = fig.add_subplot(1,1,1, projection='3d') 
        ax_p.set_xlabel('pc1', fontsize = 15)
        ax_p.set_ylabel('pc2', fontsize = 15)
        ax_p.set_zlabel('pc3', fontsize = 15)
        x = final_df['pc1']
        y = final_df['pc2']
        z = final_df['pc3']
        texts = []
        for i in range(len(x)): # plot each point + it's index as text above
            ax_p.scatter(x[i],y[i],z[i],c=color_list[i]) 
            # texts.append(ax_p.text(x[i],y[i],z[i], '%s' % (str(index_list[i])), size=12, zorder=1,  color='k'))
        #adjust_text(texts)
        ax_p.grid()
        legend_list = []
        legend_colors = []
        for entry in color_dic:
            legend_list.append(entry)
            legend_colors.append(color_dic[entry])
        custom_lines = []
        for i in legend_colors:
            custom_lines.append(Line2D([0], [0], color=i, lw=1))
        ax_p.legend(custom_lines, legend_list)
    plt.savefig('PCA_figure.pdf')


    # 2D PCA plots
    dim_a = int(np.ceil(np.sqrt(n_comp-1)))
    figs, axs = plt.subplots(dim_a, dim_a, figsize=(12, 12))
    columns = list(final_df.columns)
    z = 0
    break_out_flag = False
    for i in range(len(columns)):
        for j in range(dim_a):
            x = final_df[columns[z]]
            y = final_df[columns[z+1]]
            axs[i, j].scatter(x, y, c=color_list)
            axs[i, j].set_xlabel(columns[z])
            axs[i, j].set_ylabel(columns[z+1])
            z += 1
            if (z+1 == len(columns)):
                break_out_flag = True
                break
        if break_out_flag:
            break
    figs.tight_layout()
    #Plot the principal components
    x = []
    for i in range(n_comp):
        x.append(i+1)
    y = pca.explained_variance_
    sum_y = np.sum(y)
    y = 100*y/sum_y
    f2 = plt.figure(figsize = (4,4))
    ax2 = f2.add_subplot(111)
    ax2.scatter(x,y)
    plt.plot(x,y)
    plt.title('Variance by component plot')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance (%)')
    plt.savefig('Variance_contribution.pdf')


# In[5]:


def variance_contribution(var_contribution, eval_comp, data):
    """
    Description: 
        Calculate the variance contributed by each component up to the number of components we wnated to evaluate
    Inputs:
        -var_contribution: The variance contributed by each feature for the selected PCA.(n_comp x features dataframe)
        -eval_comp: The number of components we want to include
        -data: dataframe with all features
    Returns:
        -var_dic: Dictionary of residues and the variance they contribute
    """
    if eval_comp > len(var_contribution):
        eval_comp = len(var_contribution)
    total_pc = var_contribution[0]
    for i in range(eval_comp):
        if i != 0:
            total_pc = total_pc + var_contribution[i]
    var_dic = {}
    for i,column in enumerate(data):
        residue_val = column.split(':')[1]
        residue_val = residue_val.split('_')[0]
        if residue_val not in var_dic:
            var_dic[residue_val] = total_pc[i]
        else:
            var_dic[residue_val] += total_pc[i]
    return var_dic

