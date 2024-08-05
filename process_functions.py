import sys
sys.path.insert(1, '/cndd3/dburrows/CODE/admin_tools/')
from admin_tools import admin_functions as adm
import numpy as np
import pandas as pd

#==================================================
def rank_genes_old(data=None, names=None, gaba_mark=None, glu_mark=None):
#==================================================

    """
    Rank and identify genes related to GABAergic and glutamatergic markers for each cell.

    Parameters:
    -----------
    data : np.ndarray
        The gene expression data, typically `adata.X`, where each row corresponds to a cell and each column to a gene.
        
    names : list or np.ndarray
        The corresponding gene names for the columns in `data`. This should have the same length as the number of columns in `data`.
        
    gaba_mark : list
        A list of gene names that are markers for GABAergic neurons.
        
    glu_mark : list
        A list of gene names that are markers for glutamatergic neurons.

    Returns:
    --------
    data_l : list
        A list of tuples for each cell. Each tuple contains two arrays:
        - The first array represents the ranks (percentile) of the glutamatergic marker genes in descending order.
        - The second array represents the ranks (percentile) of the GABAergic marker genes in descending order.

    glu_l : list
        A list of arrays for each cell, containing the names of the glutamatergic marker genes found in descending order of expression.

    gaba_l : list
        A list of arrays for each cell, containing the names of the GABAergic marker genes found in descending order of expression.

    Notes:
    ------
    - The function iterates through each cell in the dataset.
    - For each cell, it sorts the genes by their expression levels.
    - It assigns numerical values to genes based on whether they are glutamatergic markers (1), GABAergic markers (2), or other genes (0).
    - It then calculates the percentile ranks for the marker genes and stores these values.
    - The function prints progress at regular intervals, approximately every 10% of the total number of cells.

    Examples:
    ---------
    >>> data = adata.X
    >>> names = adata.var_names
    >>> gaba_mark = ['GAD1', 'GAD2', 'SLC32A1']
    >>> glu_mark = ['SLC17A7', 'SLC17A6', 'CAMK2A']
    >>> data_l, glu_l, gaba_l = rank_genes(data, names, gaba_mark, glu_mark)
    """

    # Initialize lists to store the results
    data_l, glu_l, gaba_l = list(range(data.shape[0])), list(range(data.shape[0])), list(range(data.shape[0]))
    
    # Iterate over each cell
    for i in range(data.shape[0]):
        # Sort genes by expression levels
        data_ord, gene_ord = adm.sort_2list(data[i,:], names)
        
        # Assign numerical values based on gene markers
        ord_num = np.asarray([1 if x in glu_mark else 2 if x in gaba_mark else 0 for x in gene_ord])[::-1]
        ord_name = np.asarray([x if x in glu_mark else x if x in gaba_mark else 0 for x in gene_ord])[::-1]
        
        # Calculate percentile ranks for the marker genes
        data_l[i] = ((np.where(ord_num==1)[0]) / len(ord_num)) * 100, ((np.where(ord_num==2)[0]) / len(ord_num)) * 100
        
        # Extract the names of the marker genes
        glu_l[i] = ord_name[pd.Series(ord_name).isin(glu_mark)]
        gaba_l[i] = ord_name[pd.Series(ord_name).isin(gaba_mark)]
        
        # Print progress approximately every 10% of the total cells
        if i % (data.shape[0] // 10) == 0: 
            print((i / data.shape[0]) * 100)  

    glu_name = np.asarray(glu_l)
    gaba_name = np.asarray(gaba_l)
    glu_l = np.asarray([i[0] for i in data_l])
    gaba_l = np.asarray([i[1] for i in data_l])
    return(glu_l, gaba_l, glu_name, gaba_name)



#==================================================
def rank_genes(data=None, names=None, mark_l=None):
#==================================================

    """
    Rank and identify genes related to GABAergic and glutamatergic markers for each cell.

    Parameters:
    -----------
    data : np.ndarray
        The gene expression data, typically `adata.X`, where each row corresponds to a cell and each column to a gene.
        
    names : list or np.ndarray
        The corresponding gene names for the columns in `data`. This should have the same length as the number of columns in `data`.
        
    gaba_mark : list
        A list of gene names that are markers for GABAergic neurons.
        
    glu_mark : list
        A list of gene names that are markers for glutamatergic neurons.

    Returns:
    --------
    data_l : list
        A list of tuples for each cell. Each tuple contains two arrays:
        - The first array represents the ranks (percentile) of the glutamatergic marker genes in descending order.
        - The second array represents the ranks (percentile) of the GABAergic marker genes in descending order.

    glu_l : list
        A list of arrays for each cell, containing the names of the glutamatergic marker genes found in descending order of expression.

    gaba_l : list
        A list of arrays for each cell, containing the names of the GABAergic marker genes found in descending order of expression.

    Notes:
    ------
    - The function iterates through each cell in the dataset.
    - For each cell, it sorts the genes by their expression levels.
    - It assigns numerical values to genes based on whether they are glutamatergic markers (1), GABAergic markers (2), or other genes (0).
    - It then calculates the percentile ranks for the marker genes and stores these values.
    - The function prints progress at regular intervals, approximately every 10% of the total number of cells.

    Examples:
    ---------
    >>> data = adata.X
    >>> names = adata.var_names
    >>> gaba_mark = ['GAD1', 'GAD2', 'SLC32A1']
    >>> glu_mark = ['SLC17A7', 'SLC17A6', 'CAMK2A']
    >>> data_l, glu_l, gaba_l = rank_genes(data, names, gaba_mark, glu_mark)
    """

    # Initialize lists to store the results
    data_l, out_l = list(range(data.shape[0])), list(range(data.shape[0]))
    
    # Iterate over each cell
    for i in range(data.shape[0]):
        # Sort genes by expression levels
        data_ord, gene_ord = adm.sort_2list(data[i,:], names)
        
        # Assign numerical values based on gene markers
        ord_num = np.asarray([1 if x in mark_l else 0 for x in gene_ord])[::-1]
        ord_name = np.asarray([x if x in mark_l else 0 for x in gene_ord])[::-1]
        
        # Calculate percentile ranks for the marker genes
        data_l[i] = ((np.where(ord_num==1)[0]) / len(ord_num)) * 100
        
        # Extract the names of the marker genes
        out_l[i] = ord_name[pd.Series(ord_name).isin(mark_l)]
        
        # Print progress approximately every 10% of the total cells
        if i % (data.shape[0] // 10) == 0: 
            print((i / data.shape[0]) * 100)  
    
    out_name = np.asarray(out_l)
    out_ranks = np.asarray(data_l)
    return({'names':out_name, 'ranks':out_ranks})