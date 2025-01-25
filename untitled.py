import pandas as pd
import scanpy as sc
import numpy as np
import sys
import matplotlib.pyplot as plt
from tqdm import tqdm
import seaborn as sns
import process_functions as pf
sys.path.insert(1, '/cndd3/dburrows/CODE/admin_tools/')
from admin_tools import admin_functions as adm
%load_ext autoreload

#PATHS
data_path = '/datasets/Public_Datasets/KIBM_Spitzer_NeurotransmitterSwitching/'
out_path = '/cndd3/dburrows/DATA/neurotransmitter_switch/processed/3batch/'
eran_path = '/cndd/emukamel/KIBM_neuroswitch/Analyze_3batches/'


# Finalise gene panel genes
df = pd.read_excel('/cndd3/dburrows/DATA/neurotransmitter_switch/merfish/gene_panel/BICCN_merfish_panel.xlsx')
new_header = df.iloc[0]
df = df[1:]
df.columns = new_header
biccn_gene = df['Gene name'].unique()
biccn_gene = np.sort(np.asarray([g for g in biccn_gene if isinstance(g, str)]))

#PLC genes + gradient genes
plc_genes = np.asarray(['Sema3e', 'Hpcal1', 'Mapk3', 'Junb', 'Nr4a1', 'Npas4', 'S100a6', 'Grp', 'Ptn', 'Nr2f1', 'Lhx2'])

#TFs
tfs = {
    "Slc17a7": [
        {"geneID": "ENSMUSG00000043206", "name": "Tbr1"},    # TBR1
        {"geneID": "ENSMUSG00000037752", "name": "Neurod2"}, # NEUROD2
        {"geneID": "ENSMUSG00000029876", "name": "Cux1"},    # CUX1
        {"geneID": "ENSMUSG00000052172", "name": "Cux2"},    # CUX2
        {"geneID": "ENSMUSG00000038505", "name": "Satb2"},   # SATB2
        {"geneID": "ENSMUSG00000005886", "name": "Bcl11b"},  # BCL11B (CTIP2) - deep-layer excitatory
        {"geneID": "ENSMUSG00000030252", "name": "Rorb"},    # RORB - nuclear receptor TF
        {"geneID": "ENSMUSG00000070377", "name": "Fezf2"},   # FEZF2 (Fezl) - L5 subcortical projection
        {"geneID": "ENSMUSG00000029099", "name": "Tle4"},    # TLE4 (E(spI)) - deep-layer co-repressor
        {"geneID": "ENSMUSG00000024221", "name": "Foxp1"},   # FOXP1 - associated w/ L5 excitatory identity
        {"geneID": "ENSMUSG00000019920", "name": "Satb1"},   # SATB1 - can co-express w/ BCL11B in deep layer
    ],
    "Gad1/Gad2": [
        {"geneID": "ENSMUSG00000058942", "name": "Lhx6"},    # LHX6 - MGE-derived interneurons
        {"geneID": "ENSMUSG00000038246", "name": "Sox6"},    # SOX6 - interneuron subtype maintenance
        {"geneID": "ENSMUSG00000067233", "name": "Npas4"},   # NPAS4 - activity-dependent regulator
        {"geneID": "ENSMUSG00000038531", "name": "Dlx1"},    # DLX1 - GABAergic phenotype
        {"geneID": "ENSMUSG00000018285", "name": "Dlx2"},    # DLX2 - GABAergic phenotype
        {"geneID": "ENSMUSG00000020710", "name": "Dlx5"},    # DLX5 - continuing GABAergic specification
        {"geneID": "ENSMUSG00000038124", "name": "Dlx6"},    # DLX6 - same cluster with Dlx1/2/5
        {"geneID": "ENSMUSG00000005309", "name": "Mef2c"},   # MEF2C - synaptic plasticity in interneurons
        {"geneID": "ENSMUSG00000029150", "name": "Arx"},     # ARX - MGE-derived GABAergic identity
        {"geneID": "ENSMUSG00000020739", "name": "Maf"},     # MAF - can be expressed in subsets of inhibitory cells
        {"geneID": "ENSMUSG00000031318", "name": "Mafb"}     # MAFB - related to MAF, in some interneuron populations
    ]
}

#DEGs
total = np.unique(np.concatenate([biccn_gene, np.asarray([tfs['Slc17a7'][i]['name'] for i in range(len(tfs['Slc17a7']))])
                        , np.asarray([tfs['Gad1/Gad2'][i]['name'] for i in range(len(tfs['Gad1/Gad2']))]),
                        coex_genes, plc_genes])


#cluster
#load in our filtered rna data 
fdata = sc.read_h5ad(f'{out_path}KIBM_adata_downstream_dissection-definiteremoved.h5ad')
subdata = fdata.copy()
subdata.var['gene'] = subdata.var['gene'].astype(str)
subdata.var = subdata.var.set_index('gene')

#get all genes
_data = sc.read_h5ad(out_path + 'KIBM_adata_mapmycells.h5ad')
alldata = _data[subdata.obs.index]
columns_to_add = ['cluster_label', 'condition',  'is_neuron', 'neuron_type', 'leiden',
                  'dissection_class', 'uncertain_class','class_label']
alldata.obs = alldata.obs.join(subdata.obs[columns_to_add])

alldata.raw = alldata

#Normalise
##========================
sc.pp.normalize_total(alldata, target_sum=1e6)
sc.pp.log1p(alldata, base=10)

alldata.var['gene'] = alldata.var['gene'].astype(str)
alldata.var = alldata.var.set_index('gene')
alldata.var_names_make_unique()

#subset with merfish genes
int_gene=np.intersect1d(total, alldata.var_names)
newdata = alldata[:,int_gene].copy()