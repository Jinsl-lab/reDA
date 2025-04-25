import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cmocean
import phate
import scprep
import meld
import sklearn
import tempfile
import os
import scanpy as sc

# making sure plots & clusters are reproducible
np.random.seed(42)

#无需批次矫正
import anndata as ad
data = ad.read_h5ad('/.../ATACcnatest_final_seurat.h5ad')
df = pd.read_csv("/.../CRClabel.csv")
data.obs["condition"]=np.where(df["synth_labels"]=="Condition1",0,1)
data.obs["synth_samples"]=np.array(df["synth_samples"])

adata_pca=pd.DataFrame(data.obsm["X_lsi"].copy()) #load lsi matrix
adata_pca.index=data.obs.index
adata_pca=adata_pca.iloc[:,1:30]
metadata = data.obs

metadata['sample_labels']="a"
metadata.loc[metadata["synth_samples"]=="Condition2_R1",'sample_labels']="tyrA"
metadata.loc[metadata["synth_samples"]=="Condition2_R2",'sample_labels']="tyrB"
metadata.loc[metadata["synth_samples"]=="Condition2_R3",'sample_labels']="tyrC"
metadata.loc[metadata["synth_samples"]=="Condition2_R4",'sample_labels']="tyrD"
metadata.loc[metadata["synth_samples"]=="Condition2_R5",'sample_labels']="tyrE"
metadata.loc[metadata["synth_samples"]=="Condition2_R6",'sample_labels']="tyrF"
metadata.loc[metadata["synth_samples"]=="Condition2_R7",'sample_labels']="tyrG"
metadata.loc[metadata["synth_samples"]=="Condition2_R8",'sample_labels']="tyrH"
metadata.loc[metadata["synth_samples"]=="Condition1_R1",'sample_labels']="chdA"
metadata.loc[metadata["synth_samples"]=="Condition1_R2",'sample_labels']="chdB"
metadata.loc[metadata["synth_samples"]=="Condition1_R3",'sample_labels']="chdC"
metadata.loc[metadata["synth_samples"]=="Condition1_R4",'sample_labels']="chdD"
metadata.loc[metadata["synth_samples"]=="Condition1_R5",'sample_labels']="chdE"
metadata.loc[metadata["synth_samples"]=="Condition1_R6",'sample_labels']="chdF"
metadata.loc[metadata["synth_samples"]=="Condition1_R7",'sample_labels']="chdG"
metadata.loc[metadata["synth_samples"]=="Condition1_R8",'sample_labels']="chdH"

sample_cmap = {
    'chdA': '#fb6a4a',
    'chdB': '#000000',
    'chdC': '#FF0000',
    'chdD': '#00FF00',
    'chdE': '#0000FF',
    'chdF': '#FFFF00',
    'chdG': '#FF00FF',
    'chdH': '#C0C0C0',
    'tyrA': '#6baed6',
    'tyrB': '#00FFFF',
    'tyrC': '#800000',
    'tyrD': '#008000',
    'tyrE': '#000080',
    'tyrF': '#808000',
    'tyrG': '#800080',
    'tyrH': '#008080'
}

fig, ax = plt.subplots(1)

groups, counts = np.unique(metadata['sample_labels'], return_counts=True)
for i, c in enumerate(counts):
    ax.bar(i, c, color=sample_cmap[groups[i]])
    
ax.set_xticks(np.arange(i+1))
ax.set_xticklabels(groups)
ax.set_ylabel('# cells')

fig.tight_layout()

phate_op = phate.PHATE(n_jobs=-1)
adata_phate = phate_op.fit_transform(adata_pca)

scprep.plot.scatter2d(adata_phate, c=metadata['sample_labels'],cmap=sample_cmap, 
                      legend_anchor=(1,1), figsize=(6,5), s=10, label_prefix='PHATE', ticks=False)

metadata['genotype_name'] = np.where(metadata['sample_labels'].str.startswith("chd"), "chd", "tyr")
metadata['genotype'] = np.where(metadata['genotype_name'] == "chd", 1, 0)
metadata['replicate'] = metadata['sample_labels'].str[-1]

# beta is the amount of smoothing to do for density estimation
# knn is the number of neighbors used to set the kernel bandwidth
meld_op = meld.MELD(beta=67, knn=7)
sample_densities = meld_op.fit_transform(adata_pca, sample_labels=metadata['sample_labels'])

fig, axes = plt.subplots(2,8, figsize=(11,6))

for i, ax in enumerate(axes.flatten()):
    density = sample_densities.iloc[:,i]
    scprep.plot.scatter2d(adata_phate, c=density,
                          title=density.name,
                          vmin=0, 
                          ticks=False, ax=ax)
    
fig.tight_layout()

# This is a helper function to apply L1 normalization across the densities for each replicate
def replicate_normalize_densities(sample_densities, replicate):
    # Get the unique replicates
    replicates = np.unique(replicate)
    sample_likelihoods = sample_densities.copy()
    for rep in replicates:
        # Select the columns of `sample_densities` for that replicate
        curr_cols = sample_densities.columns[[col.endswith(rep) for col in sample_densities.columns]]
        curr_densities = sample_densities[curr_cols]
        # Apply L1 normalization
        sample_likelihoods[curr_cols] = sklearn.preprocessing.normalize(curr_densities, norm='l1')
    return sample_likelihoods

sample_likelihoods = replicate_normalize_densities(sample_densities, metadata['replicate'])

fig, axes = plt.subplots(1,2, figsize=(8.7,4))
experimental_samples = ['chdA','chdB','chdC','chdD','chdE','chdF','chdG','chdH']
scprep.plot.scatter2d(adata_phate, c=sample_likelihoods[experimental_samples].mean(axis=1), 
                      cmap=meld.get_meld_cmap(), vmin=0, vmax=1,
                      title='Mean', ticks=False, ax=axes[0])
scprep.plot.scatter2d(adata_phate, c=sample_likelihoods[experimental_samples].std(axis=1), vmin=0, 
                      cmap='inferno', title='St. Dev.', ticks=False, ax=axes[1])

fig.tight_layout()

#measure perturbation
metadata['chd_likelihood'] = sample_likelihoods[experimental_samples].mean(axis=1).values

