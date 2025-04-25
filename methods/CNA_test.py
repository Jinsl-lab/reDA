import numpy as np
import scipy.sparse as sp
import pandas as pd
import anndata as ad
from multianndata import MultiAnnData
import cna
import scanpy as sc
from scipy.sparse.linalg import inv  
from scipy.sparse import csr_matrix,csc_matrix
from scipy.sparse.linalg import spsolve
import warnings
np.random.seed(0)
warnings.filterwarnings('ignore') #to ignore warnings on the output.

import anndata as ad
d1 = ad.read_h5ad('/.../ATACcnatest_final_seurat.h5ad')
d1.obs.index.name='cell'


df = pd.read_csv("/.../CRClabel.csv")
d1.obs["condition"]=np.where(df["synth_labels"]=="Condition1",0,1)
d1.obs["id"]=np.array(df["synid"])

d1.obs["id"]=d1.obs["id"].astype('int')
from multianndata import MultiAnnData
d1 = MultiAnnData(d1)
print('d1.samplem has', len(d1.samplem.columns), 'columns')
d1.obs_to_sample(['condition'])
print('now d1.samplem has', len(d1.samplem.columns), 'columns')


d1.obsp["connectivities"]=d1.obsp["distances"]
# perform association test for case/ctrl status, controlling for sex as a covariate and accounting for potential batch effect
res = cna.tl.association(d1,                   #dataset
            d1.samplem.condition)            #sample-level attribute of intest (case/control status)
 

print('\nglobal association p-value:', res.p)

