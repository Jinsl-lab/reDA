#kernal python27
from pencil import *
import scanpy as sc
from matplotlib.colors import LinearSegmentedColormap

import os
os.environ['CUDA_VISIBLE_DEVICES'] = '0' #select a gpu id, otherwise set to '-1'.
#import torch
#torch.cuda.set_device(0)

import warnings
warnings.filterwarnings('ignore') #to ignore warnings on the output. 

adata = sc.read_h5ad('/.../ATACcnatest_final_seurat.h5ad')
df = pd.read_csv("/.../CRClabel.csv")
adata.obs["condition"]=np.where(df["synth_labels"]=="Condition1",0,1)

metadata = adata.obs
metadata['sample_labels']=np.where(metadata['condition']==1, "class_1", "class_2")

labels_raw = metadata.sample_labels
labels_raw = pd.Categorical(labels_raw)


from scipy.io import mmread
data = adata.obsm['X_harmony'].copy()  
data = data[:,1:30]
labels = labels_raw.codes
class_names = list(labels_raw.categories)
emd = adata.obsm['X_umap']

data_name = 'PENCIL_CRC_FINAL'
expr_id = '0.0.1'
mode = 'multi-classification'  
pencil = Pencil(mode, select_genes=True, seed=1234, data_name=data_name, expr_id=expr_id, mlflow_record=True)
with mlflow.start_run():
    pred, confidence = pencil.fit_transform(
        data, labels, 
        test=True, 
        shuffle_rate=1/4,
        lambda_L1=1e-5, 
        lambda_L2=1e-3, 
        lr=0.01,  
        class_weights=None,
        class_names=class_names, 
        emd=emd,
        plot_show=True, use_cuda=False
        )
    w = pencil.gene_weights(plot=True)
    plt.close()

colors = sns.color_palette('Set1')
num_classes = len(set(adata.obs.condition.values))
pal=colors[0:num_classes]
pal.insert(0, colors[-1]) #add grey


pred_labels = np.array(labels_raw.categories[pred])
pred_labels[confidence < 0] = 'Rejected'
adata.obs['pred_labels'] = pred_labels
adata.obs['confidence_score'] = confidence

sc.pl.umap(adata, color='pred_labels', palette=pal)
sc.pl.umap(
    adata, color='confidence_score', 
    color_map=LinearSegmentedColormap.from_list('mymap', [colors[-1], 'blue'])
    )

