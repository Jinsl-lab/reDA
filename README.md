# reDA: differential abundance testing on scATAC-seq data using random walk with restart
 reDA is a differential abundance test framework based on random walk with restart, for scATAC-seq data. To better measure the abundance of cells under different conditions, reDA introduces a random walk with restart, which can better capture the local and global information of the shared nearest neighbor (SNN) graph, mitigate the effects of information loss and information redundancy, and thus reliably identify condition-specific cell subpopulations based on association tests.
## How to install
Before installing reDA, please install the following packages manually.
### Python Dependencies
reDA depends on the following Python packages:
```
numpy 1.26.2
scipy 1.11.3
pandas 1.3.5
argparse 1.4.0
scanpy 1.9.5
packaging 23.2
anndata 0.10.2
rpy2 3.4.5
multianndata 0.0.4
```
### R Dependencies
reDA depends on the following R packages:
```
Signac 1.12.9001
Seurat 5.0.1
harmony 0.1.1
readr 2.1.5
GenomicRanges 1.50.2
SeuratDisk 0.0.0.9020
```
### Installation
reDA is developed under Python (version >= 3.9) and R (version >= 4.2). To use reDA, you can clone this repository:
```
git clone https://github.com/Jinsl-lab/reDA.git 
cd reDA
```
Then install it.
```
pip install .
```
## Quick start
```python
import numpy as np
import pandas as pd
import scanpy as sc
import reDA

#pre-processing
reDA.pp.preprocess(rpath="rpyprc.R",  
                   inpath="scATAC_data.rds",  
                   hadpath="output_data.h5Seurat", 
                   group="dataset", 
                   assay="ATAC",
                   dimsta=2,dimend=30)

#run
res = reDA.tl.association(d,                   
            d.samplem.phenotype,  
            re=0.3,                         
            covs=d.samplem[['covs']])
```
## Tutorials
We also provide a [tutorial]() for running reDA. The test dataset is available here.
