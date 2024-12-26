import rpy2.robjects as robjects
from rpy2.robjects import numpy2ri

numpy2ri.activate()

#Calculate the connectivity matrix、umap and save as mtx format
def preprocess(rpath,inpath,hadpath,group,assay="ATAC",dimsta=2,dimend=50,issave="NO",rdspath="NO"):
    robjects.r.source(rpath) #path to prc.R
    r2py = robjects.r.prc(inpath,dimsta,dimend,group,assay,hadpath,issave,rdspath)
