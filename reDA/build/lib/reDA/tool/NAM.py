import numpy as np
import scipy.sparse as sp
import pandas as pd
from argparse import Namespace
import scanpy as sc
from scipy.sparse.linalg import inv  
from scipy.sparse import csr_matrix,csc_matrix
from scipy.sparse.linalg import spsolve
import gc, warnings
import scipy.stats as st
from packaging import version
import anndata

#some are adapted from https://github.com/immunogenomics/cna
#These functions are based on https://github.com/immunogenomics/cna

def nam(data, re=0.4, covs=None, filter_samples=None,
    nsteps=None, max_frac_pcs=0.15, suffix='', ks = None,
    force_recompute=False, **kwargs):     
    def safe_same(A, B):
        if A is None: A = np.zeros(0)
        if B is None: B = np.zeros(0)
        if A.shape == B.shape:
            return np.allclose(A, B, equal_nan = True)
        else:
            return False

    # error checking
    covs = _df_to_array(data, covs)
    if filter_samples is None:
        if covs is not None:
            filter_samples = ~np.any(np.isnan(covs), axis=1)
        else:
            filter_samples = np.repeat(True, data.N)
    else:
        filter_samples = _df_to_array(data, filter_samples)

    du = data.uns
    # compute NAM
    if force_recompute or \
        'NAM.T'+suffix not in du:
        print('NAM not found; computing and saving')
        NAM = _nam(data, nsteps=nsteps, re=re)  
        NAMqc, keep = NAM.values, np.repeat(True, len(NAM.values.T)) 
        du['NAM.T'+suffix] = pd.DataFrame(NAMqc, index=NAM.index, columns=NAM.columns[keep]).T
        du['keptcells'+suffix] = keep
        

    # correct for covariates and SVD
    if force_recompute or \
        'NAM_sampleXpc'+suffix not in du or \
        not safe_same(covs, du['_covs'+suffix]) or \
        not safe_same(filter_samples, du['_filter_samples'+suffix]):

        print('covariate-adjusted NAM not found; computing and saving')
        du['_filter_samples'+suffix] = filter_samples
        NAM = du['NAM.T'+suffix].T.iloc[filter_samples]
        NAM_resid, M, r = _resid_nam(NAM.values,
                                covs[filter_samples] if covs is not None else covs)

        print('computing SVD')
        U, svs, V = _svd_nam(NAM_resid)
        npcs = min(V.shape[1], max([10]+[int(max_frac_pcs * data.N)]+[ks if ks is not None else []][0]))
        du['NAM_sampleXpc'+suffix] = pd.DataFrame(U,
            index=NAM.index,
            columns=['PC'+str(i) for i in range(1, len(U.T)+1)])
        du['NAM_svs'+suffix] = svs
        du['NAM_varexp'+suffix] = svs / len(U) / len(V)
        du['NAM_nbhdXpc'+suffix] = pd.DataFrame(V[:,:npcs],
            index=NAM.columns,
            columns=['PC'+str(i) for i in range(1, npcs+1)])
        du['_M'+suffix] = M
        du['_r'+suffix] = r
        du['_covs'+suffix] = (np.zeros(0) if covs is None else covs)


def diffuse_stepwise(data, s, maxnsteps=15,re=0.4):   
    # find connectivity matrix
    av = anndata.__version__   
    
    if type(av) == str:
        av = version.parse(av)
    if av < version.parse("0.7.2"):
        a = data.uns["neighbors"]['distances'] 
    else:
        a = data.obsp['distances'] 
    #random walk with restart
    length = a.shape[0] 
    rowsum = np.array(a.sum(axis=1)).flatten() 
    
    D = sp.diags(rowsum,0) 
    
    W = inv(D.tocsc()).dot(a) 
    
    l = sp.diags([1] * length, 0) 
    
    F = l 
    ind = s.index 
    col = s.columns
    # do diffusion
    for i in range(maxnsteps):
        print('\ttaking step', i+1)                      
        F = (1-re)*W.dot(F)+re*l                                          
        
        mid = (F.T).dot(s)    
        mid = pd.DataFrame(mid,index=ind,columns=col) 
        
        yield mid 




def _nam(data, nsteps=None, maxnsteps=15, re=0.4): 
    def R(A, B):
        return ((A - A.mean(axis=0))*(B - B.mean(axis=0))).mean(axis=0) \
            / A.std(axis=0) / B.std(axis=0)

    S = pd.get_dummies(data.obs_sampleids)[data.samplem.index.values] 
    

    prevmedkurt = np.inf 

    old_s = np.zeros(S.shape) 
    
    for i, s in enumerate(diffuse_stepwise(data, S, maxnsteps=maxnsteps, re=re)): 
        
        medkurt = np.median(st.kurtosis(s/s.sum(axis=0), axis=1)) 
        R2 = R(s, old_s)**2  
        old_s = s
        print('\tmedian kurtosis:', medkurt+3)
        print('\t20th percentile R2(t,t-1):', np.percentile(R2, 20))
        if nsteps is None:   
            if prevmedkurt - medkurt < 3 and i+1 >= 2 and prevmedkurt > medkurt:    
                print('stopping after', i+1, 'steps')
                break  
            if medkurt+3 < 8:
                print('stopping after', i+1, 'steps')
                break
            prevmedkurt = medkurt
        elif i+1 == nsteps:
            break
        gc.collect() 
    
    
    
    snorm = (s / s.sum(axis=0)).T 
    snorm.index.name = data.samplem.index.name          
    return snorm



# residualizes covariates information out of NAM
def _resid_nam(NAM, covs):
    N = len(NAM)
    NAM_ = NAM - NAM.mean(axis=0)
    if covs is None:
        covs = np.ones((N, 0))
    else:
        covs = (covs - covs.mean(axis=0))/covs.std(axis=0)


    C = covs
    if len(C.T) == 0:
        M = np.eye(N)
    else:
        M = np.eye(N) - C.dot(np.linalg.solve(C.T.dot(C), C.T))
    NAM_ = M.dot(NAM_)

    return NAM_ / NAM_.std(axis=0), M, len(C.T)


#These functions come from https://github.com/immunogenomics/cna
#performs SVD of NAM
def _svd_nam(NAM):
    U, svs, UT = np.linalg.svd(NAM.dot(NAM.T))
    V = NAM.T.dot(U) / np.sqrt(svs)

    return (U, svs, V)



def _df_to_array(data, x):
    if type(x) in [pd.DataFrame, pd.Series]:
        if all(x.index == data.samplem.index):
            return x.values
        else:
            print('ERROR: index does not match index of data.samplem')
    else:
        return x