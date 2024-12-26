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
from .NAM import nam, _df_to_array

#some are adapted from https://github.com/immunogenomics/cna
#These functions are based on https://github.com/immunogenomics/cna

def _association(NAMsvd, M, r, y, ks=None, Nnull=1000, force_permute_all=False,
                    local_test=True, seed=None):
    if seed is not None:
        np.random.seed(seed)
    
    batches = np.ones(len(y))

    # prep data
    (U, sv, V) = NAMsvd
    y = (y - y.mean())/y.std()   
    n = len(y)

    if ks is None:
        incr = max(int(0.02*n), 1)
        maxnpcs = min(4*incr, int(n/5))
        ks = np.arange(incr, maxnpcs+1, incr)

    def _reg(q, k):
        Xpc = U[:,:k]
        beta = Xpc.T.dot(q) 
        qhat = Xpc.dot(beta)
        return qhat, beta

    def _stats(yhat, ycond, k):
        ssefull = (yhat - ycond).dot(yhat - ycond)
        ssered = ycond.dot(ycond)
        deltasse =  ssered - ssefull
        f = (deltasse / k) / (ssefull/n)
        p = st.f.sf(f, k, n-(1+r+k)) # F test
        r2 = 1 - ssefull/ssered
        return p, r2

    def _minp_stats(z):
        zcond = M.dot(z)
        zcond = zcond / zcond.std()
        ps, r2s = np.array([
            _stats(
                _reg(zcond, k)[0],
                zcond,
                k)
            for k in ks
        ]).T
        k_ = np.argmin(ps)
        return ks[k_], ps[k_], r2s[k_]

    # get non-null f-test p-value
    k, p, r2 = _minp_stats(y)
    if k == max(ks):
        warnings.warn(('data supported use of {} NAM PCs, which is the maximum considered. '+\
            'Consider allowing more PCs by using the "ks" argument.').format(k))

    # compute coefficients and r2 with chosen model
    ycond = M.dot(y)
    ycond /= ycond.std()
    yhat, beta = _reg(ycond, k)
    r2_perpc = (beta / np.sqrt(ycond.dot(ycond)))**2

    # get neighborhood scores with chosen model
    ncorrs = (np.sqrt(sv[:k])*beta/n).dot(V[:,:k].T)   

    # compute final p-value using Nnull null f-test p-values
    y_ = conditional_permutation(batches, y, Nnull)
    nullminps, nullr2s = np.array([_minp_stats(y__)[1:] for y__ in y_.T]).T
    pfinal = ((nullminps <= p+1e-8).sum() + 1)/(Nnull + 1)
    if (nullminps <= p+1e-8).sum() == 0:
        warnings.warn('global association p-value attained minimal possible value. '+\
                'Consider increasing Nnull')

    # get neighborhood fdrs if requested
    fdrs, fdr_5p_t, fdr_10p_t = None, None, None
    if local_test:
        print('computing neighborhood-level FDRs')
        Nnull = min(1000, Nnull)
        y_ = y_[:,:Nnull]
        ycond_ = M.dot(y_)
        ycond_ /= ycond_.std(axis=0)
        gamma_ = U[:,:k].T.dot(ycond_)
        nullncorrs = np.abs(V[:,:k].dot(np.sqrt(sv[:k])[:,None]*(gamma_ / n))) 

        maxcorr = np.abs(ncorrs).max()
        fdr_thresholds = np.arange(maxcorr/4, maxcorr, maxcorr/400)
        fdr_vals = empirical_fdrs(ncorrs, nullncorrs, fdr_thresholds)

        fdrs = pd.DataFrame({
            'threshold':fdr_thresholds,
            'fdr':fdr_vals,
            'num_detected': [(np.abs(ncorrs)>t).sum() for t in fdr_thresholds]})

        # find maximal FDR<5% and FDR<10% sets
        if np.min(fdrs.fdr)>0.05:
            fdr_5p_t = None
        else:
            fdr_5p_t = fdrs[fdrs.fdr <= 0.05].iloc[0].threshold
        if np.min(fdrs.fdr)>0.1:
            fdr_10p_t = None
        else:
            fdr_10p_t = fdrs[fdrs.fdr <= 0.1].iloc[0].threshold

        del gamma_, nullncorrs

    del y_

    res = {'p':pfinal, 'nullminps':nullminps, 'k':k, 'ncorrs':ncorrs, 'fdrs':fdrs,
            'fdr_5p_t':fdr_5p_t, 'fdr_10p_t':fdr_10p_t,
			'yresid_hat':yhat, 'yresid':ycond, 'ks':ks, 'beta':beta,
            'r2':r2, 'r2_perpc':r2_perpc,
            'nullr2_mean':nullr2s.mean(), 'nullr2_std':nullr2s.std()}
    return Namespace(**res)

def association(data, y, re=0.4, covs=None, nsteps=None, suffix='',
    force_recompute=False, **kwargs):

    # formatting and error checking
    covs = _df_to_array(data, covs)
    y = _df_to_array(data, y)
    if y.shape != (data.N,):
        raise ValueError(
            'y should be an array of length data.N; instead its shape is: '+str(y.shape))

    if covs is not None:
        filter_samples = ~(np.isnan(y) | np.any(np.isnan(covs), axis=1))
    else:
        filter_samples = ~np.isnan(y)

    du = data.uns
    nam(data, re=re, covs=covs, filter_samples=filter_samples,
                    nsteps=nsteps, suffix=suffix,
                    force_recompute=force_recompute, **kwargs)
    NAMsvd = (
        du['NAM_sampleXpc'+suffix].values,
        du['NAM_svs'+suffix],
        du['NAM_nbhdXpc'+suffix].values
        )

    print('performing association test')
    res = _association(NAMsvd, du['_M'+suffix], du['_r'+suffix],
        y[du['_filter_samples'+suffix]], **kwargs)

    # add info about kept cells
    vars(res)['kept'] = du['keptcells'+suffix]

    return res

#These functions come from https://github.com/immunogenomics/cna


def conditional_permutation(B, Y, num):
    """
    Permutes Y conditioned on B num different times.
    """
    batchind = np.array([
        np.where(B == b)[0] for b in np.unique(B)
        ])
    
    ix = np.concatenate([
        bi[np.argsort(np.random.randn(len(bi), num), axis=0)]
        for bi in batchind
        ])

    bix = np.zeros((len(Y), num)).astype(int)
    bix[np.concatenate(batchind)] = ix
    result = Y[bix]
    return Y[bix]

def empirical_fdrs(z, znull, thresholds):
    """
    znull is assumed to be of shape len(z) x k, where k is the number of
        null instantiations.
    """
    # get tail counts
    if znull.shape[0] != len(z):
        print('ERROR: znull is shape', znull.shape, 'and z is shape', z.shape)
    tails = tail_counts(thresholds, znull)
    ranks = tail_counts(thresholds, z)


    # compute FDPs
    fdp = tails / ranks
    fdr = fdp.mean(axis=0)
    
    return fdr


def tail_counts(z, znull, atol=1e-8, rtol=1e-5):
    """
    Computes the number of null z-scores of equal or greater magnitude than each supplied
        z-score, for each null instantiation.

    znull is assumed to be either of shape len(z) or len(z) x k, where k is the number of
        null instantiations.
    """
    # re-shape znull if only one set of null results is supplied
    if len(znull.shape) == 1:
        znull = znull.reshape((-1, 1))

    # square z-scores and prepare to sort/un-sort them
    z2 = z**2
    ix = np.argsort(z2); iix = np.argsort(ix)

    # ask numpy to make a histogram with sorted squared z-scores as boundaries
    bins = np.concatenate([z2[ix] - atol - rtol*z2[ix], [np.inf]])
    hist = np.array([
            np.histogram(zn2, bins=bins)[0]
        for zn2 in znull.T**2])

    # convert histogram into tail counts (in-place)
    tails = np.flip(hist, axis=1)
    np.cumsum(tails, axis=1, out=tails)
    tails = np.flip(tails, axis=1)

    # return tail counts for un-sorted z-scores
    return tails[:, iix]