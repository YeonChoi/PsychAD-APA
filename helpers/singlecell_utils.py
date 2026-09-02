# -*- coding: utf-8 -*-
"""
authors: Pengfei Dong, Yeon Choi

modified from pegasus extensions, written by Donghoon Lee:
https://github.com/DiseaseNeuroGenomics/PsychADxD/blob/main/2_taxonomy/pge.py
"""
import numpy as np
import pandas as pd
import h5py
import anndata as ad
from scipy import sparse
from anndata._core.sparse_dataset import SparseDataset
from anndata.experimental import read_elem, write_elem

def get_pseudobulk(adata,colname,mean=False):
    z = pd.get_dummies(adata.obs[colname])
    mat = pd.get_dummies(z).to_numpy()
    if mean:
        mat = mat/mat.sum(axis=0)
    mat = sparse.csr_matrix(mat.T)
    psb = mat.dot(adata.X)
    psb = pd.DataFrame(psb.T.todense(),columns=z.columns,index=adata.var.index)   
    return psb

def get_pseudobulk_from_raw(adata, colname, mean=False):
    z = pd.get_dummies(adata.obs[colname])
    mat = pd.get_dummies(z).to_numpy()
    if mean:                                                                        
        mat = mat/mat.sum(axis=0)
    mat = sparse.csr_matrix(mat.T)
    psb = mat.dot(adata.raw.X) # only difference is here
    psb = pd.DataFrame(psb.T.todense(),columns=z.columns,index=adata.var.index)
    
    return psb

def read_adata_withraw_except_x(pth):
    with h5py.File(pth) as f:
        attrs = list(f.keys())
        attrs.remove('X')
        adata = ad.AnnData(**{k: read_elem(f[k]) for k in attrs})
        print(adata.shape)
    return adata
    
def ondisk_pseudobulk(raw_path,obs,X='X',chunk_size = 500000,obs_columns=['cell_type'],index_col='_index'):
    phi = pd.get_dummies(obs[obs_columns])
    colnames = phi.columns
    phi = sparse.csr_matrix(phi.to_numpy().T)
    with h5py.File(raw_path, 'r') as f:
        csr_indptr = f[X]['indptr'][:]
        gene = byte2utf(f['var'][index_col][:])
    psb = sparse.csr_matrix((phi.shape[0],gene.shape[0]), dtype=np.int64)
    for idx in  range(0, csr_indptr.shape[0]-1  , chunk_size):
        print('Processing', idx, 'to', idx+chunk_size)
        row_start, row_end = idx, idx+chunk_size
        
        with h5py.File(raw_path, 'r') as f:
            tmp_indptr = csr_indptr[row_start:row_end+1]
            new_data = f[X]['data'][tmp_indptr[0]:tmp_indptr[-1]]
            new_indices = f[X]['indices'][tmp_indptr[0]:tmp_indptr[-1]]
            new_indptr = tmp_indptr - csr_indptr[row_start]
            new_shape = [tmp_indptr.shape[0]-1, gene.shape[0]]
            tmp_csr = sparse.csr_matrix((new_data, new_indices, new_indptr), 
                                        shape=new_shape)
            tmp_phi = phi[:,row_start:row_end]
            psb = psb+tmp_phi.dot(tmp_csr)
    psb = pd.DataFrame(psb.toarray().T,
                     columns=colnames,index=gene)
    return psb

def byte2utf(npv):
    if isinstance(npv[0], bytes):
        return np.array([x.decode('utf-8') for x in npv])
    else:
        return npv

def read_obsvar(h5ad_path,group='obs',mod=False):

    def group2vector(vgp):
        if len(vgp.keys()) == 1:
            for i in vgp.keys():
                return byte2utf(vgp[i][:])
        else:
            cat = byte2utf(vgp['categories'][:])
            res = cat[vgp['codes'][:]]
            if any(vgp['codes'][:]<0):
                res[vgp['codes'][:]<0]=np.nan
            return res
    with h5py.File(h5ad_path, mode='r') as h5file:
        M = {}
        if mod:
            fh=h5file['mod'][mod][group]
        else:
            fh = h5file[group]
        for x in fh.keys():
            if isinstance(fh[x], h5py.Dataset):
                M[x] = byte2utf(fh[x][:])
            elif isinstance(fh[x], h5py.Group):
                if 'categories' in fh[x].keys() and fh[x]['categories'][:].shape[0]==0:
                    continue
                M[x]=group2vector(fh[x])
    return pd.DataFrame(M) 


################################################################################################################
# ###############################From Donghoon###########################################################
# ###############################################################################################################

def read_everything_but_X(pth, verbose=True) -> ad.AnnData:
    # read all keys but X and raw
    with h5py.File(pth) as f:
        attrs = list(f.keys())
        attrs.remove('X')
        if 'raw' in attrs:
            attrs.remove('raw')
        adata = ad.AnnData(**{k: read_elem(f[k]) for k in attrs})
        if verbose: print(adata.shape)
    return adata

def concat_on_disk(input_pths, output_pth, pbar=None, read_verbose=True, x_dtype=None):
    """
    Params
    ------
    input_pths
        Paths to h5ad files which will be concatenated
    output_pth
        File to write as a result
    pbar
        tqdm progressbar object (should be initialized outside)
    read_verbose
        print size of each h5ad (default: True)
    """
    annotations = ad.concat([read_everything_but_X(pth, verbose=read_verbose) for pth in input_pths])
    annotations.write_h5ad(output_pth)
    n_variables = annotations.shape[1]

    del annotations

    with h5py.File(output_pth, "a") as target:
        if x_dtype:
            dummy_X = sparse.csr_matrix((0, n_variables), dtype=x_dtype)
        else: # Assumes float X
            dummy_X = sparse.csr_matrix((0, n_variables), dtype=np.float32)
        dummy_X.indptr = dummy_X.indptr.astype(np.int64) # Guarding against overflow for very large datasets
        dummy_X.indices = dummy_X.indices.astype(np.int64) # Guarding against overflow for very large datasets

        write_elem(target, "X", dummy_X)
        mtx = SparseDataset(target["X"])
        for p in input_pths:
            with h5py.File(p, "r") as src:
                mtx.append(SparseDataset(src["X"]))
                if pbar:
                    pbar.update(n=1)
                
def write_h5ad_with_new_annotation(original_h5ad, adata, new_h5ad, raw = False):
    # new annotation
    new_uns=None
    if adata.uns:
        new_uns = adata.uns
    new_obsm=None
    if adata.obsm:
        new_obsm = adata.obsm
    new_varm=None
    if adata.varm:
        new_varm = adata.varm
    new_obsp=None
    if adata.obsp:
        new_obsp = adata.obsp
    new_varp=None
    if adata.varp:
        new_varp = adata.varp

    # save obs and var first
    ad.AnnData(None, obs=adata.obs, var=adata.var, uns=new_uns, obsm=new_obsm, varm=new_varm, obsp=new_obsp, varp=new_varp).write(new_h5ad)

    # append X
    with h5py.File(new_h5ad, "a") as target:
        # make dummy
        dummy_X = sparse.csr_matrix((0, adata.var.shape[0]), dtype=np.float32)
        dummy_X.indptr = dummy_X.indptr.astype(np.int64) # Guarding against overflow for very large datasets
        dummy_X.indices = dummy_X.indices.astype(np.int64) # Guarding against overflow for very large datasets
        
        with h5py.File(original_h5ad, "r") as src:
            write_elem(target, "X", dummy_X)
            SparseDataset(target["X"]).append(SparseDataset(src["X"]))
            # append raw/X if needed
            if raw:
                write_elem(target, "raw/X", dummy_X)
                SparseDataset(target["raw/X"]).append(SparseDataset(src["raw/X"]))

def clean_unused_categories(data):
    for obs_name in data.obs.columns:
        if data.obs[obs_name].dtype=='category':
            print('Removing unused categories from',obs_name)
            data.obs[obs_name] = data.obs[obs_name].cat.remove_unused_categories()
    for var_name in data.var.columns:
        if data.var[var_name].dtype=='category':
            print('Removing unused categories from',var_name)
            data.var[var_name] = data.var[var_name].cat.remove_unused_categories()
    return data

def ondisk_subset(orig_h5ad, new_h5ad, subset_obs, subset_var = None, chunk_size = 500000, raw = False, adata = None):

    if adata is None:
        
        # read annotations only
        adata = read_everything_but_X(orig_h5ad)
        
        # subset annotation
        if subset_var is not None:
            adata = adata[subset_obs,subset_var]
        else:
            adata = adata[subset_obs,:]

        # clean unused cat
        adata = clean_unused_categories(adata)
        
    # new annotation
    new_uns=None
    if adata.uns:
        new_uns = adata.uns

    new_obsm=None
    if adata.obsm:
        new_obsm = adata.obsm

    new_varm=None
    if adata.varm:
        new_varm = adata.varm

    new_obsp=None
    if adata.obsp:
        new_obsp = adata.obsp

    new_varp=None
    if adata.varp:
        new_varp = adata.varp
    
    # save obs and var first
    ad.AnnData(None, obs=adata.obs, var=adata.var, uns=new_uns, obsm=new_obsm, varm=new_varm, obsp=new_obsp, varp=new_varp).write(new_h5ad)
    
    # initialize new_h5ad
    with h5py.File(new_h5ad, "a") as target:
        dummy_X = sparse.csr_matrix((0, adata.var.shape[0]), dtype=np.float32)
        dummy_X.indptr = dummy_X.indptr.astype(np.int64) # Guarding against overflow for very large datasets
        dummy_X.indices = dummy_X.indices.astype(np.int64) # Guarding against overflow for very large datasets
        write_elem(target, "X", dummy_X)
        if raw:
            write_elem(target, "raw/X", dummy_X)
        
    # get indptr first
    with h5py.File(orig_h5ad, 'r') as f:
        csr_indptr = f['X/indptr'][:]

    # append subset of X
    for idx in [i for i in range(0, csr_indptr.shape[0]-1, chunk_size)]:
        print('Processing', idx, 'to', idx+chunk_size)
        row_start, row_end = idx, idx+chunk_size

        if sum(subset_obs[row_start:row_end])>0:
            # X
            with h5py.File(orig_h5ad, 'r') as f:
                tmp_indptr = csr_indptr[row_start:row_end+1]
                
                new_data = f['X/data'][tmp_indptr[0]:tmp_indptr[-1]]
                new_indices = f['X/indices'][tmp_indptr[0]:tmp_indptr[-1]]
                new_indptr = tmp_indptr - csr_indptr[row_start]
                new_shape = [tmp_indptr.shape[0]-1, adata.shape[1]]
                                
                tmp_csr = sparse.csr_matrix((new_data, new_indices, new_indptr), shape=new_shape)

                if subset_var is not None:
                    tmp_csr = tmp_csr[subset_obs[row_start:row_end]][:,subset_var]
                else:
                    tmp_csr = tmp_csr[subset_obs[row_start:row_end]]

                tmp_csr.sort_indices()

            # append X
            with h5py.File(new_h5ad, "a") as target:
                mtx = SparseDataset(target["X"])
                mtx.append(tmp_csr)

            # raw/X
            if raw:
                with h5py.File(orig_h5ad, 'r') as f:
                    tmp_indptr = csr_indptr[row_start:row_end+1]
                    
                    new_data = f['raw/X/data'][tmp_indptr[0]:tmp_indptr[-1]]
                    new_indices = f['raw/X/indices'][tmp_indptr[0]:tmp_indptr[-1]]
                    new_indptr = tmp_indptr - csr_indptr[row_start]
                    new_shape = [tmp_indptr.shape[0]-1, adata.shape[1]]
                                    
                    tmp_csr = sparse.csr_matrix((new_data, new_indices, new_indptr), shape=new_shape)

                    if subset_var is not None:
                        tmp_csr = tmp_csr[subset_obs[row_start:row_end]][:,subset_var]
                    else:
                        tmp_csr = tmp_csr[subset_obs[row_start:row_end]]

                    tmp_csr.sort_indices()

                # append raw/X
                with h5py.File(new_h5ad, "a") as target:
                    mtx = SparseDataset(target["raw/X"])
                    mtx.append(tmp_csr)


