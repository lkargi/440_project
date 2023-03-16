import pandas as pd
from scipy.io import mmread
import numpy as np
import os

def QC(data_dir, gene_count_thresh = 1500 ):
    """
    Perform QC on raw counts data by filtering cells with too few genes

    PARAM
        OUT:
         filt - filtered counts np.array (cells, genes)

        IN:
         data_dir - directory within data/raw/ to pull from 
         gene_count_thresh - number of genes required to keep cell
    """
    #Import Data
    data = 'data/raw/'+data_dir+'/matrix.mtx'
    mtx = mmread(data)
    counts = np.array(mtx.todense()).T

    #Filter out cells with too few genes
    counts_filt = np.sum(counts != 0, axis = 1 ) >gene_count_thresh
    return counts[counts_filt,:]

def log_normalize(filt, scale = 10000):
    """
    Normalize to total number of molecules per cell then scale and add 
    1 for log normalization

    PARAM
        OUT:
         norm - log normalized expression values as np.array (cells, genes)
        
        IN:
            filt - filtered counts np.array (cells, genes)
            scale - value to scale expression before log
    """
    tot_mol = np.sum(filt, axis = 1)
    return np.log10( np.divide(filt.T, tot_mol).T * scale + 1)

if __name__ == "__main__":
    data_dir = "s9_rosa26WT_p8"
    filt = QC(data_dir)
    norm = log_normalize(filt)
    if not os.path.isdir('data/filt_norm/'+data_dir):
        os.mkdir('data/filt_norm/'+data_dir)
    np.save('data/filt_norm/'+data_dir+'/norm.npy')
