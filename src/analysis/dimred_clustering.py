import numpy as np
import sys
from sklearn.decomposition import PCA

sys.path.insert(0, 'src/data/')
from scRNA_qc_filt import QC, log_normalize


def most_variable_genes(norm, num_genes = 2000):
    """
    Filter normalized data down to the 'num_genes' most variable genes for dimred and clustering

    PARAM
        OUT:
         vargenes - filtered expression including only num_genes genes (cells, genes)
        
        IN:
         norm - log normalized expression matrix (cells, genes)
         num_genes - number of genes to keep
    """
    features_var = norm.var(axis = 0)
    sort_idx = np.argsort(np.array(features_var))
    return norm[:,sort_idx[::-1][0:num_genes]]


if __name__ == "__main__":
    data_dir = "s9_rosa26WT_p8"
    filt = QC(data_dir)
    norm = log_normalize(filt)
    vargenes = most_variable_genes(norm)
