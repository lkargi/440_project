import umap
import numpy as np
from sklearn.decomposition import PCA
import sys
import matplotlib.pyplot as plt

sys.path.insert(0, 'src/data/')
from scRNA_qc_filt import QC, log_normalize

sys.path.insert(0, 'src/analysis/')
from dimred_clustering import most_variable_genes

def plot_UMAP(vargenes, PCA_dim = 10, data = "Wild Type"):
    """
       generates UMAP plot of vargenes expression data

        PARAM
            OUT:
            vargenes - filtered expression including only num_genes genes (cells, genes)
            
            IN:
            vargenes - log normalized expression matrix after filtering most variable 
                       genes (cells, genes)
            PCA_dim - number of PCA dimensions to use in UMAP
            data - data label to add to title of plot
        """

    #Run PCA to reduce dimensions of vargenes before UMAP
    pca = PCA()
    exp_pca = pca.fit_transform(vargenes)

    #Run UMAP on first 'PCA_dim' dimensions of PCA
    reduce = umap.UMAP()
    exp_umap = reduce.fit_transform(exp_pca[:,:PCA_dim])

    #Plotting UMAP
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(exp_umap[:,0],exp_umap[:,1],'o')
    
    ax.set_ylabel('UMAP2')
    ax.set_xlabel('UMAP1')
    ax.set_title(f'UMAP of {data} scRNA seq data')
    fig.savefig('fig/UMAP.png')


if __name__ == "__main__":
    data_dir = "s9_rosa26WT_p8"
    filt = QC(data_dir)
    norm = log_normalize(filt)
    vargenes = most_variable_genes(norm)
    plot_UMAP(vargenes)

