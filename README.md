# OVERVIEW

Git Repository for 20.440 final project, Analysis of cochlear hair cell reprogramming in mice. 

One of the major causes of hearing loss is damage to the hair cells in the inner ear.  These cells cannot regenerate, and therefore any damage is permanent.  Studying the gene expression cochlear cells reprogrammed into hair cells could reveal potential targets for intervention upon cochlear hair cell damage.
	
# DATA
scRNA seq of ATOH1, GFI1, and POUF4 reprogrammed cochlear hair cells from 3 genetically engineered mouse models (Iyer et al.).  Mice were bred to create 3 conditional overexpression lines (ATOH1, ATOH1+GFI1, and ATOH1+GFI1+POUF4). Tamoxofin was injected at postnatal day 1 or day 8  (P1/P8) to stimulate reprogramming and cochlear cells were purified on P8 or P15 with FACS (1 week after tamoxofin injection).  These cells were used for scRNA seq. Lastly two additional datasets are available: simulatneous scRNA ATAC seq was done on P1 wild type mouse cochlear cells and P8 wild type mouse cochlear cells.

Raw data is available at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182202

**P1-P8 scRNA seq**

- WT P1-P8 Tamoxofin: GSM5520356
- ATOH1 P1-P8 Tamoxofin: GSM5520357
- ATOH1+GFI1 P1-P8 Tamoxofin: GSM5520358
- ATOH1+GFI1+POUF4 P1-P8 Tamoxofin: GSM5520359


**P8-P15 scRNA seq**

- WT P8-P15 Tamoxofin: GSM5520360
- ATOH1 P8-P15 Tamoxofin: GSM5520361
- ATOH1+GFI1 P8-P15 Tamoxofin: GSM5520362
- ATOH1+GFI1+POUF4 P8-P15 Tamoxofin: GSM5520363

**P1 Multiome**

- P1 scRNA seq: GSM6883295
- P1 ATAC seq: GSM6883296

**P8 Multiome**

- P8 scRNA seq: GSM6883297
- P8 ATAC seq: GSM6883298



# FOLDER STRUCTURE
	
- **src/** : 	directory containing all scripts used from preprocessing to figure generation
	- **data/** :  directory containing data manipulation scripts
	- **analysis/** : directory all analysis scripts
	- **vis/** : directory with all visualization scripts
- **data/** : 	directory containing data used in analysis
    - **raw/** : raw counts and ATAC data (not included in repo due to space constraints)
- **fig/** : 	directory containing outputs from src/vis
- **notebook/** : 	directory containing notebooks used in initial data exploration 

# INSTALLATION

Analysis done on OSX 13.2.1

1. See requirements.txt for library requirements. 
    - run "pip install -r requirements.txt" to install all requirements

2. download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182202 and place into **data/raw/** directory (create if needed)
    - NOTE: for UMAP plot (figure 1) this can be done with command:  "sh figures.sh UMAP"
    - NOTE: for UMAP plot (figure 1) only GSM5520356  mtx file is needed 
    - unzip *.tar.gz to get experiment directory
    - in experiment directory unzip *.tsv.gz and *.mtx.gz to get individual data files
    



3. Use figures.sh to generate all figures.



# CITATIONS

1. Amrita A Iyer, Ishwar Hosamani, John D Nguyen, Tiantian Cai, Sunita Singh, Melissa M McGovern, Lisa Beyer, Hongyuan Zhang, Hsin-I Jen, Rizwan Yousaf, Onur Birol, Jenny J Sun, Russell S Ray, Yehoash Raphael, Neil Segil, Andrew K Groves (2022) Cellular reprogramming with ATOH1, GFI1, and POU4F3 implicate epigenetic changes and cell-cell signaling as obstacles to hair cell regeneration in mature mammals eLife 11:e79712, https://doi.org/10.7554/eLife.79712