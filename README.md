# OVERVIEW

Git Repository for 20.440 final project, Analysis of cochlear hair cell reprogramming in mice. 

The major cause of hearing loss is irreversible damage to sensory hair cells in the inner ear. One active research direction for hearing restoration focuses on the reprogramming of cochlear cells into functional hair cells. To further explore the effects of transcription factor-mediated reprogramming, we reanalyze a public scRNA-seq dataset containing combinations of three different transcription factor treatments on nascent and mature mouse cochlear cells by performing gene ontology analysis and trajectory inference. The analysis revealed that potentially detrimental off-target effects are exhibited by critical cochlear cell types and suggested possibilities to improve treatment efficacy.
	
# DATA
scRNA seq of ATOH1, GFI1, and POUF4 reprogrammed cochlear hair cells from 3 genetically engineered mouse models (Iyer et al.).  Mice were bred to create 3 conditional overexpression lines (ATOH1, ATOH1+GFI1, and ATOH1+GFI1+POUF4). Tamoxofin was injected at postnatal day 1 or day 8 to stimulate reprogramming and cochlear cells were purified on at day 8 and day 15 with FACS (1 week after tamoxofin injection, P8/P15).  These cells were used for scRNA seq.

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


# FOLDER STRUCTURE
	
- **src/** : 	directory containing all scripts used from preprocessing to figure generation
	- **data/** :  directory containing data manipulation scripts
	- **analysis/** : directory all analysis scripts
	- **vis/** : directory with all visualization scripts
- **data/** : 	directory containing data used in analysis
    - **raw/** : raw counts data (not included in repo due to space constraints)
    - **processed/** : Cell Cluster marker genes (Kolla et al.), and processed data generated in src/data/ (mostly excluded due to space constraints)
- **fig/** : 	directory containing outputs from src/vis and src/analysis (excluded due to space constraints)
- **notebook/** : 	directory containing notebooks used in initial data exploration 

# INSTALLATION

Analysis done on OSX 13.2.1

1. See requirements.txt for library requirements. 
    - install all required packages in R

2. use "sh figures.sh Data" to download data, run preprocessing and generate figures
    - figures can also be generated by manually runnning R scripts found in src/vis/ and src/analysis/

OR 

2. download data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE182202 and place into **data/raw/** directory (create if needed)
    - unzip *.tar.gz to get experiment directory
    - in experiment directory unzip *.tsv.gz and *.mtx.gz to get individual data files
    - then run "src/data/Preprocess-Integration.R", "src/data/Cluster_Labelling.R" and "src/data/Cluster_Annotation.R" to generate processed data
    - finally generate figures with "sh figures.sh"

NOTE: aside from figures
-  all files associated with cluster annotation are generated in data/processed/Celltype_Annotation/
-  inputs to Gene Ontology are generated in data/processed/Ontology_Input/

# CITATIONS

1. Amrita A Iyer, Ishwar Hosamani, John D Nguyen, Tiantian Cai, Sunita Singh, Melissa M McGovern, Lisa Beyer, Hongyuan Zhang, Hsin-I Jen, Rizwan Yousaf, Onur Birol, Jenny J Sun, Russell S Ray, Yehoash Raphael, Neil Segil, Andrew K Groves (2022) Cellular reprogramming with ATOH1, GFI1, and POU4F3 implicate epigenetic changes and cell-cell signaling as obstacles to hair cell regeneration in mature mammals eLife 11:e79712, https://doi.org/10.7554/eLife.79712

2. L. Kolla, M. C. Kelly, Z. F. Mann, A. Anaya-Rocha, K. Ellis, A. Lemons, A. T. Palermo, K. S. So, J. C. Mays, J. Orvis, J. C. Burns, R. Hertzano, E. C. Driver, and M. W. Kelley. Character- ization of the development of the mouse cochlear epithelium at the single cell level. Nature Communications, 11(1):2389, May 2020.
