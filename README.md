# Analysis of adult mouse Neuroendocrine cells in collaboration with Christin Kuo's lab

## Instructions

**1.** Clone this repo. On your terminal, navigate to a directory of your choosing *<dir>*
  
 ```
 cd <dir>
 git clone https://github.com/sdarmanis/NE_Lung.git
 ```

**2.** Download Input files from [here](https://drive.google.com/open?id=1zhrR7cPI8OeqMhNW5IgJGBNwam_Xs3M0) and move the .zip file inside the repo directory (which should be called *NE_Lung*). Unzip and remove .zip file. Now you should have a folder structure that looks like *<dir>/NE_Lung/*. In NE_Lung you should have: *input/* *plot_out/* *scripts/* *processed_input/*

 ```
 cd <dir>/NE_Lung
 unzip input.zip
 rm input.zip
 ```
  
**3.** From within Rstudio open */scripts/NE01_Import_and_create_Seurat.Rmd* and start from importing all zsGreen cells, creating metadata, Seurat objects, etc.. 

**4.** Perform the Scrublet doublet analysis (optional, but if not performed N03 needs to be adjusted to not remove doublets) using */scripts/NE02_Scrublet.Rmd* and */scripts/NE_scrublet_basics.ipynb*. The latter is a Jupyter notebook.

**5.** To continue with the analysis of the adult cells use script */scripts/NE03_Create_seurat_cluster_annotate.Rmd*. Remove annotated doublets if you so like or ignore doublet assignments and continue with the analysis of all cells (see instructions within the script).








