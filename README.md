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
  
**3.** From within Rstudio open */scripts/NE04_Seurat_analysis.Rmd* and start from importing data tables and metadata to create a Seurat object. Your options are to load: 

*i*. NE cell adult data with doublets included *(GC_table_adult_NE.tsv, Metadata_table_adult_NE.tsv)*

*ii*. NE cell adult data set without doublets *(GC_table_adult_NE_ND.tsv, Metadata_table_adult_NE_ND.tsv)*

*iii*. NE cell adult data together with all other datasets (+Kransnow,+TabulaMuris) *(GC_table_ALL_datasets.tsv, Metadata_ALL_datasets.tsv)*

*iv*. If you want to directly import the adult NE cells that we have been already analyzing with the old pipeline, use the objects
that contain the append "from_old". These would be: *GC_table_adult_NE_from_old, GC_table_adult_NE_from_old_ND, GC_table_ALL_datasets_ND* along with the corresponding metadata tables. 

**Note1.** Scrublet doublet analysis for the NE adult data is performed in *NE02_Scrublet_raw.Rmd* and *NE_scrublet_basics.ipynb*
