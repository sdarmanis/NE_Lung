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

**3.** From within Rstudio open */scripts/NE04_Seurat_analysis.Rmd* and start from importing data tables and metadata to create a Seurat object. Your options are to load any of the files as explained in **Note 1**. 


**Note 1.** If you are familiar with Seurat or other packages used to analyze single cell data then you can create the objects required from raw data tables (raw GC tables) and a basic metadata sheet. There are different versions of those within input and they are explained below: 
*"GC/Metadata"_table_zsGreen_renamed*: All zsGreen cells (fetal and adult) after renaming (pertains to all "renamed" files) to account for incorrect demuux sheet.
*"GC/Metadata"_table_adult_NE_renamed.tsv*: Same as above, only adult cells
*"GC/Metadata"_table_adult_NE_renamed_after_Scrublet.tsv*: Same as above, including scrublet analysis metadata fields
*"GC/Metadata"_table_adult_NE_renamed_ND.tsv*: Same as above, with scrublet identified doublets removed
*"GC/Metadata"_table_ALL_datasets_renamed_ND.tsv*: Same as above, with Tabula Muris and Krasnow datasets. Those have not been filtered on doublets. 

**Note 2.** Scrublet doublet analysis for the NE adult data is performed in *NE02_Scrublet_raw.Rmd* and *NE_scrublet_basics.ipynb*
