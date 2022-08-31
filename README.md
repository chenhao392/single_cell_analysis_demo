# Single_cell_analysis_demo
This is a simple demo for jointly modeling single cell ATAC-seq and RNA-seq datasets using different methods, such as Signac from Seurat package, MultiVI from scVI framework and GLUE framework. Please note that It's for demonstrating working code/pipeline only.

## Sample dataset
The sample dataset is downloaded from the NCBI GEO database (GSE151302). After decompression, the files are stored in a working folder. The demo codes just process one sample dataset (**GSM4572187**) for now.

process/  
**├── GSM4572187_Control1_filtered_peak_bc_matrix.h5  
├── GSM4572187_Control1_fragments.tsv.gz  
├── GSM4572187_Control1_fragments.tsv.gz.tbi**  
├── GSM4572188_Control2_filtered_peak_bc_matrix.h5  
├── GSM4572188_Control2_fragments.tsv.gz  
├── GSM4572188_Control2_fragments.tsv.gz.tbi  
├── GSM4572189_Control3_filtered_peak_bc_matrix.h5  
├── GSM4572189_Control3_fragments.tsv.gz  
├── GSM4572189_Control3_fragments.tsv.gz.tbi  
├── GSM4572190_Control4_filtered_peak_bc_matrix.h5  
├── GSM4572190_Control4_fragments.tsv.gz  
├── GSM4572190_Control4_fragments.tsv.gz.tbi  
├── GSM4572191_Control5_filtered_peak_bc_matrix.h5  
├── GSM4572191_Control5_fragments.tsv.gz  
├── GSM4572191_Control5_fragments.tsv.gz.tbi  
├── GSM4572192_Control1_filtered_feature_bc_matrix.h5  
├── GSM4572193_Control2_filtered_feature_bc_matrix.h5  
├── GSM4572194_Control3_filtered_feature_bc_matrix.h5  
├── GSM4572195_Control4_filtered_feature_bc_matrix.h5  
└── GSM4572196_Control5_filtered_feature_bc_matrix.h5  


## Demo

### Predicting gene activity using Signac  
``signac_and_dataFormat_manipulation_demo.R``

>This R script takes in the peak count matrices in h5 format, infers gene activities via Signac and Seurat, and return 10x MTX format files with proper genome annotations for downstream analysis. Please note that the quality control(QC) steps are not included for simplicity, since the peaks are already filtered. More information for QC can be found in the Seurat official website/tutorials.
### Joint modeling single cell ATAC-seq and RNA-seq using MultiVI
``multiVI_integration.ipynb``
>This notebook takes in the mentioned 10x MTX files, organizes them into a multiome anndata structure and trains a MultiVI model at a small scale.

### Joint modeling using MultiVI and annotating cell type via scANVI
``multiVI_integration_colab.ipynb``
>This notebook is an extension of the multiVI integration analysis above, where the full RNA-seq and ACAC-seq datasets are integrated using google Colab with GPU support. In addition, the cell types are annotated using scANVI's seed cell labeling protocol.

### Joint modeling using GLUE and gene regulatory inference
``GLUE_multimodal_integration_and_scenic_GRN_inference.ipynb``
>This notebook takes in the mentioned scATAC-seq and scRNA-seq datasets and motif information from Jasper database, builds an integrated embedding model and infers gene regulatory network using this model and scenic package.
