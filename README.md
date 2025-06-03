# Investigating sex hormone–immune interactions in reproductive success

# Overview
Sex hormones play an important role in shaping immune responses, which in turn influence reproductive outcomes. This dynamic interaction is crucial for understanding reproductive success, yet the mechanisms underlying it are not fully understood.
This project is aimed to study this interaction between sex hormones, the immune system, and reproductive outcomes, with a focus on understanding the molecular and cellular mechanisms that govern reproductive success.

In this study we utilized **single-cell RNA sequencing analysis** and **proteome analysis** to study immune responses at a systems level. Our analyses aimed to:

- Examine how sex hormones influence immune responses.
- Examine the impact of immune responses on reproductive success.

## Key Methods
- **Single-cell RNA sequencing analysis** followed standard **Seurat** workflow.
- **Sex hormone receptor expression profiling** on PBMCs and menstrual effluent samples.
- **Differential gene expression** and **gene set enrichment analysis (GSEA)**.
- **Immune proteomic profiling** of serum and follicular fluid from IVF patients.

## Software and package versions
- **R version 4.4.2**
- **Seurat version 5.1.0**

# Repository Structure

## `PBMC_analysis/`
This folder contains the code used to analyze PBMCs (Peripheral Blood Mononuclear Cells) using single-cell RNA sequencing.

The dataset used is taken from **Oelen et al. (2022)** :  
> Oelen, R. et al. *Single-cell RNA-sequencing reveals widespread personalized, context-specific gene expression regulation in immune cells*. Nat Commun **13**, 2532 (2022).  
> [Paper Link](https://www.nature.com/articles/s41467-022-30893-5)  
> [Dataset Link](https://eqtlgen.org/sc/datasets/1m-scbloodnl.html)

### Dataset Details

- Captured and sequenced data using **10x Genomics v2 reagents**.
- Includes the following files:
  - `barcodes.tsv.gz`
  - `features.tsv.gz`
  - `matrix.mtx.gz`
- Also includes **cell type annotations** and **stimulation information** and **cell barcode** annotations from the original study.

### Analysis Workflow
- **`0_Cell_type_frequency_analysis.Rmd`**:  
  - The script applies **manual annotation** using known immune cell markers from literature and compares these to the original annotations from Oelen et al.
  - **Seurat** is used for standard quality control and preprocessing.
  - **Sex prediction** is performed using the `cellXY` function to assign each cell a biological sex (male/female).
  - Cell type proportions are calculated by using **propeller test** within the speckle package to compare the frequency of each cell
  - Cell type proportions are **stratified by sex**, and the differences in immune cell composition between males and females are explored and visualized.

- **`1_sex_hormone_receptor_analysis.Rmd`**:  
  - Filters the dataset to include only female samples.
  - Annotates immune cells for expression of key sex hormone receptors (*ESR1, ESR2, PGR, AR, GPER1, FSHR, LHR, GNRHR*). Searches for aliases to capture all receptor variants.
  - Calculates the percentage of cells expressing each receptor across different immune cell types.
  - Visualizes the percentage of receptor positive cells across immune cell types using violon plots.
  
- **`2_DE_analysis.Rmd`**:  
  - Looks at the number of receptor positive vs receptor negative cells across different immune cell types.
  - Performs **differential expression analysis** of all the cell types, positive vs negative, for all receptors, using DESeq2 within Seurat.
  - Visualizes the significant (p.adjust < 0.05) upregulated and downregulated genes for each receptor.
  
- **`3_GSEA_analysis.Rmd`**:  
  - Performs **Gene Set Enrichment Analysis (GSEA)** using the BTM (Blood transcription modules) database.
  - Visualizes the results using dot plots for significantly (p.adjust < 0.05) enriched pathways.

## `Menstrual_effluent_analysis/`

This folder contains the R markdown script that analyzes the datatest downloaded from Shih et al. 2022:
> Shih, T. et al. *Immune and endometrial alterations in women with endometriosis and infertility*. BMC Med **20**, 264 (2022).  
> [Paper Link](https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-022-02500-3)  
> [Dataset Link (GEO Accession: GSE203191)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203191)

- The script performs **quality control** and **preprocessing** of the menstrual effluent dataset, and filters to include control samples only.
- It performs **cell type annotation** using known immune cell markers to compare these to the original annotations from Shih et al.
- It annotates immune cells for expression of key sex hormone receptors (*ESR1, ESR2, PGR, AR, GPER1, FSHR, LHR, GNRHR*). Searches for aliases to capture all receptor variants.
- It calculates the percentage of cells expressing each receptor across different immune cell types.
- It does **differential expression analysis** of all the cell types, positive vs negative, for all receptors, using DESeq2 within Seurat.
- It performs **Gene Set Enrichment Analysis (GSEA)** using the BTM (Blood transcription modules) database.


## `Proteome_analysis/`

This folder contains the script that analyzes the proteomic dataset requested from Nenonen et al. 2024:
> Nenonen, H., Kondic, A., Henic, E., & Hjelmér, I. *Recurrent implantation failure and inflammatory markers in serum and follicle fluid of women undergoing assisted reproduction*. J Affect Disord **340**, 18–25 (2024).  
> [Paper Link](https://www.sciencedirect.com/science/article/pii/S0165037824000184)

- The dataset was provided by the authors upon request, and the script is used to analyze the data and visualize the results.
- The script performs exploratory analysis of the proteomic dataset, that includes serum and follicular fluid measurements of 10 cytokines.
- It visualzes the immune signatures associated with recurrent implantation failure (RIF) and successful pregnancy.
- It performs correlation analysis of the cytokines within the control and RIF groups in serum and follicular fluid.
- It performs logistic regression analysis to identify the factors associated with pregnancy outcomes, including protein concentrations in serum, follicular fluid
  and clinical characteristics like, age, BMI, menstrual cycle length.

Some instructions and comments are also present within all the scripts, along with, session information.
