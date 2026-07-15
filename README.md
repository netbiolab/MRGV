

<img width="60%" height="60%" alt="Logo" src="https://github.com/user-attachments/assets/661876a0-836d-4320-a82d-5d3118e0b3a5" /> 

# MRGV: Mouse Reference Gut Virome
MRGV provides 109,778 high-confidence viral genomes representing 28,824 species-level vOTUs, together with over 46% of 1.3 million non-redundant viral protein sequences annotated using structure-informed PHROG assignments.

**You can access and browse all MRGV data and information in https://doi.org/10.5281/zenodo.20115989**
#### Citation
>   Kim, H.J. et al (2026). A genomic catalog of the mouse gut virome reveals features associated with ageing _Nature communications_



## Overview of MRGV
#### MRGV consists of 109,778 high-confidence vMAGs, represented by 28,824 species-level vOTUs
<img width="80%" height="80%" alt="Main_page3" src="https://github.com/user-attachments/assets/9a4e73f0-89bc-4585-b975-4952a53fe3cc" />

## MRGV construction pipeline
<img width="5938" height="3297" alt="Nevigation tab_notab" src="https://github.com/user-attachments/assets/f270b077-0e91-4d08-982d-4bcbe02a297f" />

## MRGV Data
### 01. Metadata
| Data |Description| Link |
| --- | --- | --- |
|**MRGV_Representative_Metadata.tsv**|Metadata for 28,824 representative vMAGs|[Click to download (8.5MB)](https://zenodo.org/records/20115989/files/MRGV_Representative_Metadata.tsv?download=1)|
|**MRGV_METADATA_ALL_GENOMES.csv**|Metadata for 109,778 All vMAGs|[Click to download (23.3MB)](https://zenodo.org/records/20115989/files/MRGV_All_Genomes_Metadata.tsv?download=1)|

### 02. MRGV Genomes
| Data |Description| Link |
| --- | --- | --- |
|**MRGV_Representative_Genomes.tar.gz**|28,824 Representative vMAGs|[Click to download (498.9MB)](https://zenodo.org/records/20115989/files/MRGV_Representative_Genomes.tar.gz?download=1)|
|**MRGV_All_Genomes.tar.gz**|109,778 vMAGs of MRGV All vMAGs|[Click to download (1.4GB)](https://zenodo.org/records/20115989/files/MRGV_All_Genomes.tar.gz?download=1)|

### 03. MRGV Protein clusters
| Data |Description| Link |
| --- | --- | --- |
|**MRGV_PC_ID100.tar.gz**|A total of 1,376,499 CDS and metadata, clusterd with 100% AAI|[Click to download (234.3MB)](https://zenodo.org/records/20115989/files/MRGV_PC_ID100.tar.gz?download=1)|
|**MRGV_PC_ID90 DB.tar.gz**|A total of 954,585 CDS and metadata, clusterd with 90% AAI|[Click to download (154.6MB)](https://zenodo.org/records/20115989/files/MRGV_PC_ID90.tar.gz?download=1)|
|**MRGV_PC_ID70 DB.tar.gz**|A total of 746,733 CDS and metadata, clusterd with 70% AAI|[Click to download (121.6MB)](https://zenodo.org/records/20115989/files/MRGV_PC_ID70.tar.gz?download=1)|
|**MRGV_PC_ID50 DB.tar.gz**|A total of 652,176 CDS and metadata, clusterd with 70% AAI|[Click to download (107.0MB)](https://zenodo.org/records/20115989/files/MRGV_PC_ID50.tar.gz?download=1)|
|**MRGV_PC_ID30 DB.tar.gz**|A total of 625,774 CDS and metadata, clusterd with 70% AAI|[Click to download (102.1MB)](https://zenodo.org/records/20115989/files/MRGV_PC_ID30.tar.gz?download=1)|

### 04. Kraken2 DB
| Data |Description| Link |
| --- | --- | --- |
|**MRGV_Repr_Kraken2DB.tar.gz**|Kraken2 DB for 28,824 representative vMAGs|[Click to download (1.5GB)](https://zenodo.org/records/20115989/files/MRGV_Repr_Kraken2DB.tar.gz?download=1)|
|**MRGV_All_Variant_kraken2DB.tar.gz**|Kraken2 DB for 109,778 All vMAGs|[Click to download (3.4GB)](https://zenodo.org/records/20115989/files/MRGV_All_Variant_kraken2DB.tar.gz?download=1)|
|**MRGV_MGBC_Combined_kraken2DB.tar.gz**|Kraken2 DB for 28,824 MRGV Repr vMAG + 26,640 MGBC Strain vMAGs|[Click to download (8.5GB)](https://zenodo.org/records/20115989/files/MGBC_26640_CONCAT_MRGV_REPR_kraken2DB.tar.gz?download=1)|

## Data related to Aginig prediction for MRGV anlaysis (Directory: Data)
### 0.Data
>  * DO_2997_MICE_METADATA.csv : All metadata for 2997 mice metagenomic samples
>  * MRGV_Aging_XGBoostRegressor.yml : Conda enviroment configure to reproduce Mouse Aging clock models
### 1.Aging_predictions
>  * READEME.md : All information for Mouse Aginig clock models and its reproducibility
>  * README_XGBoostRegressor_CageGrouped.md : Description of XGBoostRegressor model training scripts

## Scripts for MRGV anlaysis (Directory: Codes)
### 0.QualityControl
>  * HumanDecontamination.py : Removal human reads using bowtie2
>  * Trimmomatic.py : Trimming adaptors and filter low qualited reads using Trimmomatic
### 1.Assembly
>  * MEGAHIT.py : Running MEGAHIT for read assembly
>  * MetaSPAdes.py : Running MetaSPAdes fro read assembly
### 2.ViralContigPrediction
>  * DeepVirFinder.py : Running DeepVirFinder and filtering confident viral contigs
>  * Phigaro.py : Running Phigaro to predict Prophage from assemblies
>  * VIBRANT.py : Running VIBRANT to predict viral contigs and lifestyle
### 3.MergingViralContigs
>  * Vclust.py : Running VClust for sample-wise deduplication of viral contigs from DeepVirFinder, Phigaro and VIBRANT, using UCLUST
### 4.ContigRevalidation
>  * GeNomad.py : Running GeNomad on the deduplicated viral contigs for revalidation
>  * VirRep.py : Running VirRep on the deduplicated viral contigs for revalidation
### 5.ViralBinning
>  * GetCoverage.py : Computing sample-wise read coverage profile using bowtie2
>  * GetMetabat2Depth.py : Generating Metabat2 Depth format tables
>  * GenerateCovTable.py : Generating vRhyme coverage table from Metabat2 Depth table
>  * MetaBat2.py : Running Metabat2 for viral binning on viral contigs
>  * Semibin2.py : Running Semibin2 for viral binning on viral contigs
>  * vRhyme.py : Running vRhyme for viral binning on viral contigs
### 6.BinConsolidation
>  * BinConsolidate.py : Sample-wise consolidation of bins from Metabat2, Semibin2 and vRhyme
### 7.ProteinAnnotation
>  * Pharokka.py : Running Pharokka to generate initial annotated GenBank table
>  * PholdPredict.py : Running Phold Predict to predict 3Di embeddings using FrostT5 model
>  * PholdCompare.py : Running Phold Compare to find the hits using foldseek
>  * LinClust.py : Running Linclust in MMSeq2 to generate protein clusters
### 8.ViralReadAlignment
>  * Minimap2.py : Running minimap to align short reads to viral genomes
>  * CoverM.py : Running CoverM to calculate alignment coverage
### 9.ClusteringOTUs
>  * UPGMA.rs : Conduct UPGMA clustering of genomes based on taxonomic rank delineation criteria
### 10.AgingPrediction
>  * KendallTau.py : Compute Kendall Tau and pvalue, and generating Kendall distance matrix
>  * Uniqueness.py : Calculate Uniqueness based on distance matrix, with/without cage mates
>  * Maaslin2.R : Running Masslin2 to extract significantly differential viral taxa
>  * XGBoostRegressor.py : Running XGBoostRegressor to predict mice ages using viral genus abundance table
