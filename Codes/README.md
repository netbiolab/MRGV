# Scripts for MRGV anlaysis
## 0.QualityControl
* HumanDecontamination.py : Removal human reads using bowtie2
* Trimmomatic.py : Trimming adaptors and filter low qualited reads using Trimmomatic
## 1.Assembly
* MEGAHIT.py : Running MEGAHIT for read assembly
* MetaSPAdes.py : Running MetaSPAdes fro read assembly
## 2.ViralContigPrediction
* DeepVirFinder.py : Running DeepVirFinder and filtering confident viral contigs
* Phigaro.py : Running Phigaro to predict Prophage from assemblies
* VIBRANT.py : Running VIBRANT to predict viral contigs and lifestyle
## 3.MergingViralContigs
* Vclust.py : Running VClust for sample-wise deduplication of viral contigs from DeepVirFinder, Phigaro and VIBRANT, using UCLUST
## 4.ContigRevalidation
* GeNomad.py : Running GeNomad on the deduplicated viral contigs for revalidation
* VirRep.py : Running VirRep on the deduplicated viral contigs for revalidation
## 5.ViralBinning
* GetCoverage.py : Computing sample-wise read coverage profile using bowtie2
* GetMetabat2Depth.py : Generating Metabat2 Depth format tables
* GenerateCovTable.py : Generating vRhyme coverage table from Metabat2 Depth table
* MetaBat2.py : Running Metabat2 for viral binning on viral contigs
* Semibin2.py : Running Semibin2 for viral binning on viral contigs
* vRhyme.py : Running vRhyme for viral binning on viral contigs
## 6.BinConsolidation
* BinConsolidate.py : Sample-wise consolidation of bins from Metabat2, Semibin2 and vRhyme
## 7.ProteinAnnotation
* Pharokka.py : Running Pharokka to generate initial annotated GenBank table
* PholdPredict.py : Running Phold Predict to predict 3Di embeddings using FrostT5 model
* PholdCompare.py : Running Phold Compare to find the hits using foldseek
* LinClust.py : Running Linclust in MMSeq2 to generate protein clusters
## 8.ViralReadAlignment
* Minimap2.py : Running minimap to align short reads to viral genomes
* CoverM.py : Running CoverM to calculate alignment coverage
## 9.ClusteringOTUs
* UPGMA.rs : Conduct UPGMA clustering of genomes based on taxonomic rank delineation criteria
## 10.AgingPrediction
* KendallTau.py : Compute Kendall Tau and pvalue, and generating Kendall distance matrix
* Uniqueness.py : Calculate Uniqueness based on distance matrix, with/without cage mates
* Maaslin2.R : Running Masslin2 to extract significantly differential viral taxa
* XGBoostRegressor.py : Running XGBoostRegressor to predict mice ages using viral genus abundance table
