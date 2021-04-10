## Detection and classification of differential translation-efficiency genes with DTEG.R

Citation: Chothani, S., Adami, E., Ouyang, J. F., Viswanathan, S., Hubner, N., Cook, S. A., Schafer, S., & Rackham, O. J. L. (2019). deltaTE: Detection of translationally regulated genes by integrative analysis of Ribo-seq and RNA-seq data. Current Protocols in Molecular Biology, 129, e108. doi.org/10.1002/cpmb.108

Link to publication
https://currentprotocols.onlinelibrary.wiley.com/doi/full/10.1002/cpmb.108


**Preparing input files**  
Calculating differentially-TE genes (DTEGs) requires the count matrices for both Ribo-seq and RNA-seq. These should be the raw counts obtained from feature counts (or any other read counting software), they should not be normalized or batch corrected. Each row should represent a gene and each column represents a sample. The matrix should have a header as shown below.

Ribo-seq count matrix (RPFs): 

 | Gene ID | Sample 1 | Sample 2 | Sample 3 | Sample 4 |
 | --------|----------|----------|----------|----------|
 | Gene 1  | 1290     | 130      | 2	   | 1000     |
 | Gene 2  | 2	     | 10	| 5	   | 1	      |
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | Gene Z  | 200	     | 140	| 15	   | 11	      |


RNA-seq count matrix (mRNA counts): 

 | Gene ID | Sample 5 | Sample 6 | Sample 7 | Sample 8 |
 | --------|----------|----------|----------|----------|
 | Gene 1  | 1290     | 130      | 2	   | 1000     |
 | Gene 2  | 2	     | 10	| 5	   | 1	      |
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | Gene Z  | 200	     | 140	| 15	   | 11	      |


Next, we also need to prepare a also requires a sample information file which should follow the same sample order as count matrices. This file outlines the condition, sequencing type and batch for each sample. This script requires the same header names as shown below (case sensitive). If your experiment has more covariates it is recommended to use the Alternate protocol.

Sample information file:

 | SampleID | Condition | SeqType | Batch |
 | --------|----------|----------|----------|
 | Sample 1  | 1     | RIBO      | 1	   | 
 | Sample 2  | 1     | RIBO      | 2	   | 
 | Sample 3  | 2     | RIBO      | 3	   | 
 | Sample 4  | 2     | RIBO      | 4	   | 
 | Sample 5  | 1     | RNA      | 1	   | 
 | Sample 6  | 1     | RNA      | 2	   | 
 | Sample 7  | 2     | RNA      | 3	   | 
 | Sample 8  | 2     | RNA      | 4	   | 


**Running script DTEG.R for detection and classification of DTGs and DTEGs**  
Once the input files are ready the script DTEG.R can be executed on the bash shell prompt as follows:

$ Rscript --vanilla DTEG.R arg1 arg2 arg3 arg4 arg5 arg6  
where,  
      Argument 1 (arg1): Ribo-seq count matrix  
      Argument 2 (arg2): RNA-seq count matrix  
      Argument 3 (arg3): Sample information file  
      Argument 4 (arg4): Batch effect covariate yes=1, or no=0  
      Argument 5 (arg5): Default = 1, Save Rdata  
      Argument 6 (arg6): Default = 0, Verbose mode  

Example: 
$ Rscript --vanilla DTEG.R ribo_counts.txt rna_counts.txt sample_info.txt 1


**Note: Log fold shrinkage function usage is updated now to adapt to latest DESeq2 version updates**
