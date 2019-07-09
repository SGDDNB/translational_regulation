Calculating differentially-TE genes (DTEGs) requires the count matrices for both Ribo-seq and RNA-seq. These should be the raw counts obtained from feature counts (or any other read counting software), they should not be normalized or batch corrected. Each row should represent a gene and each column represents a sample. The matrix should have a header as shown below.

Ribo-seq counts (RPFs)

 | Gene ID | Sample 1 | Sample 2 | Sample 3 | Sample 4 |
 | --------|----------|----------|----------|----------|
 | Gene 1  | 1290     | 130      | 2	   | 1000     |
 | Gene 2  | 2	     | 10	| 5	   | 1	      |
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | ..	  |	     | 		|	   |	      |	
 | Gene Z  | 200	     | 140	| 15	   | 11	      |



Next, we also need to prepare a also requires a sample information file which should follow the same sample order as count matrices. This file outlines the condition, sequencing type and batch for each sample. This script requires the same header names as shown below (case sensitive). If your experiment has more covariates it is recommended to use the Alternate protocol.

Running script DTEG.R to detect

Once the input files are ready the script DTEG.R can be executed on bash command line as follows:
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

