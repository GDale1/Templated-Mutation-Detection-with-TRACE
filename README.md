# Templated-Mutation-Detection-with-TRACE

## Description:

TRACE is a novel computational script written and executed in MATLAB (v.2018b) that detects templated mutations in somatically mutated sequences. The script utilizes sequential nested Monte Carlo simulations to test if sets of somatic mutations in an IGHV subsequence match a subsequence in a reference database above an empirically derived background. 

## Overview:

Included in these scripts are several auxiliary scripts that correctly format files for analysis with TRACE. Using the “2_IMGT-gapped-nt-sequences.xlsx” file as an input (see below), the script “TRACE_IMGT_Formatting.m” creates a series of subdirectories with fasta files composed of sequences from the IMGT sequence file that is matched with a reference sequence from the IMGT database. These fasta files must contain >50 sequences and are processed to remove gapped sequences and trim edges. 
Within each subdirectory, the script “AutomatedAnalysisTRACE.m” can be run. The script will perform the TRACE pipeline and return an output file named “FinalCleanedResults.mat” that contains information of the subsequences that underwent a templated mutation event, the template used, and associated information on each. 
Data from each subdirectory can then be loaded into a single file, using the “Load_All_Results.m” script. This outputs a file “AllDataByIgHV.mat” which contains a nested cell structure that contains all the data for each subfolder, grouped by IGHV sequence. 

## How to run TRACE:

Before running the TRACE pipeline, the IMGT output file “2_IMGT-gapped-nt-sequences.txt” must be converted to a .xlsx (Microsoft Excel) file format. 
*Note: the .txt file can be imported to Microsoft Excel and subsequently saved as .xlsx. See the example located in the Test Dataset folder.
Next, use the script “TRACE_IMGT_Formating.m” which has the following input arguments: “ExcelFile”, “ReferenceFasta”, and “Clean_InputHandle”. The “ExcelFile” argument should be the name of the .xlsx formatted IMGT output file. Generally, this is “2_IMGT-gapped-nt-sequences.xlsx”. The “ReferenceFasta” argument should reference the included fasta file “Human_IgHV_Reference_total.fasta” which includes all the human IGHV reference sequences used by IMGT. The “Clean_InputHandle” argument is the file extension for fasta input. Should be made to be '*.fas'.

After the completion of this step, multiple subdirectories will be created, each containing a fasta for a respective IGHV. Within each folder, the script “AutomatedAnalysisTRACE.m” should be run. Input arguments for this script are listed below:

•	InputFileHandle - File extension for fasta input. Usually '*.fas', '*.fa', or '*.fasta'. "*" required by MATLAB and serves as a wildcard. 

•	InRange - the Motif size to be used in the TRACE pipeline. As configured in the paper this value is 38.

•	NumMuts - the number of mutations within a given InRange value to qualify as a mutation cluster. As configured in the paper this value is 8.

•	NumReps - the number of background datasets to generate. As configured in the paper this value is 10.

•	BlastDB - the FASTA file for which a local BLAST database is generated from, and which will be used for the subsequent TRACE analysis. As configured in the paper this was the hg38 Human Reference Genome. 

•	AnalysisName - User defined name for the TRACE run. Will be stored in output .mat file

•	species - can be either 'Human' or 'Mouse'. Controls chromosome number and gene coordinates in the TRACE_PostProcessing step

•	MaxNumMuts - controls the number of bins in the output Histogram figures created in the Mutation analysis steps. As configred in the paper, this value is 30.

•	Path_RequiredFiles - Path to the RequiredFiles Folder

•	Path_BlastDB - Path to the BlastDB folder where the BlastDB fasta file is located

Upon completion of each folder’s TRACE run* the data can be summed by utilizing the “Load_All_Results.m” script which has no input arguments. This outputs a file “AllDataByIgHV.mat” which contains a nested cell structure with the variable name “AllFinalResults” that contains all the data for each subfolder, grouped by IGHV sequence. 

*Note: TRACE is a computationally intensive script and requires a long runtime. On average, each IgHV analysis takes 1 day to 1 week. **Running in parallel is highly recommended**
