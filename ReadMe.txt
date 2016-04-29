The script "Execution script-parse_blast.R" is used to automatically read the blast results
downloaded from the blast search for hits on the V3kl V6r region.

The execution script and the function script ("function-parse_blast.R") should be saved in
the same place.

Step 1: Set your working directory
The full path to where your R scripts are saved. The path needs to be in quotes. You should
use forwardslashes (/) between folders in your path, or use a double backslash (\\) between folders.

Step 2: sourcing required functions
No action required here. Simply run the code in the execution script as is.

Step 3: Use the parse_blast function:
This function parses the xml output of ncbi blast of V3kl V6r regions
It returns query cover, evalue, percent identity, matched organism(s), 
accession number, and query sequence

It returns hits with the highest total score (as assigned by blast). For example, ngroups=1 gives
all hits with the highest total score, ngroups=2 gives all hits with the highest and second highest
total scores. The total scores correspond with percent identity and percent coverage of the hits.
Therefore, higher total scores correspond with greater percent identity and coverage.

The function takes the following arguments:
filename = name of xml file containing blast search results for V3klV6r sanger 
           sequencing (exported from ncbi blast); full path to file should be included;
	   This file should be obtained from a ncbi blast search, 
	   where you go to Download > XML to download the search results as one XML file.

ngroups= default 1; the number of identity-coverage groups to be returned

tophit = default FALSE; if TRUE, only returns the top hit for each query

output = default 'blast_result.csv'; name of output file; will be in csv format;
          file extension required; will be saved in same location as input file.