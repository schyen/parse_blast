#////////////////////////////////////////////////////////////////////////////////////
# Use this script to read the blast results of the V3klVr6 sequences 
# This script automatically parses all the results into one csv file.

# To execute this code, fill out the required information according to the 
# instructions in the ReadMe file.
# When ready, highlight the whole script and click "Run"

# Step 1: Set your working directory-------------------------------------------------

path <- '[your/path/here]'
setwd(path)

# Step 2: Sourcing required functions------------------------------------------------

func_path <- file.path(path, 'function-parse_blast.R')
source('function-parse_blast.R')

# Step 3: Use the parse_blast function------------------------------------------------

parse_blast(filename=NULL, ngroups=1, tophit=FALSE, output='blast_result.csv')
