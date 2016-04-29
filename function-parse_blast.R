#///////////////////////////////////////////////////////////////////////////////////
# Function to parse xml output of ncbi blast of V3kl V6r regions
# returns query cover, evalue, percent identity, matched organism(s), 
# accession number, query sequence
# returns hits in percent identity-coverage groups, starting from the group 
# with the highest scores in both those measures, and in decreasing order

# input:
# filename = name of xml file containing blast search results for V3klV6r sanger 
#            sequencing (exported from ncbi blast); full path to file should be included
# ngroups = default 1; the number of identity-coverage groups to be returned
# tophit = default FALSE; if TRUE, only returns the top hit for each query
# output = default 'blast_result.csv'; name of output file; will be in csv format;
#          will be saved in same location as input file.

parse_blast <- function(filename=NULL, ngroups=1, tophit=FALSE,
                        output='blast_result.csv') {
  
  # loading required packages
  if(require(xml2)==FALSE) install.packages('XML')
  if(require(stringr)==FALSE) install.packages('stringr')
  if(require(hash)==FALSE) install.packages('hash')
  
  library(xml2)
  library(stringr)
  library(hash)
  
  # make sure folder defined---------------------------------------------------------
  if(is.null(filename)) stop('Name of file not supplied. Please check the ReadMe to make sure it is entered correctly', call.=FALSE)
  
  # making sure folder path is in correct syntax
  path <- str_replace_all(filename, '\\\\', '/')
  
  # Reading file
  xmlfile <- read_xml(filename)
  
  blast_node <- xml_find_all(xmlfile, '//Iteration')
  blast_list <- as_list(blast_node)
  
  out <- c()
  # going through one query (sample) at a time
  for(i in 1:length(blast_list)){

    # sample
    query_sample <- xml_text(xml_find_all(blast_list[[i]], 'Iteration_query-def'))
    query_num <- xml_text(xml_find_all(blast_list[[i]], 'Iteration_iter-num'))
    query_len <- xml_text(xml_find_all(blast_list[[i]], 'Iteration_query-len'))
    
    #print a message to mark progress
    msg <- sprintf('Reading hit results for file "%s"', query_sample)
    message(msg)
    
    # hit
    hit_node <- xml_find_all(blast_list[[i]], 'Iteration_hits/Hit')
    hit_list <- as_list(hit_node)
    
    # create hash map to store keys and values for quick look up
    hmap <- hash()
    key_order <- c()
    
    # initializing hits
    hit <- c()
    
    # going through one hit at a time
    for(j in 1:length(hit_list)) {
      
      # hit values of interest
      fullmatch <- xml_text(xml_find_one(hit_list[[j]], 'Hit_def'))
      
      if(length(unlist(str_extract_all(fullmatch, 'strain')))>1) {
        match <- str_extract_all(fullmatch, '\\w+ [a-z]+(?= strain)')
        match <- unlist(match)
        match <- paste(match, collapse=', ')
      }
      else match <- str_extract(fullmatch, '.*(?= strain)')
      
      access <- xml_text(xml_find_one(hit_list[[j]], 'Hit_accession'))
      hit_num <- xml_text(xml_find_one(hit_list[[j]], 'Hit_num'))
      
      totscore <- xml_text(xml_find_all(hit_list[[j]], 'Hit_hsps/Hsp/Hsp_bit-score'))
      alignlen <- xml_text(xml_find_all(hit_list[[j]], 'Hit_hsps/Hsp/Hsp_align-len'))
      hsp_ident <- xml_text(xml_find_all(hit_list[[j]], 'Hit_hsps/Hsp/Hsp_identity'))
      qfrom <- xml_text(xml_find_all(hit_list[[j]], 'Hit_hsps/Hsp/Hsp_query-from'))
      qto <- xml_text(xml_find_all(hit_list[[j]], 'Hit_hsps/Hsp/Hsp_query-to'))
      evalue <- xml_text(xml_find_all(hit_list[[j]], 'Hit_hsps/Hsp/Hsp_evalue'))
      qseq <- xml_text(xml_find_all(hit_list[[j]], 'Hit_hsps/Hsp/Hsp_qseq'))
      
      # converting values to numeric
      totscore <- as.numeric(totscore)
      alignlen <- as.numeric(alignlen)
      hsp_ident <- as.numeric(hsp_ident)
      qfrom <- as.numeric(qfrom)
      qto <- as.numeric(qto)
      evalue <- as.numeric(evalue)
      
      identity <- round(hsp_ident/alignlen*100,2)
      coverage <- round(alignlen/qto*100,2)
      
      # use identity-coverage combo as hash key
      hkey <-sprintf('q%s_%0.2f%0.2f', query_num, identity, coverage)

      # if hash does not contain key, add key to hash
      if(has.key(hkey, hmap)==FALSE) {
        .set(hmap, hkey, c())
        
        # keep track of key order
        key_order <- c(key_order, hkey)
      }
      
      # if key is in hash, build value entry
      if(has.key(hkey, hmap)) {
        
        # first, retrieve key's value
        hkey_val <- hmap[[hkey]]
        
        # info from hit to be entered into key value
        val_entry <- data.frame(query_num=query_num, sample_name=query_sample,
                                seq_len=query_len, hit_num=hit_num, 
                                match=match, accession=access, total_score=totscore,
                                perc_cover=coverage, perc_ident=identity, 
                                evalue=evalue, match_description=fullmatch,
                                query_seq=qseq)
        
        hkey_val <- rbind(hkey_val, val_entry)
        
        # assign value to key
        hmap[[hkey]] <- hkey_val
      }
      if(tophit==TRUE) break
    } # end of hit loop
    
    # retrieving all hits in the specified number of groups
    retrieve_key <- key_order[1:ngroups]
    
    # converting results to dataframe
    retrieve <- hmap[retrieve_key]
    retrieve <- as.data.frame(as.list(retrieve))
    
    cname <- colnames(retrieve)
    cname <- str_extract_all(cname, '\\w*$', simplify=TRUE)
    cname <- cname[,1]
    colnames(retrieve) <- cname

    out <- rbind(out, retrieve)
    
  } # end of query loop
  
  # export blast results
  if(tophit==TRUE)  message('Only the top hit of each query has been returned')
  out_path <- str_extract(path, '.*(?=/)')
  out_fname <- file.path(out_path, output)
  
  write.csv(file=out_fname, out)
}