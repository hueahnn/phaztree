# 02/06/2025 | HA | appending data pulled from uniprot

# load libraries
library(jsonlite)
library(httr)

protein_ids <- filtered_by_score$domain_name # list of just the ids from jackhmmer output table

# appending domains
domain_databases <- c("InterPro", "NCBIfam", "Pfam") # databases to search for
domains_df <- data.frame(uniprot_id = id, domains = I(do)) # dataframe of ids and domains

# iterating through every protein in the jackhmmer file and extracting the id + domains from filtered databases
# for (i in 1:length(protein_ids)) {
for (i in 38885:length(protein_ids)) {
  # extracting id and formatting url to obtain data from json file
  uniprot_id <- sub("UniRef90_", "", protein_ids[i])
  url <- paste0("https://www.uniprot.org/uniprot/", uniprot_id, ".json")
  response <- GET(url)
  # json file for protein
  protein_metadata <- fromJSON(content(response, "text", encoding = "UTF-8"))
  # extracting all the databases cited for the protein
  database <- protein_metadata$uniProtKBCrossReferences$database
  # creating a list to store the domains
  d <- c()
  #iterating through every database extracted to obtain domain info
  for (j in 1:length(database)) {
    if (!is.null(database) && database[j] %in% domain_databases) {
      location <- protein_metadata$uniProtKBCrossReferences
      d <- c(d,location$properties[[j]]$value[1])
    }
  }
  if (i %% 25 == 0) print(i) # sanity check output
  domains_df <- rbind(domains_df, data.frame(unitprot_id = uniprot_id, domains = I(list(d))))
}
result <- result %>% mutate(domains = domains_df$domains)



# appending full AA sequences and domain ids
more_info_df <- tibble()
# for (i in 1:10) {
for (i in 21189:length(protein_ids)) {
  # extracting id and formatting url to obtain data from json file
  uniprot_id <- sub("UniRef90_", "", protein_ids[i])
  url <- paste0("https://www.uniprot.org/uniprot/", uniprot_id, ".json")
  response <- GET(url)
  # json file for protein
  protein_metadata <- fromJSON(content(response, "text", encoding = "UTF-8"))
  
  # extracting all the info cited for the protein
  database <- protein_metadata$uniProtKBCrossReferences$database # databases
  sequence <- protein_metadata$sequence$value
  if (is.null(sequence)) {
    sequence <- "missing"
  }
  d <- tibble() # df to store domains
  
  #iterating through every database extracted to obtain domain info
  for (j in 1:length(database)) {
    if (!is.null(database) && database[j] %in% domain_databases) {
      location <- protein_metadata$uniProtKBCrossReferences
      d <- bind_rows(d, tibble(id = location$id[j], domains = location$properties[[j]]$value[1]))
    }
  }
  
  if (i %% 25 == 0) print(i) # sanity check output
  more_info_df <- bind_rows(more_info_df, tibble(id = uniprot_id, full_AA_seq = sequence, domains = list(d)))
  
}


# obtain missing seqs
protein_ids <- more_info_df$id
for (i in 18425:length(protein_ids)) {
  curr <- more_info_df[i,]
  if(curr$full_AA_seq == "missing") {
    # extract uniParc id from uniprot
    uniParcId <- curr$id # if uniRef id is already the uniParc id
    # extract full length seq from uniParc
    url <- paste0("https://www.uniprot.org/uniparc/", uniParcId, ".json")
    response <- GET(url)
    protein_metadata <- fromJSON(content(response, "text", encoding = "UTF-8"))
    full_seq <- protein_metadata$sequence$value
    more_info_df <- more_info_df %>% mutate(full_AA_seq = replace(full_AA_seq, id == curr$id, full_seq))
  }
  if (i %% 25 == 0) print(i) # sanity check output
}

# appending organism + more specific domain info + gene name + synonym
protein_ids <- tree_seqs$id
organisms <- c()
domain_databases <- c()
gene_names <- c()
gene_synonyms <- c()
for (i in 1:length(protein_ids)) {
  # extracting id and formatting url to obtain data from json file
  url <- paste0("https://www.uniprot.org/uniprot/", protein_ids[i], ".json")
  response <- GET(url)
  # json file for protein
  protein_metadata <- fromJSON(content(response, "text", encoding = "UTF-8"))
  
  # extracting all the info cited for the protein
  organism <- protein_metadata$organism$scientificName
  domain_database <- NULL
  if ("NCBIfam" %in% protein_metadata$uniProtKBCrossReferences$database) {
    domain_database <- protein_metadata$uniProtKBCrossReferences$id[protein_metadata$uniProtKBCrossReferences$database == "NCBIfam"]
  }
  gene_name <- protein_metadata$genes$geneName$value
  gene_synonym <- protein_metadata$genes$synonyms[[1]]$value
  
  append_or_null <- function(lst, val) {
    append(lst, if (is.null(val)) "null" else val)
  }
  
  organisms <- append_or_null(organisms, organism)
  domain_databases <- append_or_null(domain_databases, domain_database)
  gene_names <- append_or_null(gene_names, gene_name)
  gene_synonyms <- append_or_null(gene_synonyms, gene_synonym)
  
  if (i %% 25 == 0) print(i) # sanity check output
}
tree_seqs <- tree_seqs %>% mutate(organism = organisms, gene_name = gene_names, gene_synonym = gene_synonyms)

# fix TIGRfam id pulling
domain_databases2 <- c()
for (i in 5516:length(protein_ids)) {
  # extracting id and formatting url to obtain data from json file
  url <- paste0("https://www.uniprot.org/uniprot/", protein_ids[i], ".json")
  response <- GET(url)
  # json file for protein
  protein_metadata <- fromJSON(content(response, "text", encoding = "UTF-8"))
  
  domain_database <- NULL
  if ("NCBIfam" %in% protein_metadata$uniProtKBCrossReferences$database) {
    domain_database <- protein_metadata$uniProtKBCrossReferences$id[protein_metadata$uniProtKBCrossReferences$database == "NCBIfam"]
  }
  domain_databases2 <- append(domain_databases2, if (is.null(domain_database)) "null" else list(domain_database))
  
  
  if (i %% 25 == 0) print(i) # sanity check output
}
tree_seqs <- tree_seqs %>% mutate(tigrfam = domain_databases2)
