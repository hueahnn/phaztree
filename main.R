
setwd("C:/Users/stera/Desktop/arkinlab/phaz")
wd <- getwd()

# load libraries
library(rhmmer)
library(ggplot2)
library(dplyr)
library(Biostrings)
library(pointr)
library(tidyverse)
library(conflicted)
library(jsonlite)
library(httr)
library(stringr)

# reading in the jackhmmer output table from a json file
setwd("C:/Users/stera/Desktop/arkinlab/phaz/json files")
result <- read_json("031225_result.json", simplifyVector = TRUE)
result_backup <- data.frame(result)


# iteration vs outputs plot
jack.out <- list.files('jackhmmer_outputs/')
jack.out
rep.size <- NULL
jack.list <- NULL
pb2 = txtProgressBar(min = 0, max = length(jack.out), initial = 0, style = 3) 
for (x in 1:length(jack.out)){
  jack.res <- read_tblout(file = paste0("jackhmmer_outputs/", jack.out[x]))
  rep.size <- append(rep.size, nrow(jack.res))
}
rep.size
iterations <- c(1,2,3,4,5)
plot_data <- data.frame(iterations, rep.size)
iter_vs_size <- ggplot(plot_data, aes(x = iterations, y = rep.size)) + geom_point() + geom_line()
iter_vs_size

# appending AA lengths and s value
result <- result %>% mutate(sequence_svalue = -log2(sequence_evalue)) %>% arrange(domain_name)
AA_String <- readAAStringSet('filtered_sequences.fasta')
AA_String
lengths <- data.frame(names(AA_String), width(AA_String))
lengths <- lengths %>% arrange(names.AA_String.)
result <- result %>% mutate(AA_length = lengths$width.AA_String.)

# color by "cluster"
cluster1 <- result %>% filter(sequence_svalue > 520 & AA_length > 350 & AA_length < 600)
cluster2 <- result %>% filter(sequence_svalue > 450 & sequence_svalue < 520 & AA_length > 350 & AA_length < 600)
cluster3 <- result %>% filter(sequence_svalue > 100 & sequence_svalue < 500 & AA_length > 50 & AA_length < 300)
cluster4 <- result %>% filter(sequence_svalue > 20 & sequence_svalue < 150 & AA_length > 320 & AA_length < 500)
cluster5 <- result %>% filter(sequence_svalue > 20 & sequence_svalue < 80 & AA_length > 550 & AA_length < 650)

# plot of s value vs length
sval_vs_length <- ggplot() + geom_point(data = result, aes(x = sequence_svalue, y = AA_length), alpha = 0.15) + ylim(0,1100)
sval_vs_length
# plot colored by "clusters"
sval_vs_length_clustered <- ggplot() + geom_point(data = result, aes(x = sequence_svalue, y = AA_length), alpha = 0.15) + 
  geom_point(data = cluster1, aes(x = sequence_svalue, y = AA_length), alpha = 0.15, colour = "purple") + 
  geom_point(data = cluster2, aes(x = sequence_svalue, y = AA_length), alpha = 0.15, colour = "green") +
  geom_point(data = cluster3, aes(x = sequence_svalue, y = AA_length), alpha = 0.15, colour = "yellow") +
  geom_point(data = cluster4, aes(x = sequence_svalue, y = AA_length), alpha = 0.15, colour = "magenta") +
  geom_point(data = cluster5, aes(x = sequence_svalue, y = AA_length), alpha = 0.15, colour = "cyan") +
  ylim(0,1100)
sval_vs_length_clustered


# analyzing each cluster by randomly sampling 10 proteins from each cluster
cluster_df <- data.frame()
for (i in 1:5) {
  str = paste("cluster",i, sep = "")
  cluster = ptr("cluster", str)
  samples <- sample.int(length(cluster$domain_name), 10)
  ids = c()
  d = c()
  for (j in 1:10) {
    ids <- c(ids, cluster$domain_name[j])
    if (is.null(cluster$domains[j])) {
      d <- c(d, "n/a")
    }
    else {
      d <- c(d, cluster$domains[j])
    }
  }
  cluster_df <- rbind(cluster_df, data.frame(unitprot_id = I(list(ids)), domains = I(list(d))))
}

# plot by domain
# filter: "PHB_depoly_PhaZ"
phaz <- result %>% filter(map_lgl(domains, ~"PHB_depoly_PhaZ" %in% .))
phaz_plot <- ggplot() + geom_point(data = result, aes(x = sequence_svalue, y = AA_length), alpha = 0.2) +
  geom_point(data = phaz, aes(x = sequence_svalue, y = AA_length), alpha = 0.2, colour = "violet") + ylim(0,1100)
phaz_plot
phaz_only_plot <- ggplot() + 
  geom_point(data = phaz, aes(x = sequence_svalue, y = AA_length), alpha = 0.2) + ylim(0,1100)
phaz_only_plot


# filter: "PHB_depo_C"
phac <- result %>% filter(map_lgl(domains, ~"PHB_depo_C" %in% .))
phac_plot <- ggplot() + geom_point(data = result, aes(x = sequence_svalue, y = AA_length), alpha = 0.2) +
  geom_point(data = phac, aes(x = sequence_svalue, y = AA_length), alpha = 0.15, colour = "palegreen") + ylim(0,1100)
phac_plot
phac_only_plot <- ggplot() +
  geom_point(data = phac, aes(x = sequence_svalue, y = AA_length), alpha = 0.2) + ylim(0,1100)
phac_only_plot

# filter: "PHA/PHB_synthase"
pha_phb_synthase <- result %>% filter(map_lgl(domains, ~"PHA/PHB_synthase" %in% .))
pha_phb_synthase_plot <- ggplot() + geom_point(data = result, aes(x = sequence_svalue, y = AA_length), alpha = 0.2) +
  geom_point(data = pha_phb_synthase, aes(x = sequence_svalue, y = AA_length), alpha = 0.15, colour = "steelblue") + ylim(0,1100)
pha_phb_synthase_plot

# combine the 3 above
combined_plot <- ggplot() + geom_point(data = result, aes(x = sequence_svalue, y = AA_length), alpha = 0.2) +
  geom_point(data = phaz, aes(x = sequence_svalue, y = AA_length), alpha = 0.15, colour = "red") +
  geom_point(data = phac, aes(x = sequence_svalue, y = AA_length), alpha = 0.1, colour = "yellow") + 
  geom_point(data = pha_phb_synthase, aes(x = sequence_svalue, y = AA_length), alpha = 0.05, colour = "steelblue") + ylim(0,1100)
combined_plot

# filter: "AB_hydrolase_fold"
ab_hydrolase_fold <- result %>% filter(map_lgl(domains, ~"AB_hydrolase_fold" %in% .))
ab_hydrolase_fold_plot <- ggplot() + geom_point(data = result, aes(x = sequence_svalue, y = AA_length), alpha = 0.2) +
  geom_point(data = ab_hydrolase_fold, aes(x = sequence_svalue, y = AA_length), alpha = 0.15, colour = "mediumorchid") + ylim(0,1100)
ab_hydrolase_fold_plot

# filter: "TIGR01849"
tigr_domain <- result %>% filter (map_lgl(domains, ~"TIGR01849" %in% .))
tigr_domain

# filter out sequences
filtered_by_score <- result %>% filter(sequence_svalue > 25 & AA_length < 900 & AA_length > 150)
filtered_plot <- ggplot() + geom_point(data = filtered_by_score, aes(x = sequence_svalue, y = AA_length), alpha = 0.15) + ylim(0,1100)
filtered_plot


# how many full length seqs are missing
more_info_df <- fromJSON("031225_more_info_df.json")
missing_full_seq <- more_info_df %>% filter(full_AA_seq == "missing")
missing_full_seq

# filtering out unique seqs and outputting a fasta file
full_seqs_only <- more_info_df %>% filter(full_AA_seq != "missing")
length(unique(full_seqs_only$full_AA_seq)) # don't need this, all seqs are unique
PhaZ.blastp.uni <- full_seqs_only[!duplicated(full_seqs_only$full_AA_seq),,]
PhaZ.blastp.fasta.unique <- AAStringSet(full_seqs_only$full_AA_seq)
names(PhaZ.blastp.fasta.unique) <- paste(full_seqs_only$id)

out <- tempfile()
writeXStringSet(PhaZ.blastp.fasta.unique, out)
file.copy(out, "031325_phaZ_blastP.fa")
file_path <- "031325_phaZ_blastP.fa"
file.access(file_path, mode = 2) 
file.access(getwd(), mode = 2)


setwd("C:/Users/stera/Downloads")
file.copy(out, "031325_phaZ_blastP.fa")


# making mapping file for tree by domain (04/15/25)
file.create("041525_tree_mapping.txt")
result <- data.frame(result_backup)
result <- result %>% mutate(domain_name = str_replace(domain_name, "UniRef90_", ""))
# restructuring and appending data for phaz domain
phaz_full_seq <- full_seqs_only %>% filter(id %in% phaz$domain_name)
phaz_map <- phaz %>% mutate(name = phaz_full_seq$id)
phaz_map <- phaz_map %>% select(name)
color <- rep("pthc_red", 3522)
phaz_map <- phaz_map %>% mutate(branch_color = color)
write.table(phaz_map, "042225_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# phaz, phac, pha/phb synthase
phac_full_seq <- full_seqs_only %>% filter(id %in% phac$domain_name)
phac_map <- phac_full_seq %>% mutate(name = phac_full_seq$id) %>% select(name)
color <- rep("pthc_yellow", 3694)
phac_map <- phac_map %>% mutate(branch_color = color)
write.table(phac_map, "042225_phac_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)

synthase_full_seq <- full_seqs_only %>% filter(id %in% pha_phb_synthase$domain_name)
synthase_map <- synthase_full_seq %>% mutate(name = synthase_full_seq$id) %>% select(name)
color <- rep("pthc_blue", 14535)
synthase_map <- synthase_map %>% mutate(branch_color = color)
write.table(synthase_map, "042225_synthase_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)

combo <- rbind(phaz_map, phac_map, synthase_map)
write.table(combo, "042225_combined_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# making mapping files for tree by database (4/22/25)
phaz1 <- result %>% filter(domain_name == "Q0KCI0")
phaz1 <- data.frame(name = "Q0KCI0", branch_color = "pthc_yellow")
write.table(phaz1, "042725_phaz1_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)
tigrfam = c()
for (i in (1:3522)) {
  curr = phaz_full_seq$domains[[i]]
  if (any(curr$id %in% c("TIGR01849"))) {
    tigrfam = c(tigrfam,phaz_full_seq$id[i])
    next
  }
}
tigr_phaz <- data.frame(name = tigrfam, branch_color = rep("ptv_orange", 3508))
write.table(tigr_phaz, "042825_tigrphaz_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)
interpro = c()
for (i in (1:3522)) {
  curr = phaz_full_seq$domains[[i]]
  if (any(curr$id %in% c("IPR010915"))) {
    interpro = c(interpro,phaz_full_seq$id[i])
    next
  }
}
ipr_phaz <- data.frame(name = interpro, branch_color = rep("ptv_teal", 3522))
write.table(ipr_phaz, "042825_iprphaz_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)


# appending all info for tree into one table
result <- result %>% mutate(id = sub("UniRef90_", "", domain_name))
tree_seqs <- result %>% filter(id %in% full_seqs_only$id)
# filter only cupriavidus necator related sequences
cupriavidus <- tree_seqs %>% filter(str_detect(organism, "Cupriavidus necator"))
cup_mapping <- data.frame(name = cupriavidus$id,leaf_label_color = rep("ptv_pink",19))
write.table(cup_mapping, "052925_cupriavidus_necator_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)

TIGR01838 <- tree_seqs %>% filter(tigrfam == "TIGR01838") 
map_01838 <- data.frame(name = TIGR01838$id, leaf_label_color = rep("ptv_teal", 3327))
write.table(map_01838, "053025_01838_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)

TIGR01839 <- tree_seqs %>% filter(tigrfam == "TIGR01839") 
map_01839 <- data.frame(name = TIGR01839$id, leaf_label_color = rep("ptb_purple", 124))
write.table(map_01839, "053025_01839_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)

TIGR01836 <- tree_seqs %>% filter(tigrfam == "TIGR01836") 
map_01836 <- data.frame(name = TIGR01836$id, leaf_label_color = rep("ptb_green", 1051))
write.table(map_01836, "053025_01836_tree_mapping.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# exporting data frames to json files and importing json files back into R in a readable format
write_json(tree_seqs, "060225_tree_seqs.json", pretty = TRUE)
write_json(full_seqs_only, "060225_full_seqs_only.json", pretty = TRUE)
t <- read_json("060225_tree_seqs.json", simplifyVector = TRUE)
