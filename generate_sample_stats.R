# Script for generating information on the samples analyzed for the
# morphotype barcoding paper

library(data.table)
data_path <- "~/Documents/Manuscripts/2023_IBA_data_paper/git-uppmax_IBA_data_April2023/data/Sweden/"

sample_field_labels <- samples <- c("SPCEJE","SUKEY2","SVLS5M","S7CRVD","S7536S","SUSTFJ","SGNIVM","SXUEFH","SEKD4C","SMBAHB")
samples_metadata <- read.table(paste0(data_path, "CO1_sequencing_metadata_SE_2019.tsv"), header=TRUE, sep="\t")

count_col_names <- samples_metadata$sampleID_LAB[match(samples,samples_metadata$sampleID_FIELD)]

D <- read.table <- read.table(paste0(data_path, "CO1_cleaned_nochimera_cluster_counts_SE_2019.tsv"), header=TRUE, sep="\t")
T <- read.table <- read.table(paste0(data_path, "CO1_cleaned_nochimera_cluster_taxonomy_SE_2019.tsv"), header=TRUE, sep="\t")
colnames(T)[1] <- "cluster"
spikeins <- read.table(paste0(data_path, "malaise_bio_spikeins_SE_2019.tsv"),header=TRUE,sep=" ")

# Extract the cluster read information
D <- D[ ,c("cluster",count_col_names)]

cluster_reads <- rep(0,times=nrow(D))
for (col in count_col_names) {
    cluster_reads <- cluster_reads + D[,col]
}

D <- D[cluster_reads!=0 & !(D$cluster %in% spikeins$cluster),]
D <- D[T$Class[match(D$cluster,T$cluster)]=="Insecta" | T$Class[match(D$cluster,T$cluster)]=="Collembola",]

asv_counts  <- fread <- read.table(paste0(data_path, "CO1_asv_counts_SE_2019.tsv"), header=TRUE, sep="\t")
asv_taxonomy <- read.table <- read.table(paste0(data_path, "CO1_asv_counts_SE_2019.tsv"), header=TRUE, sep="\t")

# Extract the ASV information
asv_counts <- asv_counts[ ,count_col_names]

asv_reads <- rep(0,times=nrows(asv_counts))
for (col in count_col_names) {
    asv_reads <- asv_reads + asv_counts[,col]
}

asv_counts <- asv_counts[asv_reads!=0 & !(asv_counts$ASV %in% spikeins$ASV),]
asv_counts <- asv_counts[T$Class[match(asv_counts$ASV,T$ASV)]=="Insecta" | T$T$Class[match(asv_counts$ASV,T$ASV)]=="Entognatha",]

write.table(asv_counts,"ASV_counts.tsv")
write.table(D,"cluster_counts.tsv")
