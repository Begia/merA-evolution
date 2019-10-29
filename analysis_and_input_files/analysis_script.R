if (!requireNamespace("BiocManager")) {
  install.packages("BiocManager")
}

#devtools::install_github("lme4/lme4")

packages <- c("ggplot2", "plyr", "scales", "svglite", "ggthemes", "tidyverse", "biomformat", "phyloseq", 
              "limma", "reshape2", "vegan", "ape", "metagenomeSeq", "Rtsne", "dbscan", "ranger", "edarf",
              "caret", "seqinr", "gplots", "ggtree", "tidytree", "car", "lme4", "sp", "maptools","cshapes",
              "maps", "rgdal", "rgeos", "Rcpp", "segmented")

is.installed <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, ask=F, version= "3.8")
  }
  sapply(pkg, require, character.only = TRUE)
}
is.installed(packages)

CSS_otu_table <- function(data) {
data.metagenomeSeq <- newMRexperiment(t(data), featureData=NULL, libSize=NULL, normFactors=NULL)
p <- cumNormStat(data.metagenomeSeq) #default is 0.5
data.cumnorm <- cumNorm(data.metagenomeSeq, p=p)
data.CSS <-  t(MRcounts(data.cumnorm, norm=TRUE, log=TRUE)) 
return(data.CSS)
}


theme_set(theme_tufte(base_family = "sans", base_size = 18) + theme(panel.border = element_rect(colour = "black", fill = NA), 
                                                                    axis.text = element_text(colour = "black", size = 18)))

#veganifyOTU function extracted from internal phyloseq function
veganifyOTU <- function(physeq){
  if(taxa_are_rows(physeq)){physeq <- t(physeq)}
  return(as(otu_table(physeq), "matrix"))
}

#qc the conserved residues Tyr605/Phe605 (Phe confers lower activity) for all protein sequences
merA_prot_qc <- read.fasta("E:/merA_2018/merA_proteins_for_qc.fasta")
merA_prot_qc <- do.call(rbind, merA_prot_qc)
rownames(merA_prot_qc) <- sub("\\(translated\\)", "", rownames(merA_prot_qc))
merA_prot_qc <- merA_prot_qc[-grep("ENA", rownames(merA_prot_qc)),]
to_remove_OTUs <- rownames(merA_prot_qc)[-(which(merA_prot_qc[,95] %in% c("y","f")))]
write.table(to_remove_OTUs, file = "E:/merA_2018/merA_proteins_to_remove.list", col.names = F, row.names = F, quote = F)

#subset the alignments for both rpoB and merA
#first rpoB
alignment <- read.alignment("E:/rpoB_2018/rpoB_trimmed_translatorx_mafft.nt_ali.fasta", format="fasta")

#calculate the codons at each position in frame +1
alignment_codons <- list(NULL)
for (i in 1:alignment$nb) {
  alignment_codons[[i]] <- splitseq(s2c(alignment$seq[[i]]), frame=0, word=3)
}
alignment_codons <- do.call(rbind, alignment_codons)

#calculate consensus codons (most common) at each position
consensus_codons <- apply(alignment_codons, 2, function(x) names(which.max(table(x))))

#check which positions have a gap as the most common codon, but limit the tails of sequences from the examined region
gap_columns <- which(grepl("---", consensus_codons[1:343]))

#find sequences that are bridging these "gap positions"
bad_sequences <- list(NULL)
i <- 1
for (co in gap_columns) {
  bad_sequences[[i]] <- alignment$nam[which(!grepl("---", alignment_codons[,co]))]
  i <- i+1
}

#write out the bridging sequences for gaps which are bridged by < 10% of all sequences 
bridging_sequences <- unique(unlist(bad_sequences[which(lengths(bad_sequences) <= 0.1*alignment$nb)]))
write.table(bridging_sequences, "E:/rpoB_2018/rpoB_bridging_sequences.list", row.names=F, col.names=F, quote=F)

# following code run at CAC cluster to realign with mafft and remove outliers with treeshrink
# faSomeRecords ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed_translatorx_mafft.nt_ali.fasta ~/matti/rpoB_2018/rpoB_all/rpoB_bridging_sequences.list ~/matti/rpoB_2018/rpoB_all/rpoB_bridging_sequences.fasta
# faSomeRecords -exclude ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed_translatorx_mafft.nt_ali.fasta ~/matti/rpoB_2018/rpoB_all/rpoB_bridging_sequences.list ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed_refined_translatorx_mafft.nt_ali.fasta
# sed -i 's/-//g' ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed_refined_translatorx_mafft.nt_ali.fasta
# perl ~/.local/bin/translatorx_vLocal.pl -i ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed_refined_translatorx_mafft.nt_ali.fasta -o ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed_refined_v2_translatorx_mafft -c 11 -p F
# ~/.local/bin/trimal -in ~/matti/rpoB_2018/rpoB_all/rpoB_trimmed_refined_v2_translatorx_mafft.nt_ali.fasta -out ~/matti/rpoB_2018/rpoB_all/rpoB_trimal.fasta -gappyout
# awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ~/matti/rpoB_2018/rpoB_all/rpoB_trimal.fasta | grep "-" -B1 | grep ">" | sed 's/>//g' > ~/matti/rpoB_2018/rpoB_all/rpoB_indel_seq.list
# faSomeRecords -exclude ~/matti/rpoB_2018/rpoB_all/rpoB_trimal.fasta ~/matti/rpoB_2018/rpoB_all/rpoB_indel_seq.list ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_noindel.fasta
# sed -i 's/\s.\+$//g' ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_noindel.fasta
# ~/.local/bin/FastTreeMP -gtr -gamma -nt ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_noindel.fasta > ~/matti/rpoB_2018/rpoB_all/rpoB_trimal.tre
# treeshrink.py -q 0.01 -i ~/matti/rpoB_2018/rpoB_all/rpoB_trimal.tre
# tr -s "\t" "\n" < ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_treeshrink/rpoB_trimal_shrunk_RS_0.01.txt > ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_treeshrink/rpoB_outliers.list
# faSomeRecords -exclude ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_noindel.fasta ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_treeshrink/rpoB_outliers.list ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_v2.fasta
# ~/.local/bin/FastTreeMP -gtr -gamma -nt ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_v2.fasta > ~/matti/rpoB_2018/rpoB_all/rpoB_trimal_v2.tre

#then merA
alignment <- read.alignment("E:/merA_2018/merA_translatorx_mafft.nt_ali.fasta", format="fasta")

#calculate the codons at each position in frame +1
alignment_codons <- list(NULL)
for (i in 1:alignment$nb) {
  alignment_codons[[i]] <- splitseq(s2c(alignment$seq[[i]]), frame=0, word=3)
}
alignment_codons <- do.call(rbind, alignment_codons)

#calculate consensus codons (most common) at each position
consensus_codons <- apply(alignment_codons, 2, function(x) names(which.max(table(x))))

#check which positions have a gap as the most common codon, but limit the tails of sequences from the examined region
gap_columns <- which(grepl("---", consensus_codons))

#find sequences that are bridging these "gap positions"
bad_sequences <- list(NULL)
i <- 1
for (co in gap_columns) {
  bad_sequences[[i]] <- alignment$nam[which(!grepl("---", alignment_codons[,co]))]
  i <- i+1
}

#write out the bridging sequences for gaps which are bridged by < 10% of all sequences 
bridging_sequences <- unique(unlist(bad_sequences[which(lengths(bad_sequences) <= 0.1*alignment$nb)]))
write.table(bridging_sequences, "E:/merA_2018/merA_bridging_sequences.list", row.names=F, col.names=F, quote=F)

# following code run at CAC cluster to realign with mafft and remove outliers with treeshrink
# faSomeRecords ~/matti/merA_2018/merA_all/merA_translatorx_mafft.nt_ali.fasta ~/matti/merA_2018/merA_all/merA_bridging_sequences.list ~/matti/merA_2018/merA_all/merA_bridging_sequences.fasta
# faSomeRecords -exclude ~/matti/merA_2018/merA_all/merA_translatorx_mafft.nt_ali.fasta ~/matti/merA_2018/merA_all/merA_bridging_sequences.list ~/matti/merA_2018/merA_all/merA_refined_translatorx_mafft.nt_ali.fasta
# sed -i 's/-//g' ~/matti/merA_2018/merA_all/merA_refined_translatorx_mafft.nt_ali.fasta
# perl ~/.local/bin/translatorx_vLocal.pl -i ~/matti/merA_2018/merA_all/merA_refined_translatorx_mafft.nt_ali.fasta -o ~/matti/merA_2018/merA_all/merA_refined_v2_translatorx_mafft -c 11 -p F
#  ~/.local/bin/trimal -in ~/matti/merA_2018/merA_all/merA_refined_v2_translatorx_mafft.nt_ali.fasta -out ~/matti/merA_2018/merA_all/merA_trimal.fasta -gappyout
# awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' ~/matti/merA_2018/merA_all/merA_trimal.fasta | grep "-" -B1 | grep ">" | sed 's/>//g' > ~/matti/merA_2018/merA_all/merA_indel_seq.list 
# faSomeRecords -exclude ~/matti/merA_2018/merA_all/merA_trimal.fasta ~/matti/merA_2018/merA_all/merA_indel_seq.list ~/matti/merA_2018/merA_all/merA_trimal_noindel.fasta
# sed -i 's/\s.\+$//g' ~/matti/merA_2018/merA_all/merA_trimal_noindel.fasta
# ~/.local/bin/FastTreeMP -gtr -gamma -nt ~/matti/merA_2018/merA_all/merA_trimal_noindel.fasta > ~/matti/merA_2018/merA_all/merA_trimal.tre
# treeshrink.py -q 0.01 -i ~/matti/merA_2018/merA_all/merA_trimal.tre
# tr -s "\t" "\n" < ~/matti/merA_2018/merA_all/merA_trimal_treeshrink/merA_trimal_shrunk_RS_0.01.txt > ~/matti/merA_2018/merA_all/merA_trimal_treeshrink/merA_outliers.list
# faSomeRecords -exclude ~/matti/merA_2018/merA_all/merA_trimal_noindel.fasta ~/matti/merA_2018/merA_all/merA_trimal_treeshrink/merA_outliers.list ~/matti/merA_2018/merA_all/merA_trimal_v2.fasta
# ~/.local/bin/FastTreeMP -gtr -gamma -nt ~/matti/merA_2018/merA_all/merA_trimal_v2.fasta > ~/matti/merA_2018/merA_all/merA_trimal_v2.tre

merA_biom <- import_biom("E:/merA_2018/swarm_otu_table_no_singletons.biom")
merA_mapping <- import_qiime_sample_data("E:/merA_2018/mapping_corrected_concatenated.txt")
merA_fasttree <- read_tree("E:/merA_2018/merA_trimal_v2.tre")
#test_mapping <- import_qiime_sample_data(data.frame(rbind(a_mapping, b_mapping)))

#since the tree is constricted to only valid merA sequences, the phyloseq object is automatically subset to them
merA_phylo <- phyloseq(merA_biom, sample_data=merA_mapping, merA_fasttree)
sample_names(merA_phylo) <- c("AQT01", "AQT03", "AQT04", "Extraction", "HAR01", "HAR02", "HAR41", "KEV01", "KEV10", "KEV08", "LAH01", 
"LAH11", "LAH07", "OJA01", "OJA15", "OJA21", "OJA45", "PAI01", "PAI04", "PAI09", "POK01", "POK17", "POK21", "POK29", "POK38", "POK09", "PUL01", "PUL10", "PUL06", "VUO01", "VUO05", "VUO07")

#import tot-Hg data
hgdata <- read.csv("E:/merA_2018/finland_sediment_totHg_overview.csv", sep = "\t")
hgdata$Site <- substring(hgdata$NAME,1,3)
hgdata$Country <- ifelse(hgdata$Site %in% c("POK", "LAH", "AQT"), "Canada", "Finland")
hgdata$NAME <- gsub("-", "", hgdata$NAME)
hgdata$NAME[which(grepl("[A-Z]{3}[0-9]{1}$", hgdata$NAME))] <- paste0(substr(hgdata$NAME[which(grepl("[A-Z]{3}[0-9]{1}$", hgdata$NAME))], 1, 3), 0, substr(hgdata$NAME[which(grepl("[A-Z]{3}[0-9]{1}$", hgdata$NAME))], 4, 4))
hgdata$Site <- revalue(hgdata$Site, c("HAR" = "Kokemäenjoki / Harjavalta", "OJA" = "Öjanjärvi", "PAI" = "Päiväjärvi", 
                                      "PUL" = "Pulmankijärvi", "VUO" = "Vuolimus Cieskuljavri", "KEV" = "Kevojärvi", "POK" = "Pocket Lake", "LAH" = "Lake Hazen", "AQT" = "Aquatuk Lake"))
lake_levels <- c("Kokemäenjoki / Harjavalta", "Pocket Lake", "Öjanjärvi", "Päiväjärvi", "Aquatuk Lake", "Pulmankijärvi", "Kevojärvi", "Lake Hazen", "Vuolimus Cieskuljavri")
hgdata$Site <- factor(hgdata$Site, lake_levels)
hgdata$PCR.BAND <- factor(hgdata$PCR.BAND, levels=c("Yes", "No", "Not_sampled"))

#improve the sample mapping with tot-Hg, and lake data
new_mapping <- data.frame(Sample=sample_names(merA_phylo), Site=substring(sample_names(merA_phylo),1,3))
new_mapping$Site <- revalue(new_mapping$Site, c("HAR" = "Kokemäenjoki / Harjavalta", "OJA" = "Öjanjärvi", "PAI" = "Päiväjärvi", 
                                                "PUL" = "Pulmankijärvi", "VUO" = "Vuolimus Cieskuljavri", "KEV" = "Kevojärvi", "POK" = "Pocket Lake", "LAH" = "Lake Hazen", "AQT" = "Aquatuk Lake"))
new_mapping$Site <- factor(new_mapping$Site, lake_levels)

#extrapolate and intrapolate the dating data for all the samples
dating_data <- read.csv("E:/merA_2018/Core_dating_raw_data.csv", sep = ",")
hgdata2 <- hgdata[c("NAME", "Site", "CONC..ug.kg.", "MIDDLE.DEPTH", "Country")]
colnames(hgdata2) <- c("Sample", "Site", "Hg.ug.kg", "MIDDLE.DEPTH", "Country")

new_dating_list <- list(NULL)
extrapolation_models_list <- list(NULL)
for (sitevar in levels(dating_data$Site)) {
  indexvar <- which(levels(dating_data$Site) %in% sitevar)
  #don't analyze sites where all available samples have been dated
  if (sitevar %in% c("Kevojärvi", "Kokemäenjoki / Harjavalta")) {
    new_data <- hgdata2[which(hgdata2$Site %in% sitevar),]
    new_dating_data <- dating_data[which(dating_data$Site %in% sitevar),]
    new_data <- merge(new_data, new_dating_data, by="MIDDLE.DEPTH", all.x=T)
    new_data <- new_data[,!grepl("Site.y", colnames(new_data))]
    colnames(new_data) <- gsub("\\.x|\\.y", "", colnames(new_data))
  } else {
      new_dating_data <- dating_data[which(dating_data$Site %in% sitevar),]
      #last_points <- new_dating_data[c(nrow(new_dating_data)-1, c(nrow(new_dating_data))),]
      date_depth_model <- lm(DATE ~ poly(MIDDLE.DEPTH, degree=2), new_dating_data)
      extrapolation_models_list[[indexvar]] <- summary(date_depth_model)
      new_data <- hgdata2[which(hgdata2$Site %in% sitevar),]
      not_dated <- new_data[which(new_data$MIDDLE.DEPTH %in% setdiff(new_data$MIDDLE.DEPTH, new_dating_data$MIDDLE.DEPTH)),]
      if (nrow(not_dated) > 0) {
        not_dated[,"DATE"] <- NA
        interpolation_dates <- approx(new_dating_data$MIDDLE.DEPTH, new_dating_data$DATE, xout=not_dated$MIDDLE.DEPTH[not_dated$MIDDLE.DEPTH < max(new_dating_data$MIDDLE.DEPTH)])$y
        extrapolation_dates <- predict(date_depth_model, newdata = data.frame(MIDDLE.DEPTH=not_dated$MIDDLE.DEPTH[not_dated$MIDDLE.DEPTH > max(new_dating_data$MIDDLE.DEPTH)]))
        not_dated$DATE[not_dated$MIDDLE.DEPTH < max(new_dating_data$MIDDLE.DEPTH)] <- interpolation_dates
        not_dated$DATE[not_dated$MIDDLE.DEPTH > max(new_dating_data$MIDDLE.DEPTH)] <- extrapolation_dates
        added_data <- merge(new_dating_data, not_dated, all=T)[c("MIDDLE.DEPTH", "DATE", "CRS.ERROR")]
        new_data <- merge(new_data, added_data, by=c("MIDDLE.DEPTH"), all.x=T) 
      } else {
          new_data <- dating_data[which(dating_data$Site %in% sitevar),]
      }
      if (sitevar %in% c("Pocket Lake", "Lake Hazen")) {
        new_data$DATE <- new_data$DATE + 2
      } 
  }
  new_dating_list[[indexvar]] <- new_data
}

new_mapping <- do.call(rbind, new_dating_list)
names(extrapolation_models_list) <- levels(dating_data$Site)

new_mapping$Hg.ug.kg <- log10(new_mapping$Hg.ug.kg)
rownames(new_mapping) <- new_mapping$Sample
new_mapping <- new_mapping[,!grepl("Sample", colnames(new_mapping))]

new_hgdata <- hgdata[,!grepl("CONC..ug.kg.", colnames(new_data))]
new_hgdata <- merge(new_hgdata, new_mapping, by=c("Site", "Country", "MIDDLE.DEPTH"))

#if contamination needs to be looked at, here's some code for a bar plot
qc_merA <- data.frame(otu_table(merA_phylo))
qc_merA_cont <- rownames(qc_merA)[which(qc_merA$Extraction > 0 )]
qc_merA <- rownames(qc_merA)[which(qc_merA$Extraction == 0 )]

merA_phylo_norm <- transform_sample_counts(merA_phylo, function(x) 100 * x/sum(x))
merA_phylo_cont <- prune_taxa(qc_merA_cont, merA_phylo_norm)
bars <- data.frame(t(otu_table(merA_phylo_cont)))
sample_code_levels <- read.csv("E:/merA_2018/sample_coding_levels.csv", sep="\t", header=F)
bars <- bars[which(rownames(bars) %in% sample_code_levels$V1),]
rownames(bars) <- sample_code_levels$V2[match(rownames(bars), sample_code_levels$V1)]
bars <- reshape2::melt(as.matrix(bars), varnames=c("Sample", "Variant"), value.name="Abundance")


bars$Variant <- factor(bars$Variant)
bars$Sample <- factor(bars$Sample, levels=rev(sample_code_levels$V2))

image <- ggplot(bars, aes(x=factor(Sample), y=Abundance)) + geom_bar(stat="identity") +
  theme(legend.position="none", axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(),
        plot.title = element_text(size=28, hjust = 0.5), axis.ticks.y = element_blank(), axis.title.y = element_blank()) +
  ylab("Abundance (%)") + ggtitle("Contaminating merA sequences") + ylab("Abundance (% of total reads after QC)") + coord_flip()
ggsave(file = "E:/merA_2018/results/merA_contamination.svg", plot=image, units="mm", width=250, height=250)

qc_merA_phylo_norm <- prune_taxa(qc_merA, merA_phylo)
#remove the test samples "Extraction" and "HAR2" to make datasets 100% comparable
qc_merA_phylo_norm <- subset_samples(qc_merA_phylo_norm, !(Description %in% c("Extraction", "HAR2")))
qc_merA_phylo_norm <- prune_taxa(taxa_sums(qc_merA_phylo_norm) > 0, qc_merA_phylo_norm)
sample_data(qc_merA_phylo_norm) <- new_mapping

rpoB_biom <- import_biom("E:/rpoB_2018/swarm_otu_table_no_singletons.biom")
rpoB_fasttree <- read_tree("E:/rpoB_2018/rpoB_trimal_v2.tre")
rpoB_phylo <- phyloseq(rpoB_biom, rpoB_fasttree)
sample_names(rpoB_phylo) <- c("AQT01", "AQT03", "AQT04", "HAR01", "HAR41", "KEV01", "KEV10", "KEV08", "LAH01", 
                              "LAH11", "LAH07", "OJA01", "OJA15", "OJA21", "OJA45", "PAI01", "PAI04", "PAI09", "POK01", "POK17", "POK21", "POK29", "POK38", "POK09", "PUL01", "PUL10", "PUL06", "VUO01", "VUO05", "VUO07")
qc_rpoB_phylo_norm <- merge_phyloseq(rpoB_phylo, sample_data(new_mapping))

#calculate alpha-diversity for rpoB
dataset_names <- c("merA", "rpoB")

  avgrichness <- list(NULL)
  for (j in 1:10) {
    set.seed(42*j)
    rarefied <- rarefy_even_depth(eval(parse(text = paste0("qc_",dataset_names[i],"_phylo_norm"))))
    richness1 <- estimate_richness(rarefied)
    Goods.cover <- NULL
    for (k in 1:nsamples(rarefied)) {
    Goods.cover[k] <- 1-(length(which(otu_table(rarefied)[,k] == 1))/sample_sums(rarefied)[k])
    }
    richness1 <- cbind(richness1, Goods.cover=Goods.cover)
    richness1$sample <- sample_names(eval(parse(text = paste0("qc_",dataset_names[i],"_phylo_norm"))))
    rownames(richness1) <- seq(1+nrow(richness1)*(j-1),nrow(richness1)*j)
    richness1$set <- dataset_names[i]
    avgrichness[[j]] <- richness1
  }

richness <- do.call("rbind",avgrichness)
richness <- cbind(richness, Site=gsub(".*\\.(.*\\..*)\\.[0-9]", "\\1", richness$sample))
richness <- richness[which(colnames(richness) %in% c("Chao1", "Shannon", "InvSimpson", "Goods.cover", "sample", "Site"))]
richness <- plyr::rename(richness, replace = c("InvSimpson" = "Simpson's dominance", "Goods.cover" = "Good's coverage"))
richness <- melt(richness)


ggplot(richness, aes(x=factor(sample), y=value)) + geom_boxplot() + 
  facet_wrap(~variable, scales = "free_y") + theme(axis.text.x = element_text(angle = 90, vjust = 0.25, hjust = 1, size = 8)) + xlab("Sample") + ylab("Value")


#plot number of reads at different time points
merA_time_reads <- cbind(n.reads=sample_sums(qc_merA_phylo_norm),sample_data(qc_merA_phylo_norm)[,-c(1,3,6)])
merA_time_reads$Site <- factor(merA_time_reads$Site, levels=lake_levels)
merA_time_reads <- merA_time_reads[-which(merA_time_reads$Site %in% "Kokemäenjoki / Harjavalta"),]

lake_colors <- c("#bb6100", #Pocket Lake
                 "#9ca500", #Öjanjärvi
                 "#194dad", #Päiväjärvi
                 "#ff67a4", #Aquatuk Lake
                 "#a9b3ff", #Pulmankijärvi
                 "#01a981", #Kevojärvi
                 "#0084dc", #Lake Hazen
                 "#aad191") #Vuolimus Cieskuljavri

image <- ggplot(merA_time_reads, aes(x=DATE, y=n.reads, colour=Site, shape=Country)) + geom_point(size=6) + 
  theme(legend.position = "none", axis.text.y = element_text(hjust = 0, size = 22), plot.title = element_text(size=28, hjust = 0.5)) + 
  scale_color_manual(values=lake_colors, name = "Site") + ggtitle("Number of merA reads over time") +
  scale_shape_manual(values=c(18,20), name = "Country") + xlab("Date (CE)") + ylab("Number of reads per sample") + 
  scale_y_log10(limits=c(0.9,150000), breaks=c(1,10,100,1000,10000,100000))
ggsave(file = "E:/merA_2018/results/merA_qc_reads.svg", plot=image, units="mm", width=250, height=200)

rpoB_time_reads <- cbind(n.reads=sample_sums(qc_rpoB_phylo_norm),sample_data(qc_rpoB_phylo_norm)[,-c(1,3,6)])
rpoB_time_reads$Site <- factor(rpoB_time_reads$Site, levels=lake_levels)
rpoB_time_reads <- rpoB_time_reads[-which(rpoB_time_reads$Site %in% "Kokemäenjoki / Harjavalta"),]

image <- ggplot(rpoB_time_reads, aes(x=DATE, y=n.reads, colour=Site, shape=Country)) + geom_point(size=6) + 
  theme(legend.position = "none", axis.text.y = element_text(hjust = 0, size = 22), plot.title = element_text(size=28, hjust = 0.5)) + 
  scale_color_manual(values=lake_colors, name = "Site") + ggtitle("Number of rpoB reads over time") +
  scale_shape_manual(values=c(18,20), name = "Country") + xlab("Date (CE)") + ylab("Number of reads per sample") + 
  scale_y_log10(limits=c(0.9,100000), breaks=c(1,10,100,1000,10000,100000))
ggsave(file = "E:/merA_2018/results/rpoB_qc_reads.svg", plot=image, units="mm", width=250, height=200)


#normalize otu tables and rename the variants by sample (concatenated), also change names for fastas
OTU_mapping_to_new_names <- list(NULL)
name_mapping <- list(NULL)
for (i in 1:2){
  if (i == 1) {
    new_otu_table <- otu_table(CSS_otu_table(data.frame(otu_table(qc_merA_phylo_norm))),taxa_are_rows=T)
    new_name_fasta <- read.fasta("E:/merA_2018/merA_trimal_v2.fasta", set.attributes = F, as.string = T)
    #this code used for NCBI submissions
    new_name_fasta2 <- read.fasta("E:/merA_2018/merA_swarm_cluster_seeds.fasta", set.attributes = F, as.string = T)
  } else {
    new_otu_table <- otu_table(CSS_otu_table(data.frame(otu_table(qc_rpoB_phylo_norm))),taxa_are_rows=T)
    new_name_fasta <- read.fasta("E:/rpoB_2018/rpoB_trimal_v2.fasta", set.attributes = F, as.string = T)
    #this code used for NCBI submissions
    new_name_fasta2 <- read.fasta("E:/rpoB_2018/rpoB_swarm_cluster_seeds.fasta", set.attributes = F, as.string = T)
  }
  new_sequence_names <- list(NULL)
  for (j in 1:ncol(new_otu_table)) {
    variant_abundances <- as.vector(new_otu_table[,j])
    new_sequence_names[[j]] <- ifelse((variant_abundances > 0), colnames(new_otu_table)[j] , "")
  }
  new_sequence_names <- do.call(cbind,new_sequence_names)
  new_sequence_names <- apply(new_sequence_names,1, paste, collapse="_")
  new_sequence_names <- gsub("_+","_",new_sequence_names)
  new_sequence_names <- gsub("^_|_$","",new_sequence_names)
  n_categories <- table(new_sequence_names)
  for (j in 1:length(n_categories)) {
    new_sequence_names[which(new_sequence_names %in% names(n_categories)[j])] <- paste(names(n_categories)[j],seq.int(length(new_sequence_names[which(new_sequence_names %in% names(n_categories)[j])])),sep="_")
  }
  rownames(new_otu_table) <- new_sequence_names
   if (i == 1) {
     old_sequence_names <- taxa_names(qc_merA_phylo_norm)
     taxa_names(qc_merA_phylo_norm) <- new_sequence_names
     otu_table(qc_merA_phylo_norm) <- new_otu_table
  } else {
    old_sequence_names <- taxa_names(qc_rpoB_phylo_norm)
    taxa_names(qc_rpoB_phylo_norm) <- new_sequence_names
    otu_table(qc_rpoB_phylo_norm) <- new_otu_table
  }
  name_mapping_1 <- data.frame(old=old_sequence_names, new=new_sequence_names)
  new_name_fasta <- new_name_fasta[names(new_name_fasta) %in% name_mapping_1$old]
  name_mapping[[i]] <- name_mapping_1[order(match(name_mapping_1$old,names(new_name_fasta))),]
  names(new_name_fasta) <- name_mapping_1$new
  OTU_mapping_to_new_names[[i]] <- new_name_fasta
}

#write the new renamed fasta files and trees
write.fasta(sequences=OTU_mapping_to_new_names[[1]], names=names(OTU_mapping_to_new_names[[1]]), file.out="E:/merA_2018/merA_final_renamed.fasta", open="w", nbchar=300, as.string=F)
write.fasta(sequences=OTU_mapping_to_new_names[[2]], names=names(OTU_mapping_to_new_names[[2]]), file.out="E:/rpoB_2018/rpoB_final_renamed.fasta", open="w", nbchar=300, as.string=F)
write.tree(phy=phy_tree(qc_merA_phylo_norm), file="E:/merA_2018/merA_final_renamed.tre")
write.tree(phy=phy_tree(qc_rpoB_phylo_norm), file="E:/rpoB_2018/rpoB_final_renamed.tre")

#write NCBI sequences
merA_ncbi = read.fasta("E:/merA_2018/merA_swarm_cluster_seeds.fasta")
merA_ncbi = data.frame(Fragments=names(merA_ncbi), Seqs=unlist(getSequence(merA_ncbi, as.string=T)))
merA_ncbi$Seqs = gsub('-', '', merA_ncbi$Seqs)
merA_ncbi = merge(merA_ncbi, name_mapping[[1]], by.x='Fragments', by.y='old')
merA_ncbi$Fragments = paste0(">merA_", row.names(merA_ncbi), "- ", "[organism=freshwater sediment metagenome] ", merA_ncbi$new, " amplicon variant from from eDNA.")
merA_ncbi$new = NULL
merA_ncbi <- do.call(rbind, lapply(seq(nrow(merA_ncbi)), function(i) t(merA_ncbi[i, ])))
write.table(merA_ncbi, file="E:/merA_2018/SI files/merA_ncbi.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)

rpoB_ncbi = read.fasta("E:/rpoB_2018/rpoB_swarm_cluster_seeds.fasta")
rpoB_ncbi = data.frame(Fragments=names(rpoB_ncbi), Seqs=unlist(getSequence(rpoB_ncbi, as.string=T)))
rpoB_ncbi$Seqs = gsub('-', '', rpoB_ncbi$Seqs)
rpoB_ncbi = merge(rpoB_ncbi, name_mapping[[2]], by.x='Fragments', by.y='old')
rpoB_ncbi$Fragments = paste0(">rpoB_", row.names(rpoB_ncbi), "- ", "[organism=freshwater sediment metagenome] ", rpoB_ncbi$new, " amplicon variant from from eDNA.")
rpoB_ncbi$new = NULL
rpoB_ncbi <- do.call(rbind, lapply(seq(nrow(rpoB_ncbi)), function(i) t(rpoB_ncbi[i, ])))
write.table(rpoB_ncbi, file="E:/merA_2018/SI files/rpoB_ncbi.fasta", row.names = FALSE, col.names = FALSE, quote = FALSE)

#scale the OTUs to 100% for the bar graphs
qc_bars_merA <- transform_sample_counts(qc_merA_phylo_norm, function(x) 100 * x/sum(x))
qc_bars_rpoB <- transform_sample_counts(qc_rpoB_phylo_norm, function(x) 100 * x/sum(x))

qc_bars_merA <- data.frame(t(otu_table(qc_bars_merA)))
qc_bars_merA <- reshape2::melt(as.matrix(qc_bars_merA), varnames=c("Sample", "Variant"), value.name="Abundance")
qc_bars_merA$Variant <- factor(qc_bars_merA$Variant)

ggplot(qc_bars_merA, aes(x=factor(Sample), y=Abundance, fill=factor(Variant)))+ geom_bar(stat="identity") + 
  theme(legend.position="none",axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + ggtitle("Variants of merA") + ylab("Abundance (%)") + coord_flip()

qc_bars_rpoB <- data.frame(t(otu_table(qc_bars_rpoB)))
qc_bars_rpoB <- reshape2::melt(as.matrix(qc_bars_rpoB), varnames=c("Sample", "Variant"), value.name="Abundance")
qc_bars_rpoB$Variant <- factor(qc_bars_rpoB$Variant)

ggplot(qc_bars_rpoB, aes(x=factor(Sample), y=Abundance, fill=factor(Variant)))+ geom_bar(stat="identity") + 
  theme(legend.position="none",axis.text.y = element_text(hjust = 0, size = 22), panel.border = element_blank(), plot.title = element_text(size=28, hjust = 0.5)) +
  xlab("Sample") + ylab("Abundance (%)") + ggtitle("Variants of rpoB") + ylab("Abundance (%)") + coord_flip()

#ordinations
qc_merA_phylo_norm2 <- subset_samples(qc_merA_phylo_norm, Site != "Kokemäenjoki / Harjavalta")
#this line used for sensitivity testing for low read number samples
#qc_merA_phylo_norm2 <- subset_samples(qc_merA_phylo_norm2, !(sample_names(qc_merA_phylo_norm2) %in% lown_samples))
qc_merA_phylo_norm2 <- prune_taxa(names(taxa_sums(qc_merA_phylo_norm2)[taxa_sums(qc_merA_phylo_norm2)>0]), qc_merA_phylo_norm2)
merA_distance <- DPCoA(qc_merA_phylo_norm2)
set.seed(42)
ordu_merA <- ordinate(merA_phylo, "NMDS", distance = merA_distance$RaoDis)

qc_rpoB_phylo_norm2 <- subset_samples(qc_rpoB_phylo_norm, Site != "Kokemäenjoki / Harjavalta")
qc_rpoB_phylo_norm2 <- prune_taxa(names(taxa_sums(qc_rpoB_phylo_norm2)[taxa_sums(qc_rpoB_phylo_norm2)>0]), qc_rpoB_phylo_norm2)
rpoB_distance <- DPCoA(qc_rpoB_phylo_norm2)
set.seed(42)
ordu_rpoB <- ordinate(rpoB_phylo, "NMDS", distance = rpoB_distance$RaoDis)

continuous_sample_data <- data.frame(sample_data(qc_merA_phylo_norm2))[-c(2, 4, 6)]
model_sample_data <-  data.frame(sample_data(qc_merA_phylo_norm2))[,-6]

datasets <- c(rep("merA",ncol(model_sample_data)),rep("rpoB",ncol(model_sample_data)))
variables <- rep(colnames(model_sample_data), 2)
datasets_ord <- c(rep("ordu_merA",ncol(model_sample_data)),rep("ordu_rpoB",ncol(model_sample_data)))
modelstorun_frame <- data.frame(set=datasets, model=variables, dataset_ord=datasets_ord)

#modelstorun_frame <- modelstorun_frame[-c(5:8),]

env_distance <- vegdist(continuous_sample_data, method="euclidean")
merA_mantel <- mantel(merA_distance$RaoDis, env_distance, method="pearson", permutations=10000)
ef_merA <- envfit(ordu_merA,model_sample_data,permu=10000)

rpoB_mantel <- mantel(rpoB_distance$RaoDis, env_distance, method="pearson", permutations=10000)
ef_rpoB <- envfit(ordu_rpoB,model_sample_data,permu=10000)

ordisurfs <- list(NULL)
pvalues <- c(NULL)
for (j in 1:nrow(modelstorun_frame)) {
  if (modelstorun_frame$model[j] %in% colnames(continuous_sample_data)) {
    ordisurfs[[j]] <- ordisurf(eval(parse(text = paste0(modelstorun_frame$dataset_ord[j],"~",modelstorun_frame$model[j]))), data = data.frame(eval(parse(text = paste0("sample_data(qc_",modelstorun_frame$set[j],"_phylo_norm2)")))),
                               method = "REML", select = TRUE, penalty= 1.4, plot = F, scaling = 3, w =NULL, permu = 10000)
    pvalues[j] <- summary(ordisurfs[[j]])$p.table[1,4]
  } else {
    ordisurfs[[j]] <- NULL
    if (modelstorun_frame$set[j] %in% "merA") {
      pvalues[j] <- unname(ef_merA$factors$pvals[names(ef_merA$factors$pvals) %in% modelstorun_frame$model[j]])
    } else {
      pvalues[j] <- unname(ef_rpoB$factors$pvals[names(ef_rpoB$factors$pvals) %in% modelstorun_frame$model[j]])
    }
  }
}

#format and save ordinations with fitted category centroids
NMDS_merA_Site <- data.frame(NMDS1 = ordu_merA$points[,1], NMDS2 = ordu_merA$points[,2], group=model_sample_data$Site)
ord_merA_Site <- ordiellipse(ordu_merA, model_sample_data$Site, display = "sites", 
                               kind = "se", conf = 0.95, label = T, draw = "none")

df_ell_merA_Site <- data.frame()
for(g in levels(NMDS_merA_Site$group)){
  df_ell_merA_Site <- rbind(df_ell_merA_Site, cbind(as.data.frame(with(NMDS_merA_Site[NMDS_merA_Site$group==g,],
                                                                           vegan:::veganCovEllipse(ord_merA_Site[[g]]$cov,ord_merA_Site[[g]]$center,ord_merA_Site[[g]]$scale))),group=g))
}

ordi_merA_Hg <- ordisurfs[[3]] #fetch the ordisurf object
ordi_grid_merA_Hg <- ordi_merA_Hg$grid #extracts the ordisurf object
ordi_a_merA_Hg <- expand.grid(x = ordi_grid_merA_Hg$x, y = ordi_grid_merA_Hg$y) #get x and ys
ordi_a_merA_Hg$z <- rescale(as.vector(ordi_grid_merA_Hg$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_merA_Hg <- data.frame(na.omit(ordi_a_merA_Hg)) #gets rid of the nas
ordi_a_na_merA_Hg$category <- "Hg"

ordi_merA_Depth <- ordisurfs[[1]] #fetch the ordisurf object
ordi_grid_merA_Depth <- ordi_merA_Depth$grid #extracts the ordisurf object
ordi_a_merA_Depth <- expand.grid(x = ordi_grid_merA_Depth$x, y = ordi_grid_merA_Depth$y) #get x and ys
ordi_a_merA_Depth$z <- rescale(as.vector(ordi_grid_merA_Depth$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_merA_Depth <- data.frame(na.omit(ordi_a_merA_Depth)) #gets rid of the nas
ordi_a_na_merA_Depth$category <- "Depth"

ordi_merA_Year <- ordisurfs[[5]] #fetch the ordisurf object
ordi_grid_merA_Year <- ordi_merA_Year$grid #extracts the ordisurf object
ordi_a_merA_Year <- expand.grid(x = ordi_grid_merA_Year$x, y = ordi_grid_merA_Year$y) #get x and ys
ordi_a_merA_Year$z <- rescale(as.vector(ordi_grid_merA_Year$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_merA_Year <- data.frame(na.omit(ordi_a_merA_Year)) #gets rid of the nas
ordi_a_na_merA_Year$category <- "Year"

ordisurf_NMDS_data_merA <- model_sample_data #there are other ways of doing this. But this is the way I do it for ease of plotting
ordisurf_NMDS_data_merA$NMDS1 <- ordu_merA$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
ordisurf_NMDS_data_merA$NMDS2 <- ordu_merA$points[ ,2]


NMDS_rpoB_Site <- data.frame(NMDS1 = ordu_rpoB$points[,1], NMDS2 = ordu_rpoB$points[,2], group=model_sample_data$Site)
ord_rpoB_Site <- ordiellipse(ordu_rpoB, model_sample_data$Site, display = "sites", 
                             kind = "se", conf = 0.95, label = T, draw = "none")

df_ell_rpoB_Site <- data.frame()
for(g in levels(NMDS_rpoB_Site$group)){
  df_ell_rpoB_Site <- rbind(df_ell_rpoB_Site, cbind(as.data.frame(with(NMDS_rpoB_Site[NMDS_rpoB_Site$group==g,],
                                                                       vegan:::veganCovEllipse(ord_rpoB_Site[[g]]$cov,ord_rpoB_Site[[g]]$center,ord_rpoB_Site[[g]]$scale))),group=g))
}

ordi_rpoB_Hg <- ordisurfs[[8]] #fetch the ordisurf object
ordi_grid_rpoB_Hg <- ordi_rpoB_Hg$grid #extracts the ordisurf object
ordi_a_rpoB_Hg <- expand.grid(x = ordi_grid_rpoB_Hg$x, y = ordi_grid_rpoB_Hg$y) #get x and ys
ordi_a_rpoB_Hg$z <- rescale(as.vector(ordi_grid_rpoB_Hg$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_rpoB_Hg <- data.frame(na.omit(ordi_a_rpoB_Hg)) #gets rid of the nas
ordi_a_na_rpoB_Hg$category <- "Hg"

ordi_rpoB_Depth <- ordisurfs[[6]] #fetch the ordisurf object
ordi_grid_rpoB_Depth <- ordi_rpoB_Depth$grid #extracts the ordisurf object
ordi_a_rpoB_Depth <- expand.grid(x = ordi_grid_rpoB_Depth$x, y = ordi_grid_rpoB_Depth$y) #get x and ys
ordi_a_rpoB_Depth$z <- rescale(as.vector(ordi_grid_rpoB_Depth$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_rpoB_Depth <- data.frame(na.omit(ordi_a_rpoB_Depth)) #gets rid of the nas
ordi_a_na_rpoB_Depth$category <- "Depth"

ordi_rpoB_Year <- ordisurfs[[10]] #fetch the ordisurf object
ordi_grid_rpoB_Year <- ordi_rpoB_Year$grid #extracts the ordisurf object
ordi_a_rpoB_Year <- expand.grid(x = ordi_grid_rpoB_Year$x, y = ordi_grid_rpoB_Year$y) #get x and ys
ordi_a_rpoB_Year$z <- rescale(as.vector(ordi_grid_rpoB_Year$z), to=c(0,100)) #unravel the matrix for the z scores
ordi_a_na_rpoB_Year <- data.frame(na.omit(ordi_a_rpoB_Year)) #gets rid of the nas
ordi_a_na_rpoB_Year$category <- "Year"

ordisurf_NMDS_data_rpoB <- model_sample_data #there are other ways of doing this. But this is the way I do it for ease of plotting
ordisurf_NMDS_data_rpoB$NMDS1 <- ordu_rpoB$points[ ,1] #this puts the NMDS scores for the plots into a new dataframe. you could put them into an existing one if you preferred.
ordisurf_NMDS_data_rpoB$NMDS2 <- ordu_rpoB$points[ ,2]


contours <- list(c("#ff67a4", #Aquatuk Lake
                   "#000000", #Hg
                   "#01a981", #Kevojärvi
                   #"#f82e55", #Kokemäenjoki
                   "#0084dc", #Lake Hazen
                   "#9ca500", #Öjanjärvi
                   "#194dad", #Päiväjärvi
                   "#bb6100", #Pocket Lake
                   "#a9b3ff", #Pulmankijärvi
                   "#aad191"), #Vuolimus Cieskuljavri
                 c("#ff67a4", #Aquatuk Lake
                   "#000000", #Depth
                   "#01a981", #Kevojärvi
                   #"#f82e55", #Kokemäenjoki
                   "#0084dc", #Lake Hazen
                   "#9ca500", #Öjanjärvi
                   "#194dad", #Päiväjärvi
                   "#bb6100", #Pocket Lake
                   "#a9b3ff", #Pulmankijärvi
                   "#aad191"), #Vuolimus Cieskuljavri
                 c("#ff67a4", #Aquatuk Lake
                   "#01a981", #Kevojärvi
                   #"#f82e55", #Kokemäenjoki
                   "#0084dc", #Lake Hazen
                   "#9ca500", #Öjanjärvi
                   "#194dad", #Päiväjärvi
                   "#bb6100", #Pocket Lake
                   "#a9b3ff", #Pulmankijärvi
                   "#aad191", #Vuolimus Cieskuljavri
                   "#000000"), #Year
                 c("#ff67a4", #Aquatuk Lake
                   "#000000", #Hg
                   "#01a981", #Kevojärvi
                   #"#f82e55", #Kokemäenjoki
                   "#0084dc", #Lake Hazen
                   "#9ca500", #Öjanjärvi
                   "#194dad", #Päiväjärvi
                   "#bb6100", #Pocket Lake
                   "#a9b3ff", #Pulmankijärvi
                   "#aad191"), #Vuolimus Cieskuljavri
                 c("#ff67a4", #Aquatuk Lake
                   "#000000", #Depth
                   "#01a981", #Kevojärvi
                   #"#f82e55", #Kokemäenjoki
                   "#0084dc", #Lake Hazen
                   "#9ca500", #Öjanjärvi
                   "#194dad", #Päiväjärvi
                   "#bb6100", #Pocket Lake
                   "#a9b3ff", #Pulmankijärvi
                   "#aad191"), #Vuolimus Cieskuljavri
                 c("#ff67a4", #Aquatuk Lake
                   "#01a981", #Kevojärvi
                   #"#f82e55", #Kokemäenjoki
                   "#0084dc", #Lake Hazen
                   "#9ca500", #Öjanjärvi
                   "#194dad", #Päiväjärvi
                   "#bb6100", #Pocket Lake
                   "#a9b3ff", #Pulmankijärvi
                   "#aad191", #Vuolimus Cieskuljavri
                   "#000000") #Year
)

plots_merA <- c("ordi_a_na_merA_Hg", "ordi_a_na_merA_Depth", "ordi_a_na_merA_Year", "ordi_a_na_rpoB_Hg", "ordi_a_na_rpoB_Depth", "ordi_a_na_rpoB_Year")
ggplot_title <-  c("Phylogenetic diversity of merA over [THg]", "Phylogenetic diversity of merA over sediment depth", "Phylogenetic diversity of merA over dating", 
                   "Phylogenetic diversity of rpoB over [THg]", "Phylogenetic diversity of rpoB over sediment depth", "Phylogenetic diversity of rpoB over dating")

for (m in 1:6) {
  if (m < 4) {
    ordisurf_data <- ordisurf_NMDS_data_merA
    ellipse_data <- df_ell_merA_Site
  } else {
    ordisurf_data <- ordisurf_NMDS_data_rpoB
    ellipse_data <- df_ell_rpoB_Site
  }
image <- ggplot(data = ordisurf_data, aes(NMDS1, NMDS2)) + 
    geom_point(aes(colour = Site, shape = Country), size=6) +
    #geom_text(aes(label=round(ordisurf_data$Hg.ug.kg,2), colour=Site), size= 5) +
    geom_path(data=ellipse_data, aes(x=NMDS1, y=NMDS2, colour=group), size=0.5, linetype=3) +
    geom_contour(data = eval(parse(text=plots_merA[m])), aes(x = x, y = y, z = z, colour=category, alpha=..level..),bins = 20, size=1) + 
    theme(legend.position="none", axis.text.y = element_text(hjust = 0, size = 22), plot.title = element_text(size=28, hjust = 0.5)) + 
    scale_color_manual(values=contours[[m]], name = "Site") + ggtitle(ggplot_title[m]) +
    scale_shape_manual(values=c(18,20), name = "Country")
ggsave(file = paste0("E:/merA_2018/results/", plots_merA[m], "2.svg"), plot=image,
       units="mm", width=250, height=200)
}

#format sequences for evolutionary analyses
merA_sequences <- read.fasta("E:/merA_2018/merA_final_renamed.fasta", set.attributes = F, as.string = T, forceDNAtolower = F)
rpoB_sequences <- read.fasta("E:/rpoB_2018/rpoB_final_renamed.fasta", set.attributes = F, as.string = T, forceDNAtolower = F)
sites <- unique(model_sample_data$Site)

#first for individual cores for merA and rpoB
all_merA_OTUs <- list(NULL)
merA_site_phyloseqs <- list(NULL)
for (i in 1:length(sites)) {
  subset_phylo <- subset_samples(qc_merA_phylo_norm, Site %in% sites[i])
  subset_phylo <- prune_taxa(names(taxa_sums(subset_phylo)[taxa_sums(subset_phylo)>0]), subset_phylo)
  subset_phylo2 <- data.frame(otu_table(subset_phylo))
  site_merA_OTUs <- list(NULL)
  subset_sequence_names <- list(NULL)
  for (j in 1:ncol(subset_phylo2)) {
    subsetlist <- rownames(subset_phylo2)[which(rowSums(subset_phylo2) == subset_phylo2[,j])]
    subset_sequences <- merA_sequences[c(which(names(merA_sequences) %in% subsetlist))]
    site_merA_OTUs[[j]] <- subset_sequences
    names(site_merA_OTUs[[j]]) <- paste(substr(colnames(subset_phylo2)[1], 1,3), 1:length(subsetlist), round(new_mapping$DATE[which(rownames(new_mapping) %in% colnames(subset_phylo2[j]))], 0), sep="_")
    taxa_names(subset_phylo)[which(taxa_names(subset_phylo) %in% names(subset_sequences))] <- names(site_merA_OTUs[[j]])
  }
  site_sequences <- do.call(rbind, lapply(site_merA_OTUs, function (x) do.call(rbind, x)))
  write.fasta(as.list(site_sequences), names=rownames(site_sequences), paste0("E:/merA_2018/BEAST/merA_", substr(colnames(subset_phylo2)[1], 1,3), ".fasta"), open="w")
  all_merA_OTUs[[i]] <- data.frame(Sample=colnames(subset_phylo2), unique=lengths(site_merA_OTUs))
  subset_phylo <- prune_taxa(!grepl("[A-Z]{3}[0-9]{2}", taxa_names(subset_phylo)), subset_phylo)
  merA_site_phyloseqs[[i]] <- subset_phylo
}
all_merA_OTUs <- do.call(rbind, all_merA_OTUs)
all_merA_OTUs <- merge(all_merA_OTUs, new_mapping, by.x = "Sample", by.y = "row.names")


all_rpoB_OTUs <- list(NULL)
rpoB_site_phyloseqs <- list(NULL)
for (i in 1:length(sites)) {
  subset_phylo <- subset_samples(qc_rpoB_phylo_norm, Site %in% sites[i])
  subset_phylo <- prune_taxa(names(taxa_sums(subset_phylo)[taxa_sums(subset_phylo)>0]), subset_phylo)
  subset_phylo2 <- data.frame(otu_table(subset_phylo))
  site_rpoB_OTUs <- list(NULL)
  subset_sequence_names <- list(NULL)
  for (j in 1:ncol(subset_phylo2)) {
    subsetlist <- rownames(subset_phylo2)[which(rowSums(subset_phylo2) == subset_phylo2[,j])]
    subset_sequences <- rpoB_sequences[c(which(names(rpoB_sequences) %in% subsetlist))]
    site_rpoB_OTUs[[j]] <- subset_sequences
    names(site_rpoB_OTUs[[j]]) <- paste(substr(colnames(subset_phylo2)[1], 1,3), 1:length(subsetlist), round(new_mapping$DATE[which(rownames(new_mapping) %in% colnames(subset_phylo2[j]))], 0), sep="_")
    taxa_names(subset_phylo)[which(taxa_names(subset_phylo) %in% names(subset_sequences))] <- names(site_rpoB_OTUs[[j]])
  }
  site_sequences <- do.call(rbind, lapply(site_rpoB_OTUs, function (x) do.call(rbind, x)))
  write.fasta(as.list(site_sequences), names=rownames(site_sequences), paste0("E:/rpoB_2018/BEAST/rpoB_", substr(colnames(subset_phylo2)[1], 1,3), ".fasta"), open="w")
  all_rpoB_OTUs[[i]] <- data.frame(Sample=colnames(subset_phylo2), unique=lengths(site_rpoB_OTUs))
  subset_phylo <- prune_taxa(!grepl("[A-Z]{3}[0-9]{2}", taxa_names(subset_phylo)), subset_phylo)
  rpoB_site_phyloseqs[[i]] <- subset_phylo
}
all_rpoB_OTUs <- do.call(rbind, all_rpoB_OTUs)
all_rpoB_OTUs <- merge(all_rpoB_OTUs, new_mapping, by.x = "Sample", by.y = "row.names")

image <- ggplot(all_merA_OTUs, aes(x=factor(MIDDLE.DEPTH), y=unique)) + geom_bar(stat="identity") + facet_wrap("Site", scales="free") + ylab("Unique sequences among horizons") + ggtitle("N merA variants") +
  theme_tufte(base_family = "sans", base_size = 14) + theme(panel.border = element_rect(colour = "black", fill = NA), axis.text = element_text(colour = "black", size = 12), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(file = "E:/merA_2018/results/merA_variants.svg", plot=image,
       units="mm", width=400, height=200)

image <- ggplot(all_rpoB_OTUs, aes(x=factor(MIDDLE.DEPTH), y=unique)) + geom_bar(stat="identity") + facet_wrap("Site", scales="free") + ylab("Unique sequences among horizons") + ggtitle("N rpoB variants") +
  theme_tufte(base_family = "sans", base_size = 14) + theme(panel.border = element_rect(colour = "black", fill = NA), axis.text = element_text(colour = "black", size = 12), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5))
ggsave(file = "E:/merA_2018/results/rpoB_variants.svg", plot=image,
       units="mm", width=400, height=200)


#heatmap of variants

hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
merA_heatmap_otu_table <- as.matrix(data.frame(otu_table(qc_merA_phylo_norm)))
colnames(merA_heatmap_otu_table) <- sample_names(qc_merA_phylo_norm)
cairo_pdf(filename="E:/merA_2018/results/merA_heatmap_otu_table.pdf", width=12, height=10)
heatmap.2(merA_heatmap_otu_table, col=hmcol, labRow = F)
dev.off()

rpoB_heatmap_otu_table <- as.matrix(data.frame(otu_table(qc_rpoB_phylo_norm)))
colnames(rpoB_heatmap_otu_table) <- sample_names(qc_rpoB_phylo_norm)
cairo_pdf(filename="E:/merA_2018/results/rpoB_heatmap_otu_table.pdf", width=12, height=10)
heatmap.2(rpoB_heatmap_otu_table, col=hmcol, labRow = F)
dev.off()

#subset the sequences per horizon per core to 50 for those that have more than that
#this is based on clustering of patristic distances, and only the first seq of each cluster is preserved
#the sequences to remove from the fastas are written to a file (to be removed from the fastas on CSC)
merA_subset_trees <- readLines("E:/merA_2018/BEAST/merA_subset_trees.list")
merA_to_remove <- list(NULL)
for (i in 1:length(merA_subset_trees)) {
  subset_tree <- read.tree(paste0("E:/merA_2018/BEAST/",merA_subset_trees[i]))
  subset_tree_distances <- cophenetic.phylo(subset_tree)
  subset_tree_hclust <- hclust(as.dist(subset_tree_distances),method="ward.D2")
  subset_tree_50_clusters <- cutree(subset_tree_hclust, k=50)
  subset_tree_50_clusters <- subset_tree_50_clusters[order(as.numeric(gsub("^[A-Z]{3}_([0-9]+)_[0-9]+", "\\1", names(subset_tree_50_clusters))))]
  subset_tree_50_clusters[unique(sort(subset_tree_50_clusters))]
  merA_to_remove[[i]] <- names(subset_tree_50_clusters[duplicated(subset_tree_50_clusters)])
}
lapply(merA_to_remove, write, "E:/merA_2018/BEAST/merA_to_remove.list", append=TRUE, ncolumns=1)

rpoB_subset_trees <- readLines("E:/rpoB_2018/BEAST/rpoB_subset_trees.list")
rpoB_to_remove <- list(NULL)
for (i in 1:length(rpoB_subset_trees)) {
  subset_tree <- read.tree(paste0("E:/rpoB_2018/BEAST/",rpoB_subset_trees[i]))
  subset_tree_distances <- cophenetic.phylo(subset_tree)
  subset_tree_hclust <- hclust(as.dist(subset_tree_distances),method="ward.D2")
  subset_tree_50_clusters <- cutree(subset_tree_hclust, k=50)
  subset_tree_50_clusters <- subset_tree_50_clusters[order(as.numeric(gsub("^[A-Z]{3}_([0-9]+)_[0-9]+", "\\1", names(subset_tree_50_clusters))))]
  subset_tree_50_clusters[unique(sort(subset_tree_50_clusters))]
  rpoB_to_remove[[i]] <- names(subset_tree_50_clusters[duplicated(subset_tree_50_clusters)])
}
lapply(rpoB_to_remove, write, "E:/rpoB_2018/BEAST/rpoB_to_remove.list", append=TRUE, ncolumns=1)

#run BEAST analyses on CSC with skyline plots

#plot skyline plots (with Hg) and trees

roundUp_pop <- function(x,to=10000)
{
  to*(x%/%to + as.logical(x%%to))
}
roundUp_hg <- function(x,to=100)
{
  to*(x%/%to + as.logical(x%%to))
}

BEAST_files <- unique(gsub("[0-9]{2}", "", sample_names(qc_merA_phylo_norm2)))

skyline_data <- list(NULL)
merA_scale_factors <- list(NULL)
for (i in BEAST_files) {
  skyline <- read.csv(paste("E:/merA_2018/BEAST/ref_25_merA_", i,"_skyline.txt",sep = ""),sep = "\t",colClasses = "character")
  colnames(skyline) <-  skyline[1, ]
  skyline <- skyline[-1, ]
  thg_for_skyline <- new_hgdata[grep(i, new_hgdata$NAME),c("Site", "CONC..ug.kg.", "DATE")]
  skyline[, 1:5] <- as.numeric(unlist(skyline[, 1:5]))
  skyline_approx <- data.frame(
    Time = thg_for_skyline$DATE,
    Mean = as.numeric(approx(skyline$Time, skyline$Mean, xout=thg_for_skyline$DATE)$y),
    Median = as.numeric(approx(skyline$Time, skyline$Median, xout=thg_for_skyline$DATE)$y),
    Lower = as.numeric(approx(skyline$Time, skyline$Lower, xout=thg_for_skyline$DATE)$y),
    Upper = as.numeric(approx(skyline$Time, skyline$Upper, xout=thg_for_skyline$DATE)$y),
    Hg = thg_for_skyline$CONC..ug.kg.
  )
  skyline$Hg <- as.numeric(approx(thg_for_skyline$DATE, thg_for_skyline$CONC..ug.kg., xout=skyline$Time)$y)
  skyline <- rbind(skyline, skyline_approx)
  scale_factor <- roundUp_pop(max(skyline$Upper,na.rm=T))
  skyline$Mean.scaled <- skyline$Mean/scale_factor*100
  skyline$Median.scaled <- skyline$Median/scale_factor*100
  skyline$Lower.scaled <- skyline$Lower/scale_factor*100
  skyline$Upper.scaled <- skyline$Upper/scale_factor*100
  skyline_site <- unique(thg_for_skyline$Site)
  skyline$Site <- skyline_site
  
  hg_scale_factor <- roundUp_hg(max(skyline$Hg,na.rm=T))
  skyline$Hg.scaled <- skyline$Hg/hg_scale_factor*100
  #skyline <- melt(skyline)
  skyline_data[[i]] <- skyline
  merA_scale_factors[[i]] <- data.frame(pop=rev(seq(0,scale_factor,length.out=5)), hg=rev(seq(0,hg_scale_factor,length.out=5)))
  }
merA_skyline_data <- do.call(rbind, skyline_data)

image <- ggplot(merA_skyline_data, aes(x = Time, y = Median.scaled)) + geom_line(size = 2) +  geom_line(aes(y = Hg.scaled), color="red", size=1) +
  geom_vline(xintercept = 1800,colour = "blue",linetype = "longdash") + 
  geom_ribbon(aes(ymin = Lower.scaled, ymax = Upper.scaled), alpha = 0.2) + 
  theme(panel.background = element_rect(fill = "white", linetype = "solid",colour = "black"), legend.key = element_rect(fill = "white"), legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_blank(), text = element_text(family = "sans", size = 12), plot.title = element_text(size=20, hjust = 0.5),
        strip.text.y = element_text(angle = 0), axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black")) +
  labs(title = "Inferred population dynamics of merA\nwith [THg] overlaid") + scale_x_continuous(lim = c(1500, 2018)) + scale_y_continuous(lim = c(0,100), sec.axis = sec_axis(~.)) + facet_wrap(.~Site,nrow=4, scales="free_y")

ggsave(file = "E:/merA_2018/results/merA_skyline_plot_new.svg", plot=image, units="in", width=8.5, height=11, scale=1.5)

skyline_data <- list(NULL)
rpoB_scale_factors <- list(NULL)
for (i in BEAST_files) {
  skyline <- read.csv(paste("E:/rpoB_2018/BEAST/ref_25_rpoB_", i,"_skyline.txt",sep = ""),sep = "\t",colClasses = "character")
  colnames(skyline) <-  skyline[1, ]
  skyline <- skyline[-1, ]
  thg_for_skyline <- new_hgdata[grep(i, new_hgdata$NAME),c("Site", "CONC..ug.kg.", "DATE")]
  skyline[, 1:5] <- as.numeric(unlist(skyline[, 1:5]))
  skyline_approx <- data.frame(
    Time = thg_for_skyline$DATE,
    Mean = as.numeric(approx(skyline$Time, skyline$Mean, xout=thg_for_skyline$DATE)$y),
    Median = as.numeric(approx(skyline$Time, skyline$Median, xout=thg_for_skyline$DATE)$y),
    Lower = as.numeric(approx(skyline$Time, skyline$Lower, xout=thg_for_skyline$DATE)$y),
    Upper = as.numeric(approx(skyline$Time, skyline$Upper, xout=thg_for_skyline$DATE)$y),
    Hg = thg_for_skyline$CONC..ug.kg.
  )
  skyline$Hg <- as.numeric(approx(thg_for_skyline$DATE, thg_for_skyline$CONC..ug.kg., xout=skyline$Time)$y)
  skyline <- rbind(skyline, skyline_approx)
  scale_factor <- roundUp_pop(max(skyline$Upper,na.rm=T))
  skyline$Mean.scaled <- skyline$Mean/scale_factor*100
  skyline$Median.scaled <- skyline$Median/scale_factor*100
  skyline$Lower.scaled <- skyline$Lower/scale_factor*100
  skyline$Upper.scaled <- skyline$Upper/scale_factor*100
  skyline_site <- unique(thg_for_skyline$Site)
  skyline$Site <- skyline_site
  
  hg_scale_factor <- roundUp_hg(max(skyline$Hg,na.rm=T))
  skyline$Hg.scaled <- skyline$Hg/hg_scale_factor*100
  #skyline <- melt(skyline)
  skyline_data[[i]] <- skyline
  rpoB_scale_factors[[i]] <- data.frame(pop=rev(seq(0,scale_factor,length.out=5)), hg=rev(seq(0,hg_scale_factor,length.out=5)))
}
rpoB_skyline_data <- do.call(rbind, skyline_data)

image <- ggplot(rpoB_skyline_data, aes(x = Time, y = Median.scaled)) + geom_line(size = 2) +  geom_line(aes(y = Hg.scaled), color="red", size=1) +
  geom_vline(xintercept = 1800,colour = "blue",linetype = "longdash") + 
  geom_ribbon(aes(ymin = Lower.scaled, ymax = Upper.scaled), alpha = 0.2) + 
  theme(panel.background = element_rect(fill = "white", linetype = "solid",colour = "black"), legend.key = element_rect(fill = "white"), legend.position = "none",
        axis.title.x = element_blank(), axis.title.y = element_blank(), text = element_text(family = "sans", size = 12), plot.title = element_text(size=20, hjust = 0.5),
        strip.text.y = element_text(angle = 0), axis.text = element_text(colour = "black"), axis.ticks = element_line(colour = "black")) +
  labs(title = "Inferred population dynamics of rpoB\nwith [THg] overlaid") + scale_x_continuous(lim = c(1500, 2018)) + scale_y_continuous(lim = c(0,100), sec.axis = sec_axis(~.)) + facet_wrap(.~Site,nrow=4, scales="free_y")

ggsave(file = "E:/merA_2018/results/rpoB_skyline_plot_new.svg", plot=image, units="in", width=8.5, height=11, scale=1.5)

#subsample the trees for BaTS analysis
for (i in 1:length(BEAST_files)) {
  #merA
  tr1 <- read.nexus(paste0("E:/merA_2018/BEAST/ref_25_merA_",BEAST_files[i],"_combined.trees"))
  tr1_samp <- sample(tr1, 10000)
  bats_table <- data.frame(label=attributes(tr1_samp)$TipLabel, group=gsub("^.+[0-9]+_","",attributes(tr1_samp)$TipLabel))
  write.tree(tr1_samp, paste0("E:/merA_2018/BEAST/ref_25_merA_",BEAST_files[i],"_subsampled.trees"))
  write.table(bats_table, file=paste0("E:/merA_2018/BEAST/ref_25_merA_",BEAST_files[i],"_combined.trees.txt"), col.names=F, row.names=F, quote=F, sep=" ")
  #rpoB
  tr1 <- read.nexus(paste0("E:/rpoB_2018/BEAST/ref_25_rpoB_",BEAST_files[i],"_combined.trees"))
  tr1_samp <- sample(tr1, 10000)
  bats_table <- data.frame(label=attributes(tr1_samp)$TipLabel, group=gsub("^.+[0-9]+_","",attributes(tr1_samp)$TipLabel))
  write.tree(tr1_samp, paste0("E:/rpoB_2018/BEAST/ref_25_rpoB_",BEAST_files[i],"_subsampled.trees"))
  write.table(bats_table, file=paste0("E:/rpoB_2018/BEAST/ref_25_rpoB_",BEAST_files[i],"_combined.trees.txt"), col.names=F, row.names=F, quote=F, sep=" ")
}
#trees then converted to NEXUS format for BaTS analysis

#plot BEAST trees
better_colors <- brewer.pal(5,"Set1")
better_colors <- c(better_colors, "#8B4513")

for (i in 1:length(BEAST_files)) {
  new_ggtree <-  as_tibble(read.beast(paste0("E:/merA_2018/BEAST/ref_25_merA_",BEAST_files[i],"_combined.trees.tre")))
  groupinfo <- rev(split(new_ggtree$label, gsub("^.+[0-9]+_","",new_ggtree$label)))
  samplingdate <- paste(max(gsub("^.+[0-9]+_","",new_ggtree$label),na.rm=T), "01-01", sep="-")
  new_ggtree <- groupOTU(new_ggtree, groupinfo)
  new_ggtree$group <- factor(new_ggtree$group, levels=names(groupinfo))
  # new_ggtree$branch.length <- new_ggtree$branch.length
  new_ggtree2 <- as.treedata(new_ggtree)

  
  
  image <- ggtree(new_ggtree2, aes(color=group), ladderize = T, size=1.5, mrsd=samplingdate)  + #geom_text2(aes(subset=!isTip, label=round(as.numeric(label),2)), color="grey50", size=3) + 
    ggtitle(paste0("merA ", sites[i])) + theme_tree2() + 
    theme(plot.title = element_text(size= 24, hjust = 0.5),legend.position = "right", legend.title=element_blank(), legend.text = element_text(size = 12), 
          axis.line.x = element_line(color="black", size = 0.5), axis.ticks.x = element_line(size = 1), axis.ticks.length = unit(5,"points"),
          axis.text=element_text(size=14)) + #scale_x_continuous(limits=c(0,2016)) +  
    guides(colour = guide_legend(override.aes = list(size=3))) + scale_color_manual(values=better_colors) 
  ggsave(file = paste0("E:/merA_2018/results/merA_",BEAST_files[i],"_BEAST_tree.svg"), plot=image,
         units="mm", width=250, height=250, scale=0.75)
  new_ggtree <-  as_tibble(read.beast(paste0("E:/rpoB_2018/BEAST/ref_25_rpoB_",BEAST_files[i],"_combined.trees.tre")))
  groupinfo <- rev(split(new_ggtree$label, gsub("^.+[0-9]+_","",new_ggtree$label)))
  samplingdate <- paste(max(gsub("^.+[0-9]+_","",new_ggtree$label),na.rm=T), "01-01", sep="-")
  new_ggtree <- groupOTU(new_ggtree, groupinfo)


  new_ggtree$group <- factor(new_ggtree$group, levels=names(groupinfo))
  # new_ggtree$branch.length <- new_ggtree$branch.length
  new_ggtree2 <- as.treedata(new_ggtree)
  
  image <- ggtree(new_ggtree2, aes(color=group), ladderize = T, size=1.5, mrsd=samplingdate)  + #geom_text2(aes(subset=!isTip, label=round(as.numeric(label),2)), color="grey50", size=3) + 
    ggtitle(paste0("rpoB ", sites[i])) + theme_tree2() + 
    theme(plot.title = element_text(size= 24, hjust = 0.5),legend.position = "right", legend.title=element_blank(), legend.text = element_text(size = 12), 
          axis.line.x = element_line(color="black", size = 0.5), axis.ticks.x = element_line(size = 1), axis.ticks.length = unit(5,"points"),
          axis.text=element_text(size=14)) + #scale_x_continuous(limits=c(0,2016)) + 
    guides(colour = guide_legend(override.aes = list(size=3))) + scale_color_manual(values=better_colors) 
  ggsave(file = paste0("E:/merA_2018/results/rpoB_",BEAST_files[i],"_BEAST_tree.svg"), plot=image,
         units="mm", width=250, height=250, scale=0.75)
}


#random forest models of Skyline data
rf_data <- select(merA_skyline_data, Time, Median, Site)
rf_data <- rf_data[which(rf_data$Time > 1500),]
rf_data <- rf_data[which(complete.cases(rf_data)),]

Time_range <- range(rf_data$Time)
Time_plot_data <- approx(rf_data$Time, rf_data$Median, xout=seq(Time_range[1],Time_range[2], length.out=120))

plot_data <- data.frame(Time=Time_plot_data$x, Site=rep(levels(rf_data$Site)[-1], each=15), Median=rep(0,120))

rf_stats <- list(NULL)
pd_continuous_data <- list(NULL)
pd_categorical_data <- list(NULL)
for (i in 1:10) {
  cf_seed <- 42*i
  set.seed(cf_seed)
  trainIndex <- createDataPartition(rf_data$Site, p = .8, list = FALSE, times = 1)
  train_rf_data <- rf_data[trainIndex,]
  test_rf_data  <- rf_data[-trainIndex,]
  #run random forest model
  set.seed(cf_seed)
  classify <- ranger(Median ~ ., data = train_rf_data, num.trees=5000, importance="impurity")
  test_prediction <- predict(classify, test_rf_data)
  rf_stats[[i]] <- postResample(pred=test_prediction$predictions, obs=test_rf_data$Median)
  
  #gather partial dependence plot data
  pd_data <- partial_dependence(classify, vars=c("Time", "Site"), data=plot_data, n=c(120,nrow(plot_data)))
  pd_continuous <- pd_data[,!(colnames(pd_data) %in% names(Filter(is.factor,pd_data)))]
  pd_continuous <- pd_continuous[!(apply(pd_continuous[-ncol(pd_continuous)], 1, function(x) all(is.na(x)))),]
  types_cont <- c(rep(colnames(pd_continuous)[-ncol(pd_continuous)],as.vector(colSums(!is.na(pd_continuous))[1:(ncol(pd_continuous)-1)])))
  pd_continuous_data[[i]] <- cbind(pd_continuous[ncol(pd_continuous)], value = na.omit(unlist(pd_continuous[-ncol(pd_continuous)])), type=types_cont, rep=i)
  pd_categorical <- cbind(pd_data$Median, Filter(is.factor,pd_data))
  colnames(pd_categorical)[1] <- colnames(pd_data)[ncol(pd_data)]
  pd_categorical <- pd_categorical[!(apply(data.frame(pd_categorical[,-1]), 1, function(x) all(is.na(x)))),]
  types_cat <- c(rep(colnames(pd_categorical)[-1],as.vector(colSums(!is.na(pd_categorical))[2:ncol(pd_categorical)])))
  pd_categorical_data[[i]] <- cbind(pd_categorical[1], value = na.omit(unlist(pd_categorical[-1])), type=types_cat, rep=i)
}

rf_stats <- do.call(rbind, rf_stats)
merA_rf_stats <- round(colMeans(rf_stats),3)
pd_continuous_data <- do.call(rbind, pd_continuous_data)

pd_continuous_time_data <- pd_continuous_data[which(pd_continuous_data$type %in% "Time"),]

time_model <- lm(Median ~ value, data=pd_continuous_time_data)
seg_time_model <- segmented(time_model, seg.Z = ~ value, psi = 1800)
seg_time_fitted <- fitted(seg_time_model) 
seg_time_model_df <- data.frame(Time=pd_continuous_time_data$value[1:120], Median=seg_time_fitted[1:120], type="Time")
seg_time_lines <- data.frame(value=c(seg_time_model$psi[,2], seg_time_model$psi[,2]-seg_time_model$psi[,3]*2.58, seg_time_model$psi[,2]+seg_time_model$psi[,3]*2.58),
                             linetype=c("estimate", "99CI", "99CI"), type="Time")

Clackett_data <- read.csv("E:/merA_2018/Clackett_Hg_data.csv")
Clackett_data <- Clackett_data[which(Clackett_data$robust.biw.mean != "NaN"),]
colnames(Clackett_data) <- c("Time", "mean", "smooth")



image <- ggplot(data = pd_continuous_data, aes(x=value, y=Median, color=factor(rep))) + geom_line(size = 1, alpha=0.2) +
  geom_line(data=seg_time_model_df, aes(x=Time, y=Median), colour="black", size=0.5) + 
  geom_vline(data=seg_time_lines, aes(xintercept = value, linetype = linetype), color="royalblue3", size=0.5) +
  labs(y="Predicted effective population size", x="") + scale_linetype_manual(values=c("dotted", "twodash")) +
  theme(strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.position = "none") + scale_color_manual(values=rep("purple",10))
ggsave(file = "E:/merA_2018/results/merA_pd_continuous.svg", plot=image, units="mm", width=400, height=200)

image <- ggplot(data = Clackett_data, aes(x=Time, y=smooth)) +
  geom_line(colour="orange", size=2) +
  geom_vline(data=seg_time_lines, aes(xintercept = value, linetype = linetype), color="royalblue3", size=0.5) +
  labs(y="Smoothed tree ring [THg]", x="") + scale_linetype_manual(values=c("dotted", "twodash")) +
  theme(strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.position = "none") +
  scale_x_continuous(limits=c(1502.1, 2016)) + scale_y_continuous(labels = scales::number_format(accuracy = 0.00001, decimal.mark = '.'))
ggsave(file = "E:/merA_2018/results/merA_pd_continuous_Hg_plot.svg", plot=image, units="mm", width=400, height=200)

#gather, plot and save categorical variable partial dependency plot
pd_categorical_data <- do.call(rbind, pd_categorical_data)
pd_categorical_data$value <- factor(pd_categorical_data$value, lake_levels)

image <- ggplot(data = pd_categorical_data, aes(x=value, y=Median)) + geom_boxplot(color="purple") +
  labs(y="Predicted effective population size", x="")  + 
  theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size = 16), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file = "E:/merA_2018/results/merA_pd_categorical.svg", plot=image, units="mm", width=400, height=200)

#rpoB
rf_data <- select(rpoB_skyline_data, Time, Median, Site)
rf_data <- rf_data[which(rf_data$Time > 1500),]
rf_data <- rf_data[which(complete.cases(rf_data)),]

rf_stats <- list(NULL)
pd_continuous_data <- list(NULL)
pd_categorical_data <- list(NULL)
for (i in 1:10) {
  cf_seed <- 42*i
  set.seed(cf_seed)
  trainIndex <- createDataPartition(rf_data$Site, p = .8, list = FALSE, times = 1)
  train_rf_data <- rf_data[trainIndex,]
  test_rf_data  <- rf_data[-trainIndex,]
  #run random forest model
  set.seed(cf_seed)
  classify <- ranger(Median ~ ., data = train_rf_data, num.trees=5000, importance="impurity")
  test_prediction <- predict(classify, test_rf_data)
  rf_stats[[i]] <- postResample(pred=test_prediction$predictions, obs=test_rf_data$Median)
  
  #gather partial dependence plot data
  pd_data <- partial_dependence(classify, vars=c("Time", "Site"), data=plot_data, n=c(120,nrow(plot_data)))
  pd_continuous <- pd_data[,!(colnames(pd_data) %in% names(Filter(is.factor,pd_data)))]
  pd_continuous <- pd_continuous[!(apply(pd_continuous[-ncol(pd_continuous)], 1, function(x) all(is.na(x)))),]
  types_cont <- c(rep(colnames(pd_continuous)[-ncol(pd_continuous)],as.vector(colSums(!is.na(pd_continuous))[1:(ncol(pd_continuous)-1)])))
  pd_continuous_data[[i]] <- cbind(pd_continuous[ncol(pd_continuous)], value = na.omit(unlist(pd_continuous[-ncol(pd_continuous)])), type=types_cont, rep=i)
  pd_categorical <- cbind(pd_data$Median, Filter(is.factor,pd_data))
  colnames(pd_categorical)[1] <- colnames(pd_data)[ncol(pd_data)]
  pd_categorical <- pd_categorical[!(apply(data.frame(pd_categorical[,-1]), 1, function(x) all(is.na(x)))),]
  types_cat <- c(rep(colnames(pd_categorical)[-1],as.vector(colSums(!is.na(pd_categorical))[2:ncol(pd_categorical)])))
  pd_categorical_data[[i]] <- cbind(pd_categorical[1], value = na.omit(unlist(pd_categorical[-1])), type=types_cat, rep=i)
}

rf_stats <- do.call(rbind, rf_stats)
rpoB_rf_stats <- round(colMeans(rf_stats),3)
pd_continuous_data <- do.call(rbind, pd_continuous_data)

pd_continuous_time_data <- pd_continuous_data[which(pd_continuous_data$type %in% "Time"),]

time_model <- lm(Median ~ value, data=pd_continuous_time_data)
seg_time_model <- segmented(time_model, seg.Z = ~ value, psi = 1800)
seg_time_fitted <- fitted(seg_time_model) 
seg_time_model_df <- data.frame(Time=pd_continuous_time_data$value[1:120], Median=seg_time_fitted[1:120], type="Time")
seg_time_lines <- data.frame(value=c(seg_time_model$psi[,2], seg_time_model$psi[,2]-seg_time_model$psi[,3]*2.58, seg_time_model$psi[,2]+seg_time_model$psi[,3]*2.58),
                             linetype=c("estimate", "99CI", "99CI"), type="Time")

image <- ggplot(data = pd_continuous_data, aes(x=value, y=Median, color=factor(rep))) + geom_line(size = 1, alpha=0.2) +
  geom_line(data=seg_time_model_df, aes(x=Time, y=Median), colour="black", size=0.5) + geom_vline(data=seg_time_lines, aes(xintercept = value, linetype = linetype), color="royalblue3", size=0.5) +
  labs(y="Predicted effective population size", x="") + scale_linetype_manual(values=c("dotted", "twodash")) +
  theme(strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.position = "none") + scale_color_manual(values=rep("#99566b",10))
ggsave(file = "E:/merA_2018/results/rpoB_pd_continuous.svg", plot=image, units="mm", width=400, height=200)

###As part of the review process, an alternative Figure 4B with pCO2 data was plotted###
#Historical CO2 data from https://doi.org/10.5194/gmd-10-2057-2017

Meinshausen_data <- read.csv("E:/merA_2018/Meinshausen_CO2_data.csv")
colnames(Meinshausen_data) <- c("Time", "Global", "Northern_hemisphere", "Southern_hemisphere")

image <- ggplot(data = Meinshausen_data, aes(x=Time, y=Northern_hemisphere*1000)) + #easiest fix to make the area the same size as in the rpoB plot is to make the labels as long by multiplying them by 1000.
  geom_line(colour="turquoise3", size=2) +
  geom_vline(data=seg_time_lines, aes(xintercept = value, linetype = linetype), color="royalblue3", size=0.5) +
  labs(y="Northern hemispheric CO2 (ppm)", x="") + scale_linetype_manual(values=c("dotted", "twodash")) +
  theme(strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.position = "none") +
  scale_x_continuous(limits=c(1500.041, 2016)) + scale_y_continuous()
ggsave(file = "E:/merA_2018/results/rpoB_pd_continuous_CO2_plot.svg", plot=image, units="mm", width=400, height=200)
###Alternative Figure 4B code ends here###

#gather, plot and save categorical variable partial dependency plot
pd_categorical_data <- do.call(rbind, pd_categorical_data)
pd_categorical_data$value <- factor(pd_categorical_data$value, lake_levels)

image <- ggplot(data = pd_categorical_data, aes(x=value, y=Median)) + geom_boxplot(color="#99566b") +
  labs(y="Predicted effective population size", x="")  + 
  theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size = 16), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file = "E:/merA_2018/results/rpoB_pd_categorical.svg", plot=image, units="mm", width=400, height=200)


#the same with Hg in the models (input group means by site for missing values)
rf_data <- select(merA_skyline_data, Time, Median, Site, Hg)
rf_data <- rf_data[which(rf_data$Time > 1500),]
rf_data <- merge(rf_data, aggregate(Hg ~ Site, rf_data, mean), by="Site")
rf_data$Hg.x[which(is.na(rf_data$Hg.x))] <- rf_data$Hg.y[which(is.na(rf_data$Hg.x))]
rf_data <- rf_data[-5]
colnames(rf_data)[4] <- "Hg"
rf_data <- rf_data[which(complete.cases(rf_data)),]

Hg_range <- range(rf_data$Hg)
Hg_plot_data <- approx(rf_data$Hg, rf_data$Median, xout=seq(Hg_range[1],Hg_range[2], length.out=120))
Time_range <- range(rf_data$Time)
Time_plot_data <- approx(rf_data$Time, rf_data$Median, xout=seq(Time_range[1],Time_range[2], length.out=120))

plot_data <- data.frame(Time=Time_plot_data$x, Hg=Hg_plot_data$x, Site=rep(levels(rf_data$Site)[-1], each=15), Median=rep(0,120))
#conduct 10 random samplings
rf_stats <- list(NULL)
pd_continuous_data <- list(NULL)
pd_categorical_data <- list(NULL)
for (i in 1:10) {
  cf_seed <- 42*i
  set.seed(cf_seed)
  trainIndex <- createDataPartition(rf_data$Site, p = .8, list = FALSE, times = 1)
  train_rf_data <- rf_data[trainIndex,]
  test_rf_data  <- rf_data[-trainIndex,]
  #run random forest model
  set.seed(cf_seed)
  classify <- ranger(Median ~ ., data = train_rf_data, num.trees=5000, importance="impurity")
  test_prediction <- predict(classify, test_rf_data)
  rf_stats[[i]] <- postResample(pred=test_prediction$predictions, obs=test_rf_data$Median)
  
  #gather partial dependence plot data
  pd_data <- partial_dependence(classify, vars=c("Time", "Hg", "Site"), data=plot_data, n=c(120,nrow(plot_data)))
  pd_continuous <- pd_data[,!(colnames(pd_data) %in% names(Filter(is.factor,pd_data)))]
  pd_continuous <- pd_continuous[!(apply(pd_continuous[-ncol(pd_continuous)], 1, function(x) all(is.na(x)))),]
  types_cont <- c(rep(colnames(pd_continuous)[-ncol(pd_continuous)],as.vector(colSums(!is.na(pd_continuous))[1:(ncol(pd_continuous)-1)])))
  pd_continuous_data[[i]] <- cbind(pd_continuous[ncol(pd_continuous)], value = na.omit(unlist(pd_continuous[-ncol(pd_continuous)])), type=types_cont, rep=i)
  pd_categorical <- cbind(pd_data$Median, Filter(is.factor,pd_data))
  colnames(pd_categorical)[1] <- colnames(pd_data)[ncol(pd_data)]
  pd_categorical <- pd_categorical[!(apply(data.frame(pd_categorical[,-1]), 1, function(x) all(is.na(x)))),]
  types_cat <- c(rep(colnames(pd_categorical)[-1],as.vector(colSums(!is.na(pd_categorical))[2:ncol(pd_categorical)])))
  pd_categorical_data[[i]] <- cbind(pd_categorical[1], value = na.omit(unlist(pd_categorical[-1])), type=types_cat, rep=i)
}

rf_stats <- do.call(rbind, rf_stats)
merA_rf_stats_alt <- round(colMeans(rf_stats),3)
pd_continuous_data <- do.call(rbind, pd_continuous_data)

#fit a segmented linear regression models for date and Hg and find breakpoints
pd_continuous_time_data <- pd_continuous_data[which(pd_continuous_data$type %in% "Time"),]

time_model <- lm(Median ~ value, data=pd_continuous_time_data)
#seg_time_model <- segmented(time_model, seg.Z = ~ value, psi = 1800)
seg_time_model <- segmented(time_model, seg.Z = ~ value, psi = c(1800,2000))
seg_time_fitted <- fitted(seg_time_model) 
seg_time_model_df <- data.frame(Time=pd_continuous_time_data$value[1:120], Median=seg_time_fitted[1:120], type="Time")
seg_time_lines <- data.frame(value=c(seg_time_model$psi[,2], seg_time_model$psi[,2]-seg_time_model$psi[,3]*2.58, seg_time_model$psi[,2]+seg_time_model$psi[,3]*2.58),
                             linetype=rep(c("estimate", "99CI", "99CI"), each=2), type="Time")

pd_continuous_Hg_data <- pd_continuous_data[which(pd_continuous_data$type %in% "Hg"),]

Hg_model <- lm(Median ~ value, data=pd_continuous_Hg_data)
seg_Hg_model <- segmented(Hg_model, seg.Z = ~ value, psi = c(50,1000))
seg_Hg_fitted <- fitted(seg_Hg_model) 
seg_Hg_model_df <- data.frame(Time=pd_continuous_Hg_data$value[1:120], Median=seg_Hg_fitted[1:120], type="Hg")
seg_Hg_lines <- data.frame(value=c(seg_Hg_model$psi[,2], seg_Hg_model$psi[,2]-seg_Hg_model$psi[,3]*2.58, seg_Hg_model$psi[,2]+seg_Hg_model$psi[,3]*2.58),
                           linetype=rep(c("estimate", "99CI", "99CI"), each=2), type="Hg")#, line=rep(c("1", "2", "3"), each=3))

seg_lines <- rbind(seg_time_lines, seg_Hg_lines)
seg_model_df <- rbind(seg_time_model_df, seg_Hg_model_df)

image <- ggplot(data = pd_continuous_data, aes(x=value, y=Median, color=factor(rep))) + geom_line(size = 1, alpha=0.2) +
  geom_line(data=seg_model_df, aes(x=Time, y=Median), colour="black", size=0.5) + geom_vline(data=seg_lines, aes(xintercept = value, linetype = linetype), color="royalblue3", size=0.5) +
  labs(y="Predicted effective population size", x="") + scale_linetype_manual(values=c("dotted", "twodash")) + facet_wrap(.~type, nrow=2, scales="free_x") +
  theme(strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.position = "none") + scale_color_manual(values=rep("purple",10))
ggsave(file = "E:/merA_2018/results/merA_pd_continuous_alt.svg", plot=image, units="mm", width=400, height=200)

#gather, plot and save categorical variable partial dependency plot
pd_categorical_data <- do.call(rbind, pd_categorical_data)
pd_categorical_data$value <- factor(pd_categorical_data$value, lake_levels)
  
image <- ggplot(data = pd_categorical_data, aes(x=value, y=Median)) + geom_boxplot(color="purple") +
  labs(y="Predicted effective population size", x="")  + 
  theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size = 16), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file = "E:/merA_2018/results/merA_pd_categorical_alt.svg", plot=image, units="mm", width=400, height=200)

#rpoB
rf_data <- select(rpoB_skyline_data, Time, Median, Site, Hg)
rf_data <- rf_data[which(rf_data$Time > 1500),]
rf_data <- merge(rf_data, aggregate(Hg ~ Site, rf_data, mean), by="Site")
rf_data$Hg.x[which(is.na(rf_data$Hg.x))] <- rf_data$Hg.y[which(is.na(rf_data$Hg.x))]
rf_data <- rf_data[-5]
colnames(rf_data)[4] <- "Hg"
rf_data <- rf_data[which(complete.cases(rf_data)),]

#conduct 10 random samplings
rf_stats <- list(NULL)
pd_continuous_data <- list(NULL)
pd_categorical_data <- list(NULL)
for (i in 1:10) {
  cf_seed <- 42*i
  set.seed(cf_seed)
  trainIndex <- createDataPartition(rf_data$Site, p = .8, list = FALSE, times = 1)
  train_rf_data <- rf_data[trainIndex,]
  test_rf_data  <- rf_data[-trainIndex,]
  #run random forest model
  set.seed(cf_seed)
  classify <- ranger(Median ~ ., data = train_rf_data, num.trees=5000, importance="impurity")
  test_prediction <- predict(classify, test_rf_data)
  rf_stats[[i]] <- postResample(pred=test_prediction$predictions, obs=test_rf_data$Median)
  
  #gather partial dependence plot data
  pd_data <- partial_dependence(classify, vars=c("Time", "Hg", "Site"), data=plot_data, n=c(120,nrow(plot_data)))
  pd_continuous <- pd_data[,!(colnames(pd_data) %in% names(Filter(is.factor,pd_data)))]
  pd_continuous <- pd_continuous[!(apply(pd_continuous[-ncol(pd_continuous)], 1, function(x) all(is.na(x)))),]
  types_cont <- c(rep(colnames(pd_continuous)[-ncol(pd_continuous)],as.vector(colSums(!is.na(pd_continuous))[1:(ncol(pd_continuous)-1)])))
  pd_continuous_data[[i]] <- cbind(pd_continuous[ncol(pd_continuous)], value = na.omit(unlist(pd_continuous[-ncol(pd_continuous)])), type=types_cont, rep=i)
  pd_categorical <- cbind(pd_data$Median, Filter(is.factor,pd_data))
  colnames(pd_categorical)[1] <- colnames(pd_data)[ncol(pd_data)]
  pd_categorical <- pd_categorical[!(apply(data.frame(pd_categorical[,-1]), 1, function(x) all(is.na(x)))),]
  types_cat <- c(rep(colnames(pd_categorical)[-1],as.vector(colSums(!is.na(pd_categorical))[2:ncol(pd_categorical)])))
  pd_categorical_data[[i]] <- cbind(pd_categorical[1], value = na.omit(unlist(pd_categorical[-1])), type=types_cat, rep=i)
}

rf_stats <- do.call(rbind, rf_stats)
rpoB_rf_stats_alt <- round(colMeans(rf_stats),3)
pd_continuous_data <- do.call(rbind, pd_continuous_data)

#fit a segmented linear regression models for date and Hg and find breakpoints
pd_continuous_time_data <- pd_continuous_data[which(pd_continuous_data$type %in% "Time"),]

time_model <- lm(Median ~ value, data=pd_continuous_time_data)
#seg_time_model <- segmented(time_model, seg.Z = ~ value, psi = 1800)
seg_time_model <- segmented(time_model, seg.Z = ~ value, psi = c(1800,2000))
seg_time_fitted <- fitted(seg_time_model) 
seg_time_model_df <- data.frame(Time=pd_continuous_time_data$value[1:120], Median=seg_time_fitted[1:120], type="Time")
seg_time_lines <- data.frame(value=c(seg_time_model$psi[,2], seg_time_model$psi[,2]-seg_time_model$psi[,3]*2.58, seg_time_model$psi[,2]+seg_time_model$psi[,3]*2.58),
                             linetype=rep(c("estimate", "99CI", "99CI"), each=2), type="Time")

pd_continuous_Hg_data <- pd_continuous_data[which(pd_continuous_data$type %in% "Hg"),]

Hg_model <- lm(Median ~ value, data=pd_continuous_Hg_data)
seg_Hg_model <- segmented(Hg_model, seg.Z = ~ value, psi = c(50,1000))
seg_Hg_fitted <- fitted(seg_Hg_model) 
seg_Hg_model_df <- data.frame(Time=pd_continuous_Hg_data$value[1:120], Median=seg_Hg_fitted[1:120], type="Hg")
seg_Hg_lines <- data.frame(value=c(seg_Hg_model$psi[,2], seg_Hg_model$psi[,2]-seg_Hg_model$psi[,3]*2.58, seg_Hg_model$psi[,2]+seg_Hg_model$psi[,3]*2.58),
                           linetype=rep(c("estimate", "99CI", "99CI"), each=2), type="Hg")

seg_lines <- rbind(seg_time_lines, seg_Hg_lines)
seg_model_df <- rbind(seg_time_model_df, seg_Hg_model_df)

image <- ggplot(data = pd_continuous_data, aes(x=value, y=Median, color=factor(rep))) + geom_line(size = 1, alpha=0.2) +
  geom_line(data=seg_model_df, aes(x=Time, y=Median), colour="black", size=0.5) + geom_vline(data=seg_lines, aes(xintercept = value, linetype = linetype), color="royalblue3", size=0.5) +
  labs(y="Predicted effective population size", x="") + scale_linetype_manual(values=c("dotted", "twodash")) + facet_wrap(.~type, nrow=2, scales="free_x") +
  theme(strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16), legend.position = "none") + scale_color_manual(values=rep("#99566b",10))
ggsave(file = "E:/merA_2018/results/rpoB_pd_continuous_alt.svg", plot=image, units="mm", width=400, height=200)

#gather, plot and save categorical variable partial dependency plot
pd_categorical_data <- do.call(rbind, pd_categorical_data)
pd_categorical_data$value <- factor(pd_categorical_data$value, lake_levels)

image <- ggplot(data = pd_categorical_data, aes(x=value, y=Median)) + geom_boxplot(color="#99566b") +
  labs(y="Predicted effective population size", x="")  + 
  theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_text(size = 16), strip.text.x = element_text(size = 16), axis.title.y = element_text(size = 16))
ggsave(file = "E:/merA_2018/results/rpoB_pd_categorical_alt.svg", plot=image, units="mm", width=400, height=200)

#analyze ddPCR data
DNA_conc <- read.csv("E:/merA_2018/FIN_CAN_DNA_conc_mixed.csv", sep = ",", row.names = NULL)
DNA_conc <- ddply(DNA_conc,~SAMPLE,summarise,avg.DNA.CONC.nguL = round(mean(DNA.CONC.nguL),2))
#fix names for compatibility with hgdata
DNA_conc$NAME <- gsub("*_[C-Z]*[0-9]*", "", DNA_conc$SAMPLE)
DNA_conc$NAME <- as.character(DNA_conc$NAME)
DNA_conc$NAME[which(grepl("[A-Z]{3}[0-9]{1}$", DNA_conc$NAME))] <- paste0(substr(DNA_conc$NAME[which(grepl("[A-Z]{3}[0-9]{1}$", DNA_conc$NAME))], 1, 3), 0, substr(DNA_conc$NAME[which(grepl("[A-Z]{3}[0-9]{1}$", DNA_conc$NAME))], 4, 4))
new_hgdata2 <- merge(new_hgdata, DNA_conc, by = "NAME", all = T)

merA_ddPCR1 <- read.csv("E:/merA_2018/June_8th_qmerA_rows_1-9.csv", sep = ",", row.names = NULL)
#colnames(merA_ddPCR1) <- colnames(merA_ddPCR1)[-1]
merA_ddPCR1 <- subset(merA_ddPCR1, select=c("Well", "Concentration", "PoissonConfMin", "PoissonConfMax"))
#remove everything with less copies than the negative control
merA_ddPCR1 <- merA_ddPCR1[merA_ddPCR1$Concentration > 0.2,]
#correct concentrations for 1 µL (concentrations given this far were per µL in the reaction, which were originally 2 µL of original sample diluted into a 22 µL reaction)
merA_ddPCR1$Concentration <- merA_ddPCR1$Concentration*11
merA_ddPCR1$PoissonConfMin <- merA_ddPCR1$PoissonConfMin*11
merA_ddPCR1$PoissonConfMax <- merA_ddPCR1$PoissonConfMax*11
merA_ddPCR1$GENE <- "merA"
merA_ddPCR2 <- read.csv("E:/merA_2018/June_8th_qmerA_rows_10-17.csv", sep = ",", row.names = NULL)
colnames(merA_ddPCR2) <- colnames(merA_ddPCR2)[-1]
merA_ddPCR2 <- subset(merA_ddPCR2, select=c("Well", "Concentration", "PoissonConfMin", "PoissonConfMax"))
#remove everything with less copies than the negative control
merA_ddPCR2 <- merA_ddPCR2[merA_ddPCR2$Concentration > 0.2,]
#correct concentrations for 1 µL (concentrations given this far were per µL in the reaction, which were originally 2 µL of original sample diluted into a 22 µL reaction)
merA_ddPCR1$Concentration <- merA_ddPCR1$Concentration*11
merA_ddPCR1$PoissonConfMin <- merA_ddPCR1$PoissonConfMin*11
merA_ddPCR1$PoissonConfMax <- merA_ddPCR1$PoissonConfMax*11
merA_ddPCR2$GENE <- "merA"

glnA_ddPCR1 <- read.csv("E:/merA_2018/June_9th_glnA_rows_1-9.csv", sep = ",", row.names = NULL)
colnames(glnA_ddPCR1) <- colnames(glnA_ddPCR1)[-1]
glnA_ddPCR1 <- subset(glnA_ddPCR1, select=c("Well", "Concentration", "PoissonConfMin", "PoissonConfMax"))
#remove everything with less copies than the negative control
glnA_ddPCR1 <- glnA_ddPCR1[glnA_ddPCR1$Concentration > 0.56,]
#correct concentrations for 1 µL (concentrations given this far were per µL in the reaction, which were originally 2 µL of original sample diluted into a 22 µL reaction, and glnA in addition diluted 1:10)
glnA_ddPCR1$Concentration <- glnA_ddPCR1$Concentration*110
glnA_ddPCR1$PoissonConfMin <- glnA_ddPCR1$PoissonConfMin*110
glnA_ddPCR1$PoissonConfMax <- glnA_ddPCR1$PoissonConfMax*110
glnA_ddPCR1$GENE <- "glnA"
#glnA_ddPCR1$Concentration <- as.numeric(levels(glnA_ddPCR1$Concentration))[glnA_ddPCR1$Concentration]
glnA_ddPCR2 <- read.csv("E:/merA_2018/June_9th_glnA_rows_10-17.csv", sep = ",", row.names = NULL)
colnames(glnA_ddPCR2) <- colnames(glnA_ddPCR2)[-1]
glnA_ddPCR2 <- subset(glnA_ddPCR2, select=c("Well", "Concentration", "PoissonConfMin", "PoissonConfMax"))
#remove everything with less copies than the negative control
glnA_ddPCR2 <- glnA_ddPCR2[glnA_ddPCR2$Concentration > 0.56,]
#correct concentrations for 1 µL (concentrations given this far were per µL in the reaction, which were originally 2 µL of original sample diluted into a 22 µL reaction, and glnA in addition diluted 1:10)
glnA_ddPCR2$Concentration <- glnA_ddPCR2$Concentration*110
glnA_ddPCR2$PoissonConfMin <- glnA_ddPCR2$PoissonConfMin*110
glnA_ddPCR2$PoissonConfMax <- glnA_ddPCR2$PoissonConfMax*110
glnA_ddPCR2$GENE <- "glnA"
#glnA_ddPCR2$Concentration <- as.numeric(levels(glnA_ddPCR2$Concentration))[glnA_ddPCR2$Concentration]

ddPCR1 <- rbind(merA_ddPCR1, glnA_ddPCR1)
ddPCR1$PLATE <- 1
ddPCR2 <- rbind(merA_ddPCR2, glnA_ddPCR2)
ddPCR2$PLATE <- 2
ddPCR <- rbind(ddPCR1, ddPCR2)

new_hgdata2 <- merge(new_hgdata2, ddPCR, by = c("PLATE","Well"), all = T)
#normalize counts for ng DNA
new_hgdata2$Concentration <- new_hgdata2$Concentration/new_hgdata2$avg.DNA.CONC.nguL
new_hgdata2$PoissonConfMin <- new_hgdata2$PoissonConfMin/new_hgdata2$avg.DNA.CONC.nguL
new_hgdata2$PoissonConfMax <- new_hgdata2$PoissonConfMax/new_hgdata2$avg.DNA.CONC.nguL

new_hgdata2 <- new_hgdata2[!is.na(new_hgdata2$Site),]

max_hg <- max(new_hgdata2$CONC..ug.kg.)
max_conc <- max(new_hgdata2$Concentration[!is.na(new_hgdata2$Concentration)])

new_hgdata2$CONC..ug.kg._norm <- new_hgdata2$CONC..ug.kg. / (max_hg / max_conc)

new_hgdata_lim <- new_hgdata2[which(new_hgdata2$DATE > 1750),]

#plot tot-Hg, ddPCR results and sequencing for all the cores
image <- ggplot(new_hgdata_lim, aes(x=DATE, y=Concentration, group=Site)) + geom_point(aes(fill = GENE), size=3, shape = 21, color = "white") + geom_point(aes(y=CONC..ug.kg._norm, color = PCR.BAND, shape = Sequencing), size=4) + 
  geom_errorbar(aes(ymin=PoissonConfMin, ymax=PoissonConfMax), size=1, color = "black", width=0) +
  geom_line(aes(y=CONC..ug.kg._norm), color="black") + facet_wrap(~Site, scales="free_x") + coord_flip() +
  xlab("Date (CE)") + ylab("Copies per ng DNA") + ggtitle("[THg] profiles and gene copies in sediment cores") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_log10(sec.axis = sec_axis(trans = ~ . * (max_hg / max_conc), name = "[THg] ng / g", breaks = c(10,100,1000,2000)), breaks = c(0.1,1,10,100,1000,10000)) +
  scale_fill_manual(values = c("green", "purple"), name = "Gene") + scale_color_manual(values = c("deeppink1", "black", "grey85"), name = "PCR band") +
  scale_shape_manual(values = c(0, 15))
ggsave(file = "E:/merA_2018/results/Dating_Hg_ddPCR_sequencing_overview.svg", plot=image, units="mm", width=500, height=400)

#plot dating info for all the cores
dating_data_lim <- unique(new_hgdata_lim[c("Site", "MIDDLE.DEPTH", "DATE", "CRS.ERROR", "Sequencing")])
dating_data_lim$measured <- ifelse(is.na(dating_data_lim$CRS.ERROR), "Estimated", "Primary")
#dating_data_lim$measured[which(grepl("Aquatuk Lake",dating_data_lim$Site))] <- "Primary"
dating_data_lim$measured[which(grepl("Hazen",dating_data_lim$Site))] <- "Previous data"
dating_data_lim$measured[which(grepl("Harjavalta",dating_data_lim$Site))] <- "Estimated"
dating_data_lim$measured[which(grepl("Pocket",dating_data_lim$Site))] <- ifelse(dating_data_lim$measured[which(grepl("Pocket",dating_data_lim$Site))] %in% "Primary", "Previous data", "Estimated")

image <- ggplot(dating_data_lim, aes(y=DATE, x=MIDDLE.DEPTH, group=Site)) + geom_point(aes(color = measured, shape=Sequencing), size= 4) + coord_flip() + scale_x_reverse() + 
  geom_errorbar(aes(ymin=DATE-CRS.ERROR, ymax=DATE+CRS.ERROR), size=1, color = "black", width=0) + facet_wrap(~Site, scales="free") + 
  scale_color_manual(values=c("#f87341","#875fb4","#558d2a"), name="Measurement") + guides(color = guide_legend(reverse=T)) + scale_shape_manual(values=c(0, 15)) +
  scale_y_continuous(limits=c(1750,2020), breaks=seq(1750,2000,by=50)) + ylab("Date (CE)") + xlab("Sediment depth from surface (cm)")
ggsave(file = "E:/merA_2018/results/Dating.svg", plot=image, units="mm", width=500, height=400)

#model glnA and merA copy numbers on dating and tot-Hg (separate models)

center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

glnA_ddPCR_data <- new_hgdata_lim[which(new_hgdata_lim$GENE %in% "glnA"),c(4,6,11,15,16,19,24,20)]
glnA_ddPCR_data <- glnA_ddPCR_data[which(complete.cases(glnA_ddPCR_data)),]
#glnA_ddPCR_data2 <- new_hgdata2[which(new_hgdata2$GENE %in% "glnA"),c(4,6,11,15,16,19,24,21)]
#glnA_ddPCR_data2 <- plyr::rename(glnA_ddPCR_data2,c("PoissonConfMin"="Concentration"))
# glnA_ddPCR_data3 <- new_hgdata2[which(new_hgdata2$GENE %in% "glnA"),c(4,6,11,15,16,19,24,22)]
# glnA_ddPCR_data3 <- plyr::rename(glnA_ddPCR_data3,c("PoissonConfMax"="Concentration"))
#remove Harjavalta samples and Pulmankijärvi 9.5 and 10.5 cm samples (outliers)
glnA_ddPCR_data <- glnA_ddPCR_data[!(glnA_ddPCR_data$Site %in% "Kokemäenjoki / Harjavalta"),]
#remove outlier data points
glnA_ddPCR_data <- glnA_ddPCR_data[-which(grepl( "12294|8880",  glnA_ddPCR_data$Concentration)),]

glnA_rescaled <- data.frame(lapply(glnA_ddPCR_data[-1], function(x) scale(x, center = T, scale = max(x, na.rm = TRUE)/1)))
glnA_ddPCR_data <- cbind(glnA_ddPCR_data[1], glnA_rescaled)
# #normalized gene copy numbers need to be log-transformed: these are the best fits
# qqp(log10(glnA_ddPCR_data1$Concentration), "norm")
# qqp(log10(glnA_ddPCR_data1$CONC..ug.kg.), "norm")
# qqp(log10(2018-glnA_ddPCR_data1$DATE), "norm")
# qqp(log10(glnA_ddPCR_data1$MIDDLE.DEPTH), "norm")

glnA_ddPCR_model <- lmer(log10(Concentration) ~ CONC..ug.kg. + DATE + MIDDLE.DEPTH + ( 1 | Site), data=glnA_ddPCR_data)
glnA_ddPCR_model.null <- lmer(log10(Concentration) ~ ( 1 | Site), data=glnA_ddPCR_data)
anova(glnA_ddPCR_model, glnA_ddPCR_model.null)
summary(glnA_ddPCR_model)
qqnorm(residuals(glnA_ddPCR_model))

glnA_ddPCR_model_confint <- data.frame(confint(glnA_ddPCR_model))
glnA_ddPCR_model_confint <- glnA_ddPCR_model_confint[-(1:3),]
glnA_ddPCR_model_confint$estimate <- fixef(glnA_ddPCR_model)[-1]
glnA_ddPCR_model_confint$tested_variable <- c("[THg]", "Date", "Sediment depth from surface")

image <- ggplot(glnA_ddPCR_model_confint, aes(y=tested_variable, x=estimate, xmin=X2.5.., xmax=X97.5..)) + geom_point(shape=18, size=5) +
  geom_errorbarh(height=.1) + geom_vline(xintercept=0, color='black', linetype='dashed') + 
  theme(axis.title.y = element_blank(), plot.title = element_text(size=20, hjust = 0.5)) +
  ggtitle("Normalized effect sizes (glnA copy number)") + xlab("")
ggsave(file = "E:/merA_2018/results/rpoB_ddCPR_CI.svg", plot=image, units="mm", width=250, height=200)


#seems that sediment depth is the only significant variable
#plot middle depth effect size for rpoB

test_MIDDLE.DEPTH <- round(seq(range(glnA_ddPCR_data$MIDDLE.DEPTH)[1], range(glnA_ddPCR_data$MIDDLE.DEPTH)[2], length.out=10),2)
glnA_MIDDLE.DEPTH_newdata <- expand.grid(MIDDLE.DEPTH = scale(test_MIDDLE.DEPTH, center = T, scale = max(test_MIDDLE.DEPTH, na.rm = TRUE)/1),
                                         CONC..ug.kg. = 0,
                                         DATE = 0)
glnA_MIDDLE.DEPTH_newdata$Concentration <- predict(glnA_ddPCR_model, glnA_MIDDLE.DEPTH_newdata, re.form=NA)
mm <- model.matrix(terms(glnA_ddPCR_model), glnA_MIDDLE.DEPTH_newdata)
pvar1 <- diag(mm %*% tcrossprod(vcov(glnA_ddPCR_model),mm))
tvar1 <- pvar1+VarCorr(glnA_ddPCR_model)$Site[1]  ## must be adapted for more complex models
cmult <- 2 ## could use 1.96
glnA_MIDDLE.DEPTH_newdata <- data.frame(
  glnA_MIDDLE.DEPTH_newdata, 
  plo = glnA_MIDDLE.DEPTH_newdata$Concentration-cmult*sqrt(pvar1), 
  phi = glnA_MIDDLE.DEPTH_newdata$Concentration+cmult*sqrt(pvar1), 
  tlo = glnA_MIDDLE.DEPTH_newdata$Concentration-cmult*sqrt(tvar1), 
  thi = glnA_MIDDLE.DEPTH_newdata$Concentration+cmult*sqrt(tvar1)
)

glnA_MIDDLE.DEPTH_newdata$orig.MIDDLE.DEPTH <- test_MIDDLE.DEPTH
glnA_copies_scales <-  scale(glnA_ddPCR_data$Concentration, center = T, scale = max(glnA_ddPCR_data$Concentration, na.rm = TRUE)/1)
glnA_MIDDLE.DEPTH_newdata$orig.Concentration <- 10^glnA_MIDDLE.DEPTH_newdata$Concentration * attr(glnA_copies_scales, 'scaled:scale') + attr(glnA_copies_scales, 'scaled:center')
glnA_MIDDLE.DEPTH_newdata$orig.plo <- 10^glnA_MIDDLE.DEPTH_newdata$plo * attr(glnA_copies_scales, 'scaled:scale') + attr(glnA_copies_scales, 'scaled:center')
glnA_MIDDLE.DEPTH_newdata$orig.phi <- 10^glnA_MIDDLE.DEPTH_newdata$phi * attr(glnA_copies_scales, 'scaled:scale') + attr(glnA_copies_scales, 'scaled:center')
glnA_MIDDLE.DEPTH_newdata$orig.tlo <- 10^glnA_MIDDLE.DEPTH_newdata$tlo * attr(glnA_copies_scales, 'scaled:scale') + attr(glnA_copies_scales, 'scaled:center')
glnA_MIDDLE.DEPTH_newdata$orig.thi <- 10^glnA_MIDDLE.DEPTH_newdata$thi * attr(glnA_copies_scales, 'scaled:scale') + attr(glnA_copies_scales, 'scaled:center')

image <- ggplot(glnA_MIDDLE.DEPTH_newdata, aes(x=orig.MIDDLE.DEPTH, y=orig.Concentration)) + 
  coord_flip() + scale_x_reverse() + scale_y_log10(limits = c(1e3,1e7)) +
  geom_ribbon(aes(ymin = orig.tlo, ymax = orig.thi), colour="lightgrey", fill="lightgrey") +
  geom_line(lwd=2) + theme(axis.text.y = element_text(hjust = 0, size = 22), plot.title = element_text(size=28, hjust = 0.5)) +
  labs(title="95% CI based on FE uncertainty + RE variance") +
  xlab("Sediment depth (cm)") + ylab("Modeled glnA copies per ng DNA")
ggsave(file = "E:/merA_2018/results/rpoB_lme_effect_size.svg", plot=image, units="mm", width=250, height=200)

merA_ddPCR_data1 <- new_hgdata_lim[which(new_hgdata_lim$GENE %in% "merA"),c(4,6,11,15,16,19,24,20)]
merA_ddPCR_data1 <- merA_ddPCR_data1[which(complete.cases(merA_ddPCR_data1)),]
#merA_ddPCR_data2 <- new_hgdata2[which(new_hgdata2$GENE %in% "merA"),c(4,6,11,15,16,19,24,21)]
#merA_ddPCR_data2 <- plyr::rename(merA_ddPCR_data2,c("PoissonConfMin"="Concentration"))
# merA_ddPCR_data3 <- new_hgdata2[which(new_hgdata2$GENE %in% "merA"),c(4,6,11,15,16,19,24,22)]
# merA_ddPCR_data3 <- plyr::rename(merA_ddPCR_data3,c("PoissonConfMax"="Concentration"))
merA_ddPCR_data_alt <- merA_ddPCR_data1
merA_ddPCR_data1 <- merA_ddPCR_data1[!(merA_ddPCR_data1$Site %in% "Kokemäenjoki / Harjavalta"),]

merA_rescaled <- data.frame(lapply(merA_ddPCR_data1[-1], function(x) scale(x, center = T, scale = max(x, na.rm = TRUE)/1)))
merA_ddPCR_data1 <- cbind(merA_ddPCR_data1[1], merA_rescaled)

# #normalized gene copy numbers need to be log-transformed: these are the best fits
# qqp(log10(merA_ddPCR_data1$Concentration), "norm")
# qqp(log10(merA_ddPCR_data1$CONC..ug.kg.), "norm")
# qqp(log10(2018-merA_ddPCR_data1$DATE), "norm")
# qqp(log10(merA_ddPCR_data1$MIDDLE.DEPTH), "norm")

merA_ddPCR_model <- lmer(log10(Concentration) ~ CONC..ug.kg. + DATE + MIDDLE.DEPTH + ( 1 | Site), data=merA_ddPCR_data1)
merA_ddPCR_model.null <- lmer(log10(Concentration) ~ ( 1 | Site), data=merA_ddPCR_data1)
anova(merA_ddPCR_model, merA_ddPCR_model.null)
summary(merA_ddPCR_model)
qqnorm(residuals(merA_ddPCR_model))

merA_ddPCR_model_confint <- data.frame(confint(merA_ddPCR_model))
merA_ddPCR_model_confint <- merA_ddPCR_model_confint[-(1:3),]
merA_ddPCR_model_confint$estimate <- fixef(merA_ddPCR_model)[-1]
merA_ddPCR_model_confint$tested_variable <- c("[THg]", "Date", "Sediment depth from surface")

image <- ggplot(merA_ddPCR_model_confint, aes(y=tested_variable, x=estimate, xmin=X2.5.., xmax=X97.5..)) + geom_point(shape=18, size=5) +
  geom_errorbarh(height=.1) + geom_vline(xintercept=0, color='black', linetype='dashed') + 
  theme(axis.title.y = element_blank(), plot.title = element_text(size=20, hjust = 0.5)) +
  ggtitle("Normalized effect sizes (merA copy number)") + xlab("")
ggsave(file = "E:/merA_2018/results/merA_ddCPR_CI.svg", plot=image, units="mm", width=250, height=200)

#seems that there is no effect of either Hg concentration or deposition date on merA copy number

#do merA copy numbers in the top 5 cm per site correlate with averaged THg
merA_5cm_avg_df <- merA_ddPCR_data_alt[which(merA_ddPCR_data_alt$MIDDLE.DEPTH <= 5.0),]
merA_5cm_avg_df <- aggregate(.~ Site, merA_5cm_avg_df, mean)
merA_5cm_avg_df$Country <- ifelse(merA_5cm_avg_df$Site %in% c("Pocket Lake", "Aquatuk Lake", "Lake Hazen"), "Canada", "Finland")
merA_5cm_model <- lm(log10(Concentration) ~ CONC..ug.kg., data=merA_5cm_avg_df)
summary(merA_5cm_model)
#no effect of THg across sites

merA_5cm_avg_df$Site <- factor(merA_5cm_avg_df$Site, levels=lake_levels)

lake_colors2 <- c("#f82e55", #Kokemäenjoki
                 "#bb6100", #Pocket Lake
                 "#9ca500", #Öjanjärvi
                 "#194dad", #Päiväjärvi
                 "#ff67a4", #Aquatuk Lake
                 "#a9b3ff", #Pulmankijärvi
                 "#01a981", #Kevojärvi
                 "#0084dc", #Lake Hazen
                 "#aad191") #Vuolimus Cieskuljavri

image <- ggplot(merA_5cm_avg_df, aes(y=Concentration, x=CONC..ug.kg.)) + geom_point(aes(colour = Site, shape = Country), size=5) +
  theme(plot.title = element_text(size=20, hjust = 0.5)) + scale_color_manual(values=lake_colors2, name = "Site") +
  scale_shape_manual(values=c(18,20), name = "Country") + scale_y_log10(labels = scales::comma) + 
  ggtitle("Correlation of merA copy numbers in\ntop 5 cm of sediments with [THg]") + xlab("[THg] ng / g") + ylab("merA copy number") +
  stat_smooth(method = "lm", col = "blue") + labs(subtitle = paste("Adj R2 = ",signif(summary(merA_5cm_model)$adj.r.squared, 5),
                                                                    "Intercept =",signif(merA_5cm_model$coef[[1]],5 ),
                                                                    " Slope =",signif(merA_5cm_model$coef[[2]], 5),
                                                                    " P =",signif(summary(merA_5cm_model)$coef[2,4], 5)))
ggsave(file = "E:/merA_2018/results/merA_5cm.svg", plot=image, units="mm", width=250, height=200)


#pocket lake dates cross-calibration
old_core_data <- read.csv("E:/merA_2018/Pocket_lake_core_2014_THg_dating_data.csv", sep = "\t")
old_pocket_data <- old_core_data[,1:2]
old_pocket_data$Data <- "2014 core"
new_core_data <- hgdata[hgdata$Site %in% "Pocket Lake",c(9,6)]
new_core_data$Data <- "2016 core\n(this study)"
colnames(new_core_data) <- colnames(old_pocket_data)
new_core_data <- new_core_data[order(new_core_data$depth.cm),]
#interpolate the old data to match depth values of the new data
old_pocket_interpolation <- approx(old_pocket_data$depth.cm, old_pocket_data$tHg.ngg, xout=new_core_data$depth.cm)
old_pocket_data <- data.frame(depth.cm = old_pocket_interpolation$x, tHg.ngg = old_pocket_interpolation$y, Data = "2014 core")
combined_pocket_data <- rbind(old_pocket_data, new_core_data)

image <- ggplot(combined_pocket_data, aes(x=depth.cm, y=tHg.ngg, colour=Data)) + geom_line(size=1) + 
  scale_x_reverse(breaks=seq(0,90,5)) + coord_flip() + ylab("[THg] ng / g") + xlab("Sediment depth")
ggsave(file = "E:/merA_2018/results/Pocket_Lake_THg_profiles.svg", plot=image, units="mm", width=250, height=200)

ks.test(old_pocket_data$tHg.ngg, new_core_data$tHg.ngg, exact = T)


# code from Koen Hufkens for plotting field sites on an Arctic map
# https://khufkens.com/2017/01/18/r-polar-plots/
# function to slice and dice a map and convert it to an sp() object
require(raster)
maps2sp <-  function(xlim, ylim, l.out = 100, clip = TRUE) {
  stopifnot(require(maps))
  m <- map(xlim = xlim, ylim = ylim, plot = FALSE, fill = TRUE)
  p <- rbind(cbind(xlim[1], seq(ylim[1],ylim[2],length.out = l.out)),
            cbind(seq(xlim[1],xlim[2],length.out = l.out),ylim[2]),
            cbind(xlim[2],seq(ylim[2],ylim[1],length.out = l.out)),
            cbind(seq(xlim[2],xlim[1],length.out = l.out),ylim[1]))
  LL <- CRS("+init=epsg:4326")
  IDs <- sapply(strsplit(m$names, ":"), function(x) x[1])
  stopifnot(require(maptools))
  m <- map2SpatialPolygons(m, IDs=IDs, proj4string = LL)
  bb <- SpatialPolygons(list(Polygons(list(Polygon(list(p))),"bb")), proj4string = LL)
  
  if (!clip)
    m
  else {
    stopifnot(require(rgeos))
    gIntersection(m, bb)
  }
}

# set colours for map grid
grid.col.light <- rgb(0.5,0.5,0.5,0.8)
grid.col.dark <- rgb(0.5,0.5,0.5)

# coordinate systems
polar <- CRS("+init=epsg:3995")
longlat <- CRS("+init=epsg:4326")

# main map (blue marble next generation, 2004 may) from https://eoimages.gsfc.nasa.gov/images/imagerecords/73000/73701/world.topo.bathy.200405.3x21600x10800.png
# converted to TIFF

# read in the raster map and
# set the extent, crop to extent and reproject to polar
r <- raster::brick("E:/merA_2018/world.topo.bathy.200405.3x21600x10800.TIFF", crs=longlat)
r <- raster::setExtent(r, raster::extent(c(-180,180,-90,90)))
e <- raster::extent(c(-180,180,53,90))
r_crop <- raster::crop(r,e)

# traps NA values and sets them to 1
r_crop[is.na(r_crop)] <- 1 

r_polar <- raster::projectRaster(r_crop, crs = polar, method = "bilinear")

# some values are not valid after transformation 
# (rgb range = 1 - 255) set these back to 1
# as they seem to be the black areas
r_polar[r_polar < 1 ] <- 1

# define the graticule / grid lines by first specifying
# the larger bounding box in which to place them, and
# feeding this into the sp() gridlines function
# finally the grid lines are transformed to
# the EPSG 3995 projection
pts <- SpatialPoints(rbind(c(-180,53),c(0,53),c(180,85),c(180,85)), CRS("+init=epsg:4326"))
gl <- gridlines(pts, easts = seq(-180,180,30), norths = seq(50,85,10), ndiscr = 100)
gl.polar <- spTransform(gl, polar)

# I also create a single line which I use to mark the
# edge of the image (which is rather unclean due to pixelation)
# this line sits at 55 degrees North similar to where I trimmed
# the image
pts <- SpatialPoints(rbind(c(-180,53),c(0,53),c(180,80),c(180,80)), CRS("+init=epsg:4326"))
my_line <- SpatialLines(list(Lines(Line(cbind(seq(-180,180,0.5),rep(53,721))), ID="outer")), CRS("+init=epsg:4326"))

# crop a map object (make the x component a bit larger not to exclude)
# some of the eastern islands (the centroid defines the bounding box)
# and will artificially cut of these islands
m <- maps2sp(c(-180,200),c(53,90),clip = TRUE)

#----- below this point is the plotting routine
# set margins to let the figure "breath" and accommodate labels
par(mar=rep(1,4))

# plot the grid, to initiate the area
# plotRGB() overrides margin settings in default plotting mode
#svg("E:/merA_2018/results/sampling_large.svg", width=10, height=10)
plot(spTransform(gl, polar), lwd=2, lty=2,col="white")

# plot the blue marble raster data
raster::plotRGB(r_polar, add = TRUE)
detach("package:raster", unload=TRUE)

# plot grid lines / graticule
lines(gl.polar, add = TRUE, lwd=2, lty=2,col=grid.col.light)

# plot outer margin of the greater circle
lines(spTransform(my_line, polar), lwd = 3, lty = 1, col=grid.col.dark)

# plot continent outlines, for clarity
plot(spTransform(m, polar), lwd = 1, lty = 1, col = "transparent", border=grid.col.dark, add = T)

# plot longitude labels
crs.longlat <- CRS("+init=epsg:4326")
l <- labels(gl.polar, crs.longlat, side = 1)
l$pos <- NULL
text(l, cex = 1, adj = c( 0.5, 2 ),  col = "black")

# plot latitude labels
l <- labels(gl.polar, crs.longlat, side = 2)
l$srt <- 0
l$pos <- NULL
text(l, cex = 1, adj = c(1.2, -1), col = "white")

# After all this you can plot your own site locations etc
# but don't forget to tranform the data from lat / long
# into the arctic polar stereographic projection using
# spTransform()
site_palette <- c(
  "#f82e55", #Kokemäenjoki
  "#bb6100", #Pocket Lake
  "#9ca500", #Öjanjärvi
  "#194dad", #Päiväjärvi
  "#a9b3ff", #Pulmankijärvi
  "#01a981", #Kevojärvi
  "#0084dc", #Lake Hazen
  "#aad191", #Vuolimus Cieskuljavri
  "#ff67a4") #Aquatuk Lake

sample_coords <- read.csv("E:/merA_2018/sampling_site_coordinates.csv", sep=",")

sample_points <- SpatialPointsDataFrame(sample_coords[,c("lon", "lat")], data.frame(sample_coords[,1]))
sample_points <- SpatialPoints(sample_points, CRS("+init=epsg:4326"))
sample_points$Site <- sample_coords$Site
sample_points.polar <- spTransform(sample_points, polar)
plot(sample_points.polar, cex = 1.5,  col = "black", bg=site_palette, pch=23, add = T)
set.seed(42)
pointLabel(coordinates(sample_points.polar), labels=sample_points.polar$Site, method="GA", col = site_palette, cex=1.1, font=2)
dev.off()

sample_coords_fin <- sample_coords[-which(sample_coords$Site %in% c("POK", "LAH", "AQT")),]
sample_points_fin <- SpatialPointsDataFrame(sample_coords_fin[,c("lon", "lat")], data.frame(sample_coords_fin[,1]))
sample_points_fin <- SpatialPoints(sample_points_fin, CRS("+init=epsg:4326"))
merc <- CRS("+init=epsg:3857")
WGS84 <- CRS("+init=epsg:4326")
sample_points_fin.ll <- spTransform(sample_points_fin, WGS84)
r2 <- raster::brick("E:/merA_2018/world.topo.bathy.200405.3x21600x21600.C1.TIFF", crs=longlat)
r2 <- raster::setExtent(r2, raster::extent(c(0,90,0,90)))
e2 <- raster::extent(c(19,32.5,59.5,70.5))
r2_crop <- raster::crop(r2,e2)
pts2 <- SpatialPoints(rbind(c(19,59.5),c(19,70.5),c(32.5,59.5),c(32.5,70.5)), CRS("+init=epsg:4326"))
gl2 <- gridlines(pts2, easts = seq(20,30,5), norths = seq(60,70,2.5), ndiscr = 100)
cshp <- cshp(as.Date("2000-01-1"))
finland <- cshp[cshp$ISO1AL2 == "FI",]

laea <- CRS("+proj=laea +lat_0=65 +lon_0=26")
r2_laea <- raster::projectRaster(r2_crop, crs = laea, method = "bilinear")

par(mar=rep(1,4))
#svg("E:/merA_2018/results/sampling_finland.svg", width=10, height=10)
plot(spTransform(gl2, laea), lwd=2, lty=2,col="white")
raster::plotRGB(r2_laea, add = TRUE)
detach("package:raster", unload=TRUE)
lines(spTransform(gl2, laea), lwd=2.5, lty=2, col="black")

plot(spTransform(finland, laea), lwd = 3, lty = 1, col = "transparent", border="black", add = T)
text(labels(spTransform(gl2, laea), crs.longlat), cex=1.5)
plot(spTransform(sample_points_fin.ll, laea), cex = 2.5,  col = "black", bg=site_palette[-c(2,7,9)], pch=23, add = T)
dev.off()
