
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}

pstax2veg <- function(physeq) {
  sd <- tax_table(physeq)
  return(as(sd,"matrix"))
}

# convert the otu_table() within a phyloseq object to a vegan compatible data object
psotu2veg <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}


format_physeq_barplot <- function(physeq, taxrank) {
  phyglom <- tax_glom(physeq, taxrank, NArm = TRUE)
  sample_data(phyglom)$samgroup <- paste0(sample_data(phyglom)$Group, "_", sample_data(phyglom)$Timepoint)
  phyglom2 <- merge_samples(phyglom, "samgroup")
  sample_data(phyglom2)$Group <- factor(sample_data(phyglom2)$Group, labels = c("CTRL", "RIF", "RM"))
  sample_data(phyglom2)$Timepoint <- factor(sample_data(phyglom2)$TimepointNumeric, 
                                            labels = c("Follicular", "Ovulation", "Luteal"))
  
  phyglom3 <- transform_sample_counts(phyglom2, function(x) x/sum(x) * 100)
  y4 <- psmelt(phyglom3) # create dataframe from phyloseq object
  
  
  if (taxrank =="Phylum") {
    y4$Phylum <- as.character(y4$Phylum) #convert to character
    y4$Phylum[y4$Abundance < 3] <- "Other" #rename genera with < 5% abundance
  }
  if (taxrank =="Class") {
    y4$Class <- as.character(y4$Class) #convert to character
    y4$Class[y4$Abundance < 3] <- "Other" #rename genera with < 5% abundance
  }
  if (taxrank =="Order") {
    y4$Order <- as.character(y4$Order) #convert to character
    y4$Order[y4$Abundance < 5] <- "Other" #rename genera with < 5% abundance
  }
  
  
  if (taxrank =="Family") {
    y4$Family <- as.character(y4$Family) #convert to character
    y4$Family[y4$Abundance < 5] <- "Other" #rename genera with < 5% abundance
  }
  
  if (taxrank =="Genus") {
    y4$Genus <- as.character(y4$Genus) #convert to character
    y4$Genus[y4$Abundance < 3] <- "Other" #rename genera with < 5% abundance
  }
  y4
}

####---- addMarginalBp() Function for adding marginal boxplots to beta-diversity NMDS plots ----####
addMarginalBp <- function(bdiv_plot) {
  xbox <- axis_canvas(bdiv_plot, axis = "x", coord_flip = TRUE)+
    geom_boxplot(data = bdiv_plot$data, outlier.shape = NA, aes(y = Axis.1, x = factor(Pastegroup), color = factor(Timepoint))) +
    scale_color_manual(values = groupscale)+ 
    scale_x_discrete() + coord_flip()
  
  ybox <- axis_canvas(bdiv_plot, axis = "y") + 
    geom_boxplot(data = bdiv_plot$data, outlier.shape = NA, aes(y = Axis.2, x = factor(Pastegroup), color = factor(Timepoint))) +
    scale_color_manual(values = groupscale)+ 
    scale_x_discrete() 

  
  p1 <- insert_xaxis_grob(bdiv_plot, xbox, grid::unit(1.2, "in"), position = "top")
  p2 <- insert_yaxis_grob(p1, ybox, grid::unit(1.2, "in"), position = "right")
  p2
}

# Function for geometric mean, because of zero counts
geo_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}


# function to generate a Deseq2 data object from a phyloseq object
buildDds <- function(phyloseqobj, taxlevel) {
  
  sample_data(phyloseqobj)$nested <- factor(paste0(sample_data(phyloseqobj)$Group,
                                                   "_",
                                                   sample_data(phyloseqobj)$TimepointNumeric))
  if(taxlevel != "ASV") {
    phyloseqobj <- tax_glom(phyloseqobj, taxrank = taxlevel, NArm = TRUE)
  }
  
  dds = phyloseq_to_deseq2(phyloseqobj, ~nested)
  geoMeans = apply(counts(dds), 1, geo_mean)
  dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  dds$nested <- relevel(dds$nested, ref = "CTRL_1")
  dds = DESeq(dds,fitType = "local") # BH-corrected pvals
  dds
}



# define a function which returns dds results for group and timepoint comparisons
getddsResults <- function(deseqobj, psobj, comparison, alph = 0.05) {
  dummysigtab <- as.data.frame(cbind(NA, 0, 0, 0, 0, "NsPhyl", "NsClass", "NsOrder", "NsFamily", "NsGenus", "NsSpecies", "exclude", NA, "exclude", NA, "exclude"))
  names(dummysigtab) <- c( "ASV", "baseMean","log2FoldChange",
                           "lfcSE", "padj", "Phylum",
                           "Class", "Order", "Family", "Genus", "Species", "Group1",
                           "Time1", "Group2", "Time2",  "comptype")
  res = results(deseqobj, contrast = comparison)
  res = res[order(res$padj, na.last=NA), ]
  sigtab = res[(res$padj < alph), ]
  nsig <- dim(sigtab)[1]
  if(nsig >= 1) {
    sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(psobj)[rownames(sigtab), ], "matrix"))
    sigtab <- tibble::rownames_to_column(sigtab, var = "ASV") %>%
      select(ASV, baseMean, log2FoldChange, lfcSE, padj, Phylum:Species) %>%
      mutate(comp1 = comparison[2], comp2 = comparison[3]) %>%
      separate(comp1, into = c("Group1", "Time1"), sep = "_")%>%
      separate(comp2, into = c("Group2", "Time2"), sep = "_")%>%
      mutate(comptype = ifelse(Group1==Group2, Group1, Time1))
    sigtab
  } else {dummysigtab}
  
}

# define a function which returns formatted results for group and timepoint comparisons within and between groups
formatddsResults <- function(dds, psobj, taxlevel){
  wd <- getwd()
  #dir.create(taxlevel)
  setwd(dir = taxlevel)
  if (taxlevel != "ASV") {psobj <- tax_glom(psobj, taxrank = taxlevel)}
  # Results within groups over time
  res_ctrl2v1 <- getddsResults(dds,psobj,c("nested", "CTRL_2", "CTRL_1"))
  res_ctrl3v1 <- getddsResults(dds,psobj,c("nested", "CTRL_3", "CTRL_1"))
  res_ctrl3v2 <- getddsResults(dds,psobj,c("nested", "CTRL_3", "CTRL_2"))
  res_rm2v1 <- getddsResults(dds,psobj,c("nested", "RM_2", "RM_1"))
  res_rm3v1 <- getddsResults(dds,psobj,c("nested", "RM_3", "RM_1"))
  res_rm3v2 <- getddsResults(dds,psobj,c("nested", "RM_3", "RM_2"))
  res_rif2v1 <- getddsResults(dds,psobj,c("nested", "RIF_2", "RIF_1"))
  res_rif3v1 <- getddsResults(dds,psobj,c("nested", "RIF_3", "RIF_1"))
  res_rif3v2 <- getddsResults(dds,psobj,c("nested", "RIF_3", "RIF_2"))
  
  # Results beween groups within cycle timepoints
  res_RM1 <- getddsResults(dds,psobj,c("nested", "RM_1", "CTRL_1"))
  res_RIF1 <- getddsResults(dds,psobj,c("nested", "RIF_1", "CTRL_1"))
  res_RM2 <- getddsResults(dds,psobj,c("nested", "RM_2", "CTRL_2"))
  res_RIF2 <- getddsResults(dds,psobj,c("nested", "RIF_2", "CTRL_2"))
  res_RM3 <- getddsResults(dds,psobj,c("nested", "RM_3", "CTRL_3"))
  res_RIF3 <- getddsResults(dds,psobj,c("nested", "RIF_3", "CTRL_3"))
  
  # Merge tables
  # within groups table
  res_withingroups <- as.data.frame(rbind(res_ctrl2v1,res_ctrl3v1,res_ctrl3v2,
                                          res_rm2v1,res_rm3v1,res_rm3v2,
                                          res_rif2v1, res_rif3v1, res_rif3v2)) %>%
    select(-Group1, -Group2)
  
  
  res_betweengroups <- as.data.frame(rbind(res_RM1, res_RM2, res_RM3, 
                                           res_RIF1, res_RIF2, res_RIF3)) %>%
    select(-Time1,-Time2)
  
  rm(res_ctrl2v1,res_ctrl3v1,res_ctrl3v2,
     res_rm2v1,res_rm3v1,res_rm3v2,
     res_rif2v1, res_rif3v1, res_rif3v2,
     res_RM1, res_RM2, res_RM3, 
     res_RIF1, res_RIF2, res_RIF3)
  
  res_withingroups$padj <- as.numeric(res_withingroups$padj)
  res_betweengroups$padj <- as.numeric(res_betweengroups$padj)
  res_withingroups <- filter(res_withingroups, !is.na(ASV))
  res_betweengroups <- filter(res_betweengroups, !is.na(ASV))
  
  #Write results to file
  write.table(x = res_withingroups, 
              file = "res_withingroups.csv",
              sep = ";", row.names = F)
  write.table(x = res_betweengroups, 
              file = "res_betweengroups.csv",
              sep = ";", row.names = F)
  
  # Make wide format
  
  wide_withingroups <- res_withingroups %>%
    group_by(comptype, ASV) %>%
    mutate(Time = paste0(Time1, "_", Time2)) %>%
    select(-Time1:-Time2, -baseMean:-lfcSE, Group = comptype) %>%
    tidyr::pivot_wider(names_from = Time, names_prefix = "comp", values_from = padj, values_fill = 1)
  
  wide_betweengroups <- res_betweengroups %>%
    group_by(comptype, ASV) %>%
    mutate(Group = paste0(Group1, "_", Group2)) %>%
    select(-Group1:-Group2, -baseMean:-lfcSE, Timepoint = comptype) %>%
    tidyr::pivot_wider(names_from = Group, names_prefix = "comp", values_from = padj, values_fill = 1)
  
  
  write.table(x = wide_withingroups, 
              file = "res_withingroups_pvals.csv",
              sep = ";", row.names = F)
  write.table(x = wide_betweengroups, 
              file = "res_betweengroups_pvals.csv",
              sep = ";", row.names = F)
  setwd(wd)
}
