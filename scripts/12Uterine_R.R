library(dplyr)
library(tibble)
library(tidyr)
library(phyloseq)
library(biomformat)
library(ggplot2)
library(ggpubr)
library(vegan)
library(DESeq2)
library(cowplot)
library(RColorBrewer)
library(scales)
library(car)


wd <- "C:/Users/Simon Reider/Google Drive/UTERINE/upload github/"

setwd(wd)
source(file = "scripts/12Uterine_analysis_functions.R")

sv <- read_biom("qiime2/table_filt.biom")
sv2 <- biom_data(sv)
meta <- read.table("metadata/uterine_metadata.txt", sep = "\t", header=TRUE)
rownames(meta) <- meta$Sample.id
meta <- dplyr::select(meta,-Sample.id )
meta$Group <- factor(meta$Group, labels = c("CTRL", "RIF", "RM"))
m <- phyloseq::sample_data(meta)
taxo <- read.delim("qiime2/taxonomy.tsv")
taxo2 <- taxo %>% 
  tidyr::separate("Taxon", into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";")

rownames(taxo2) <- taxo2$Feature.ID
taxo2 <- dplyr::select(taxo2, -Feature.ID)
tax <- tax_table(as(taxo2, "matrix"))
tree <- 
ps_input <- phyloseq(otu_table(as.matrix(sv2), taxa_are_rows = TRUE), m, tax)
ps_input <- subset_samples(ps_input, Patient !="CTRLc011")

rm(m,meta,sv,sv2,taxo,taxo2,tax)


# get percentages unassigned per taxlevel
taxontable <- as.data.frame(as(tax_table(ps_input), "matrix"))
kingdom <- summary(factor(taxontable$Kingdom))
phylum <- summary(factor(taxontable$Phylum))
class <- summary(factor(taxontable$Class))
order <- summary(factor(taxontable$Order))
family <- summary(factor(taxontable$Family))
genus <- summary(factor(taxontable$Genus))

ps <- subset_taxa(ps_input, Class!=" c__Chloroplast") # remove Chloroplasts
ps <- subset_taxa(ps, Family!=" f__mitochondria") # remove Mitochondria
ntaxa(ps_input)-ntaxa(ps) # 976 taxa are removed


taxontable <- as.data.frame(as(tax_table(ps), "matrix"))
kingdom <- summary(factor(taxontable$Kingdom))
phylum <- summary(factor(taxontable$Phylum))
class <- summary(factor(taxontable$Class))
order <- summary(factor(taxontable$Order))
family <- summary(factor(taxontable$Family))
genus <- summary(factor(taxontable$Genus))

## rarefication for beta-diversity
set.seed(111)
rarecurve(t(otu_table(ps)), step=250, label = F, xlim = c(0,5000)) # based on this, rarefaction @ 2000 seems sufficient

ps.ns <- phyloseq::prune_taxa(taxa_sums(ps) > 1, ps)
ps.ns.rel <- transform_sample_counts(ps.ns, function(x) x / sum(x)*100)

ps.rare <- phyloseq::rarefy_even_depth(ps.ns, sample.size = 2000)
ps.filt <- phyloseq::filter_taxa(ps.ns, function(x) sum(x >= 3) >= (0.05*length(x)), TRUE)
ps.filt.rel <- transform_sample_counts(ps.filt, function(x) x / sum(x)*100)


# a diversity

adat4 <- phyloseq::plot_richness(ps, x = "Timepoint", measures=c("Chao1", "Shannon"))
adiv_time <- adat4$data %>%
  group_by(variable, Group, Patient) %>%
  select(Patient,Group, Timepoint, variable, value) %>%
  tidyr::pivot_wider(names_from = Timepoint, values_from = value) %>%
  mutate(diff12 = `1_follicular`-`2_ovul`, diff23 = `2_ovul`-`3_luteal`) %>%
  group_by(variable, Group) 

adiv_time2 <- adiv_time %>%
  mutate(ovulation = `2_ovul` - `1_follicular`,
         luteal = `3_luteal` - `1_follicular`) %>%
  group_by() %>%
  select(Patient,Group, variable, ovulation, luteal) %>%
  mutate(follicular = 0) %>%
  pivot_longer(cols = ovulation:follicular, names_to = "timepoint")

adiv_time2$timepoint = factor(adiv_time2$timepoint, labels = c("follicular", "ovulation", "luteal"))


adiv_time_summary <- adiv_time %>%
  select(Group, variable, diff12, diff23) %>%
  summarize_all(.funs = c("mean", "sd"))

adiv_time_summary2 <- filter(adiv_time, Patient !="CTRLc011") %>%
  select(Group, variable, diff12, diff23) %>%
  summarize_all(.funs = c("mean", "sd"))


wilcox_comparisons <- list( c("follicular", "ovulation"), c("ovulation", "luteal"), c("follicular", "luteal") )


p_overtime <- ggplot(adiv_time2,
                     aes(x = timepoint, y = value))+
  geom_point()+
  geom_line(aes(group = Patient))+
  facet_grid(variable~Group, scales = "free")+
  labs(x = "Menstrual cycle phase", y = "Alpha diversity measure")

p_overtime2 <- p_overtime
p_overtime2$data <- filter(p_overtime2$data, Patient != "CTRLc011")
p_overtime2shannon <- p_overtime2 
p_overtime2shannon$data <- filter(p_overtime2shannon$data, variable =="Shannon") 

stat.test_shannon <- p_overtime2shannon$data %>%
  group_by(Group) %>%
  rstatix::wilcox_test(value~timepoint,
                       comparisons= wilcox_comparisons,
                       p.adjust.method = "BH")
stat.test_shannon

shanplotfin <-p_overtime2shannon+
  theme_bw(base_size = 7)+
  lims(y = c(-5,5))+
  stat_compare_means(label.y = -5, size= 3,label = "p.format")+
  stat_pvalue_manual(stat.test_shannon, y.position = c(3,3.5,4),label = "p.adj.signif", tip.length = 0.01)


p_overtime2chao <- p_overtime2
p_overtime2chao$data <- filter(p_overtime2chao$data, variable =="Chao1") 

stat.test_chao <- p_overtime2chao$data %>%
  group_by(Group) %>%
  rstatix::wilcox_test(value~timepoint,
                       comparisons= wilcox_comparisons,
                       p.adjust.method = "BH")
stat.test_chao

chaoplotfin<-p_overtime2chao+
  theme_bw(base_size = 7)+
  lims(y = c(-150,150))+
  stat_compare_means(label.y = -150, size = 3, hide.ns = F, label = "p.format")+
  stat_pvalue_manual(stat.test_chao, y.position = c(100, 120, 140),label = "p.adj.signif",tip.length = 0.01)

#Levene Test for equal variances in CTRL, RIF, RM Group regardless of timepoint
levene_chao1 <- leveneTest(value~Group, data =chaoplotfin$data)
levene_shannon <- leveneTest(value~Group, data =shanplotfin$data)

# Figure 1B = Alpha diversity plots
alphaplot <- gridExtra::grid.arrange(chaoplotfin, shanplotfin, ncol =1)
ggsave(dpi = "print", filename = "reporting/11-2020/allpat_alphaovertime_withsign.pdf", plot = alphaplot, device = "pdf", width = 5, height = 5)
ggsave(dpi = "print",filename = "reporting/11-2020/allpat_alphaovertime_withsign.png", plot = alphaplot, device = "png", width = 5, height = 5)

# Export Raw data
# Paired, i.e. follicular = 0
write.table(chaoplotfin$data, file = "reporting/11-2020/adiv-rawdata_paired.csv", sep = ";", row.names = F)
write.table(shanplotfin$data, file = "reporting/11-2020/adiv-rawdata_paired.csv", sep = ";", row.names = F, col.names = F,
            append = TRUE)

# unpaired data
write.table(select(adat4$data, samples,Patient, Group, Timepoint, variable, value, se),
                   file = "reporting/11-2020/adiv-rawdata_unpaired.csv", sep = ";", row.names = F)

# Make statistics summary for Chao and Shannon per grozp and timepoint
chaostat = chaoplotfin$data %>%
  group_by(Group, timepoint) %>%
  summarise(Chao1_mean = mean(value),
            Chao1_sd = sd(value),
            Chao1_median = median(value),
            Chao1_IQR = IQR(value))

shannonstat = shanplotfin$data %>%
  group_by(Group, timepoint) %>%
  summarise(Shannon_mean = mean(value),
            Shannon_sd = sd(value),
            Shannon_median = median(value),
            Shannon_IQR = IQR(value))

astats <- inner_join(chaostat, shannonstat, by = c("Group", "timepoint"))

write.table(astats, file = "reporting/11-2020/allpat_alphaovertimestat_paired.csv")

# Taxonomy  barplots on Genus, Family, and Phylum levels for all3 Cohorts over time with connecting arrows
# for every taxlevel: use max 12 most abundant, group the rest.


mycolors <- colorRampPalette(brewer.pal(8, "Dark2"))(12)

phylumplotdat <- format_physeq_barplot(ps.ns, "Phylum")
phylumplotdat <- phylumplotdat %>%
  mutate(Phylum = factor(Phylum), 
         Phylum = forcats::fct_reorder(factor(Phylum), .x = Abundance, .fun = mean))

phylumplot <- ggplot(phylumplotdat, aes(x = Timepoint, y = Abundance, color = Phylum, fill = Phylum))+
  geom_bar(stat="identity", position = "stack", width = 0.7)+
  facet_wrap(~Group)+
  scale_color_manual(values =mycolors)+
  scale_fill_manual(values =mycolors)+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90))

classplotdat <- format_physeq_barplot(ps.ns, "Class")
classplotdat <- classplotdat %>%
  mutate(Class = factor(Class), 
         Class = forcats::fct_reorder(factor(Class), .x = Abundance, .fun = mean))

classplot <- ggplot(classplotdat, aes(x = Timepoint, y = Abundance, color = Class, fill = Class))+
  geom_bar(stat="identity", position = "stack", width = 0.7)+
  facet_wrap(~Group)+
  scale_color_manual(values =mycolors)+
  scale_fill_manual(values =mycolors)+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90))


orderplotdat <- format_physeq_barplot(ps.ns, "Order")
orderplotdat <- orderplotdat %>%
  mutate(Order = factor(Order), 
         Order = forcats::fct_reorder(factor(Order), .x = Abundance, .fun = mean))

orderplot <- ggplot(orderplotdat, aes(x = Timepoint, y = Abundance, color = Order, fill = Order))+
  geom_bar(stat="identity", position = "stack", width = 0.7)+
  facet_wrap(~Group)+
  scale_color_manual(values =mycolors)+
  scale_fill_manual(values =mycolors)+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90))

familyplotdat <- format_physeq_barplot(ps.ns, "Family")
familyplotdat <- familyplotdat %>%
  mutate(Family = factor(Family), 
         Family = forcats::fct_reorder(factor(Family), .x = Abundance, .fun = mean))

familyplot <- ggplot(familyplotdat, aes(x = Timepoint, y = Abundance, color = Family, fill = Family))+
  geom_bar(stat="identity", position = "stack", width = 0.7)+
  facet_wrap(~Group)+
  scale_color_manual(values =mycolors)+
  scale_fill_manual(values =mycolors)+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90))



genusplotdat <- format_physeq_barplot(ps.ns, "Genus")
genusplotdat <- genusplotdat %>%
  mutate(Genus = factor(Genus), 
         Genus = forcats::fct_reorder(factor(Genus), .x = Abundance, .fun = mean)) 


genusplot <- ggplot(genusplotdat, aes(x = Timepoint, y = Abundance, color = Genus, fill = Genus))+
  geom_bar(stat="identity", position = "stack", width = 0.7)+
  facet_wrap(~Group)+
  scale_color_manual(values =mycolors)+
  scale_fill_manual(values =mycolors)+
  theme_bw(base_size = 10)+
  theme(axis.text.x = element_text(angle = 90))


ggsave("reporting/11-2020/12_taxaplot-genus.png", genusplot, width = 5, height = 4)
ggsave("reporting/11-2020/12_taxaplot-family.png", familyplot, width = 5, height = 4)
ggsave("reporting/11-2020/12_taxaplot-order.png", orderplot, width = 5, height = 4)
ggsave("reporting/11-2020/12_taxaplot-class.png", classplot, width = 5, height = 4)
ggsave("reporting/11-2020/12_taxaplot-phylum.png", phylumplot, width = 5, height = 4)

ggsave("reporting/11-2020/12_taxaplot-genus.pdf", genusplot, width = 5, height = 4)
ggsave("reporting/11-2020/12_taxaplot-order.pdf", orderplot, width = 5, height = 4)
ggsave("reporting/11-2020/12_taxaplot-class.pdf", classplot, width = 5, height = 4)
ggsave("reporting/11-2020/12_taxaplot-family.pdf", familyplot, width = 5, height = 4)
ggsave("reporting/11-2020/12_taxaplot-phylum.pdf", phylumplot, width = 5, height = 4)

taxaplot <- plot_grid(
  plot_grid(
    phylumplot + theme(legend.position = "none"),
    classplot+ theme(legend.position = "none"),
    orderplot+ theme(legend.position = "none"),
    familyplot+ theme(legend.position = "none"),
    genusplot + theme(legend.position = "none"),
    ncol = 1,
    align = "hv"),
  plot_grid(
    get_legend(phylumplot),
    get_legend(classplot),
    get_legend(orderplot),
    get_legend(familyplot),
    get_legend(genusplot), 
    ncol =1, 
    align = "hv")
  , rel_widths = c(7,3)
  )
ggsave("reporting/11-2020/12_taxaplot-all.pdf", taxaplot, width = 5.5, height = 18)


x_phylum <- psmelt(tax_glom(ps.ns.rel, "Phylum"))
x_class <- psmelt(tax_glom(ps.ns.rel, "Class"))
x_order <- psmelt(tax_glom(ps.ns.rel, "Order"))
x_family <- psmelt(tax_glom(ps.ns.rel, "Family"))
x_genus <- psmelt(tax_glom(ps.ns.rel, "Genus"))

phylumstat <- x_phylum %>%
  group_by(Phylum, Timepoint, Group) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)


phylumstat2 <- psmelt(tax_glom(ps.filt.rel, taxrank ="Phylum")) %>%
  group_by(Phylum, Timepoint, Group) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)

classstat2 <- psmelt(tax_glom(ps.filt.rel, taxrank ="Class")) %>%
  group_by(Class, Timepoint, Group) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)

familystat2 <- psmelt(tax_glom(ps.filt.rel, taxrank ="Family")) %>%
  group_by(Family, Timepoint, Group) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)

genusstat2 <- psmelt(tax_glom(ps.filt.rel, taxrank ="Genus")) %>%
  group_by(Genus, Timepoint, Group) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)


classstat <- x_class %>%
  group_by(Class, Timepoint, Group) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)

orderstat <- x_order %>%
  group_by(Order, Timepoint, Group) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)

familystat <- x_family %>%
  group_by(Family, Timepoint, Group) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)

genusstat <- x_genus %>%
  group_by(Genus, Timepoint, Group) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)

phylumstat3 <- x_phylum %>%
  group_by(Phylum, Timepoint, Group, IsolatedKitLot) %>%
  summarize(mean = mean(Abundance, na.rm = T),
            median = median(Abundance, na.rm =T),
            sd = sd(Abundance, na.rm =T),
            n = length(unique(x_phylum$Sample)),
            IQR = IQR(Abundance, na.rm =T)) %>%
  filter(mean > 5)

write.table("reporting/11-2020/12_taxatable_phylum.csv", sep = ";", col.names  = T, row.names=T, x = phylumstat)
write.table("reporting/11-2020/12_taxatable_class.csv", sep = ";", col.names  = T, row.names=T, x = classstat)
write.table("reporting/11-2020/12_taxatable_order.csv", sep = ";", col.names  = T, row.names=T, x = orderstat)
write.table("reporting/11-2020/12_taxatable_family.csv", sep = ";", col.names  = T, row.names=T, x = familystat)
write.table("reporting/11-2020/12_taxatable_genus.csv", sep = ";", col.names  = T, row.names=T, x = genusstat)
write.table("reporting/11-2020/12_taxatable_phylum_filt.csv", sep = ";", col.names  = T, row.names=T, x = phylumstat2)
write.table("reporting/11-2020/12_taxatable_class_filt", sep = ";", col.names  = T, row.names=T, x = classstat2)
write.table("reporting/11-2020/12_taxatable_family_filt", sep = ";", col.names  = T, row.names=T, x = familystat2)
write.table("reporting/11-2020/12_taxatable_genus_filt.csv", sep = ";", col.names  = T, row.names=T, x = genusstat2)


####---- Beta-Diversity Analysis ----####
# Transform data to proportions as appropriate for Bray-Curtis distances
set.seed(111)
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))

sample_data(ps.prop)$Pastegroup <- factor(paste0(sample_data(ps.prop)$Group, "_", sample_data(ps.prop)$Timepoint),
levels = c("CTRL_1_follicular", "CTRL_2_ovul", "CTRL_3_luteal",
"RIF_1_follicular", "RIF_2_ovul", "RIF_3_luteal",
"RM_1_follicular", "RM_2_ovul", "RM_3_luteal"))



ord.nmds.bray <- ordinate(ps.prop, method="MDS", distance="bray")

ord.nmds.wuni <- ordinate(ps.prop, method="MDS", distance="wunifrac")
ord.nmds.uwuni <- ordinate(ps.prop, method="MDS", distance="unifrac")

groupscale <- c("1_follicular" = "#009933", 
                "2_ovul" = "#0000ff", 
                "3_luteal" = "#ff9900")

bray <- plot_ordination(ps.prop, ord.nmds.bray, shape = "Group", color="Timepoint", title="Bray MDS")+
  scale_color_manual(values = groupscale)+
  theme_bw(base_size = 12)+
  theme(legend.position = "right")+
  coord_fixed()+
  geom_point(size = 3)
braybp <- addMarginalBp(bray)

wuni <- plot_ordination(ps.prop, ord.nmds.wuni, shape = "Group", color="Timepoint", title="W. Unifrac MDS")+
  scale_color_manual(values = groupscale)+
  theme_bw(base_size = 12)+
  theme(legend.position = "right")+
  coord_fixed()+
  geom_point(size = 3)
wunibp <- addMarginalBp(wuni)

uwuni <- plot_ordination(ps.prop, ord.nmds.uwuni, shape = "Group", color="Timepoint", title="Unw. Unifrac MDS")+
  scale_color_manual(values = groupscale)+
  theme_bw(base_size = 12)+
  theme(legend.position = "right")+
  coord_fixed()+
  geom_point(size = 3)
uwunibp <- addMarginalBp(uwuni)


ggsave(bray+facet_wrap(~Group), filename = "reporting/11-2020/uterine_bdiv_bray.png", width =8, height = 4)
ggsave(bray+facet_wrap(~Group), filename = "reporting/11-2020/uterine_bdiv_bray.pdf", width =8, height = 4)

ggsave(wuni+facet_wrap(~Group), filename = "reporting/11-2020/uterine_bdiv_wuni.png", width =8, height = 4)
ggsave(wuni+facet_wrap(~Group), filename = "reporting/11-2020/uterine_bdiv_wuni.pdf", width =8, height = 4)

ggsave(uwuni+facet_wrap(~Group), filename = "reporting/11-2020/uterine_bdiv_uwuni.png", width =8, height = 4)
ggsave(uwuni+facet_wrap(~Group), filename = "reporting/11-2020/uterine_bdiv_uwuni.pdf", width =8, height = 4)

ggsave(braybp, filename = "reporting/11-2020/uterine_bdiv_braybp.png", width =6, height = 5)
ggsave(braybp, filename = "reporting/11-2020/uterine_bdiv_braybp.pdf", width =6, height = 5)

ggsave(wunibp, filename = "reporting/11-2020/uterine_bdiv_wunibp.png", width =6, height = 5)
ggsave(wunibp, filename = "reporting/11-2020/uterine_bdiv_wunibp.pdf", width =6, height = 5)

ggsave(uwunibp, filename = "reporting/11-2020/uterine_bdiv_uwunibp.png", width =6, height = 5)
ggsave(uwunibp, filename = "reporting/11-2020/uterine_bdiv_uwunibp.pdf", width =6, height = 5)


# Assess statistical significance using vegan::adonis() implementation of permanova
metadata <- as(sample_data(ps.prop), "data.frame")
permanova_bray <- adonis(phyloseq::distance(ps.prop, method = "bray") ~ Group*Timepoint, data = metadata)
permanova_uwuni <- adonis(phyloseq::distance(ps.prop, method = "unifrac") ~ Group*Timepoint, data = metadata)
permanova_wuni <- adonis(phyloseq::distance(ps.prop, method = "wunifrac") ~ Group*Timepoint, data = metadata)

readr::write_excel_csv(permanova_bray$aov.tab, path = "reporting/11-2020/permanova_bray_results.csv")
readr::write_excel_csv(permanova_uwuni$aov.tab, path = "reporting/11-2020/permanova_unw-unifrac_results.csv")
readr::write_excel_csv(permanova_wuni$aov.tab, path = "reporting/11-2020/permanova_w-unifrac_results.csv")



####---- DESEQ2 analsis of differential abudnacne
# According ot Holmes/McMurdie SOP https://bioconductor.org/packages/devel/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html#convert-to-deseq2s-deseqdataset-class

# 1) Filtering Dataset
#dir.create("reporting/09-2020/deseq2")

setwd("reporting/11-2020/deseq2")


dds_asv <- buildDds(ps.filt, "ASV")
dds_genus <- buildDds(ps.filt, "Genus")
dds_family <- buildDds(ps.filt, "Family")
dds_order <- buildDds(ps.filt, "Order")
dds_class <- buildDds(ps.filt, "Class")
dds_phylum <- buildDds(ps.filt, "Phylum")



formatddsResults(dds_asv, ps.filt, "ASV")
formatddsResults(dds_genus, ps.filt, "Genus")
formatddsResults(dds_family, ps.filt, "Family")
formatddsResults(dds_class, ps.filt, "Class")
formatddsResults(dds_order, ps.filt, "Order")
formatddsResults(dds_phylum, ps.filt, "Phylum")

setwd(wd)



#### ----Plot deseq2 results as heatmap---- ####

jBuPuFun <- colorRampPalette(brewer.pal(n = 9, "RdBu"))
paletteSize <- 256
jBuPuPalette <- jBuPuFun(paletteSize)

####---- Within Groups ----####

withinASV <- read.csv("reporting/11-2020/deseq2/ASV/res_withingroups_pvals.csv", sep = ";") 
withinGenus <- read.csv("reporting/11-2020/deseq2/Genus/res_withingroups_pvals.csv", sep = ";")
withinFamily <- read.csv("reporting/11-2020/deseq2/Family/res_withingroups_pvals.csv", sep = ";")
withinOrder <- read.csv("reporting/11-2020/deseq2/Order/res_withingroups_pvals.csv", sep = ";")
withinClass <- read.csv("reporting/11-2020/deseq2/Class/res_withingroups_pvals.csv", sep = ";")
withinPhylum <- read.csv("reporting/11-2020/deseq2/Phylum/res_withingroups_pvals.csv", sep = ";")

### Phylum level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt()  %>%
  filter(Phylum %in% withinPhylum$Phylum) %>%
  select(Phylum, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Phylum, Sample, Group, Timepoint, Abundance) %>%
  group_by(Phylum, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Phylum, Group, Timepoint, MeanAbund) %>%
  group_by(Phylum, Group) %>%
  tidyr::pivot_wider(names_from = Timepoint, values_from = MeanAbund) %>%
  mutate(ovulation = `2_ovul`/`1_follicular`,
         luteal = `3_luteal`/`1_follicular`,
         follicular = 1) %>%
  select(Phylum, Group, follicular, ovulation, luteal) %>%
  tidyr::pivot_longer(cols = follicular:luteal, names_to = "Timepoint") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("follicular", "ovulation", "luteal")))


withinPhylumplot <- ggplot(ds2, aes(y= Phylum, x = Timepoint, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                      low = jBuPuPalette[paletteSize],
                      mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),oob = squish, 
                      name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Group)

ggsave("reporting/11-2020/deseq2/phylum_within_plot.png", withinPhylumplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/phylum_within_plot.pdf", withinPhylumplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
withinPhylumSig <- left_join(ds2, withinPhylum, by = "Phylum") %>% 
  filter(Group.x == Group.y) %>% 
  select(Phylum, Group = Group.x, Timepoint, value, comp2_1, comp3_1, comp3_2)
write.table(withinPhylumSig, file = "reporting/11-2020/deseq2/phylum_final_within-sig.csv", sep= ";")

### Class level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Class") %>%
  psmelt()  %>%
  filter(Class %in% withinClass$Class) %>%
  select(Class, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Class, Sample, Group, Timepoint, Abundance) %>%
  group_by(Class, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Class, Group, Timepoint, MeanAbund) %>%
  group_by(Class, Group) %>%
  tidyr::pivot_wider(names_from = Timepoint, values_from = MeanAbund) %>%
  mutate(ovulation = `2_ovul`/`1_follicular`,
         luteal = `3_luteal`/`1_follicular`,
         follicular = 1) %>%
  select(Class, Group, follicular, ovulation, luteal) %>%
  tidyr::pivot_longer(cols = follicular:luteal, names_to = "Timepoint") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("follicular", "ovulation", "luteal")))


withinClassplot <- ggplot(ds2, aes(y= Class, x = Timepoint, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),oob = squish, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Group)

ggsave("reporting/11-2020/deseq2/class_within_plot.png", withinClassplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/class_within_plot.pdf", withinClassplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
withinClassSig <- left_join(ds2, withinClass, by = "Class") %>% 
  filter(Group.x == Group.y) %>% 
  select(Class, Group = Group.x, Timepoint, value, comp2_1, comp3_1, comp3_2)
write.table(withinClassSig, file = "reporting/11-2020/deseq2/class_final_within-sig.csv", sep= ";")

rm(ds, ds2, dssum)

### Order level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Order") %>%
  psmelt()  %>%
  filter(Order %in% withinOrder$Order) %>%
  select(Order, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Order, Sample, Group, Timepoint, Abundance) %>%
  group_by(Order, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Order, Group, Timepoint, MeanAbund) %>%
  group_by(Order, Group) %>%
  tidyr::pivot_wider(names_from = Timepoint, values_from = MeanAbund) %>%
  mutate(ovulation = `2_ovul`/`1_follicular`,
         luteal = `3_luteal`/`1_follicular`,
         follicular = 1) %>%
  select(Order, Group, follicular, ovulation, luteal) %>%
  tidyr::pivot_longer(cols = follicular:luteal, names_to = "Timepoint") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("follicular", "ovulation", "luteal")))

withinOrderplot <- ggplot(ds2, aes(y= Order, x = Timepoint, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),oob = squish, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Group)

ggsave("reporting/11-2020/deseq2/order_within_plot.png", withinOrderplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/order_within_plot.pdf", withinOrderplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
withinOrdersig <- left_join(ds2, withinOrder, by = "Order") %>% 
  filter(Group.x == Group.y) %>% 
  select(Order, Group = Group.x, Timepoint, value, comp2_1, comp3_1, comp3_2)
write.table(withinOrdersig, file = "reporting/11-2020/deseq2/order_final_within-sig.csv", sep= ";")

rm(ds, ds2, dssum)

### Family level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Family") %>%
  psmelt()  %>%
  filter(Family %in% withinFamily$Family) %>%
  select(Family, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Family, Sample, Group, Timepoint, Abundance) %>%
  group_by(Family, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Family, Group, Timepoint, MeanAbund) %>%
  group_by(Family, Group) %>%
  tidyr::pivot_wider(names_from = Timepoint, values_from = MeanAbund) %>%
  mutate(ovulation = `2_ovul`/`1_follicular`,
         luteal = `3_luteal`/`1_follicular`,
         follicular = 1) %>%
  select(Family, Group, follicular, ovulation, luteal) %>%
  tidyr::pivot_longer(cols = follicular:luteal, names_to = "Timepoint") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("follicular", "ovulation", "luteal")))

squish_any = scales::oob_squish_any

withinFamilyplot <- ggplot(ds2, aes(y= Family, x = Timepoint, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),
                       oob = squish_any, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Group)

ggsave("reporting/11-2020/deseq2/family_within_plot.png", withinFamilyplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/family_within_plot.pdf", withinFamilyplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
withinFamilySig <- left_join(ds2, withinFamily, by = "Family") %>% 
  filter(Group.x == Group.y) %>% 
  select(Family, Group = Group.x, Timepoint, value, comp2_1, comp3_1, comp3_2)
write.table(withinFamilySig, file = "reporting/11-2020/deseq2/family_final_within-sig.csv", sep= ";")

rm(ds, ds2, dssum)


### Genus level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()  %>%
  filter(Genus %in% withinGenus$Genus) %>%
  select(Genus, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Genus, Sample, Group, Timepoint, Abundance) %>%
  group_by(Genus, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Genus, Group, Timepoint, MeanAbund) %>%
  group_by(Genus, Group) %>%
  tidyr::pivot_wider(names_from = Timepoint, values_from = MeanAbund) %>%
  mutate(ovulation = `2_ovul`/`1_follicular`,
         luteal = `3_luteal`/`1_follicular`,
         follicular = 1) %>%
  select(Genus, Group, follicular, ovulation, luteal) %>%
  tidyr::pivot_longer(cols = follicular:luteal, names_to = "Timepoint") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("follicular", "ovulation", "luteal")))


withinGenusplot <- ggplot(ds2, aes(y= Genus, x = Timepoint, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),
                       oob = squish_any, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Group)

ggsave("reporting/11-2020/deseq2/genus_within_plot.png", withinGenusplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/genus_within_plot.pdf", withinGenusplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
withinGenusSig <- left_join(ds2, withinGenus, by = "Genus") %>% 
  filter(Group.x == Group.y) %>% 
  select(Genus, Group = Group.x, Timepoint, value, comp2_1, comp3_1, comp3_2)
write.table(withinGenusSig, file = "reporting/11-2020/deseq2/genus_final_within-sig.csv", sep= ";")

rm(ds, ds2, dssum)


### ASV level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  psmelt() %>%
  mutate(ASV = paste(Genus, Species, OTU, sep = "|")) %>%
  filter(OTU %in% withinASV$ASV) %>%
  select(ASV, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(ASV, Sample, Group, Timepoint, Abundance) %>%
  group_by(ASV, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(ASV, Group, Timepoint, MeanAbund) %>%
  group_by(ASV, Group) %>%
  tidyr::pivot_wider(names_from = Timepoint, values_from = MeanAbund) %>%
  mutate(ovulation = `2_ovul`/`1_follicular`,
         luteal = `3_luteal`/`1_follicular`,
         follicular = 1) %>%
  select(ASV, Group, follicular, ovulation, luteal) %>%
  tidyr::pivot_longer(cols = follicular:luteal, names_to = "Timepoint") %>%
  mutate(Timepoint = factor(Timepoint, levels = c("follicular", "ovulation", "luteal")))


withinASVplot <- ggplot(ds2, aes(y= ASV, x = Timepoint, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-4,4),
                       oob = squish_any, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Group)

ggsave("reporting/11-2020/deseq2/asv_within_plot.png", withinASVplot, width = 10, height= 6)
ggsave("reporting/11-2020/deseq2/asv_within_plot.pdf", withinASVplot, width = 10, height= 6)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
withinASVsig <- left_join(ds2, withinASV, by = "ASV") %>% 
  filter(Group.x == Group.y) %>% 
  select(ASV, Group = Group.x, Timepoint, value, comp2_1, comp3_1, comp3_2)
write.table(withinASVsig, file = "reporting/11-2020/deseq2/genus_final_within-sig.csv", sep= ";")

rm(ds, ds2, dssum)

####---- Between Groups ----####

betweenASV <- read.csv("reporting/11-2020/deseq2/ASV/res_betweengroups_pvals.csv", sep = ";") 
betweenGenus <- read.csv("reporting/11-2020/deseq2/Genus/res_betweengroups_pvals.csv", sep = ";")
betweenFamily <- read.csv("reporting/11-2020/deseq2/Family/res_betweengroups_pvals.csv", sep = ";")
betweenOrder <- read.csv("reporting/11-2020/deseq2/Order/res_betweengroups_pvals.csv", sep = ";")
betweenClass <- read.csv("reporting/11-2020/deseq2/Class/res_betweengroups_pvals.csv", sep = ";")
betweenPhylum <- read.csv("reporting/11-2020/deseq2/Phylum/res_betweengroups_pvals.csv", sep = ";")

### Phylum level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Phylum") %>%
  psmelt()  %>%
  filter(Phylum %in% betweenPhylum$Phylum) %>%
  select(Phylum, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Phylum, Sample, Group, Timepoint, Abundance) %>%
  group_by(Phylum, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Phylum, Group, Timepoint, MeanAbund) %>%
  group_by(Phylum, Timepoint) %>%
  tidyr::pivot_wider(names_from = Group, values_from = MeanAbund) %>%
  mutate(RIF = RIF/CTRL,
         RM = RM/CTRL,
         CTRL = 1) %>%
  select(Phylum, Timepoint, CTRL, RIF, RM) %>%
  tidyr::pivot_longer(cols = CTRL:RM, names_to = "Group") %>%
  mutate(Group = factor(Group, levels = c("CTRL", "RIF", "RM")))


betweenPhylumplot <- ggplot(ds2, aes(y= Phylum, x = Group, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),oob = squish, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Timepoint)

ggsave("reporting/11-2020/deseq2/phylum_between_plot.png", betweenPhylumplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/phylum_between_plot.pdf", betweenPhylumplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
betweenPhylumSig <- left_join(ds2, betweenPhylum, by = "Phylum") %>%
  mutate(Timepoint.y = factor(Timepoint.y, levels = c(1,2,3), labels = c("1_follicular", "2_ovulation", "3_luteal"))) %>%
  filter(Timepoint.x == Timepoint.y) %>%
  select(Phylum, Timepoint = Timepoint.x, Group, value, compRM_CTRL, compRIF_CTRL)
write.table(betweenPhylumSig, file = "reporting/11-2020/deseq2/phylum_final_between-sig.csv", sep= ";")
rm(ds, ds2, dssum)


### Class level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Class") %>%
  psmelt()  %>%
  filter(Class %in% betweenClass$Class) %>%
  select(Class, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Class, Sample, Group, Timepoint, Abundance) %>%
  group_by(Class, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Class, Group, Timepoint, MeanAbund) %>%
  group_by(Class, Timepoint) %>%
  tidyr::pivot_wider(names_from = Group, values_from = MeanAbund) %>%
  mutate(RIF = RIF/CTRL,
         RM = RM/CTRL,
         CTRL = 1) %>%
  select(Class, Timepoint, CTRL, RIF, RM) %>%
  tidyr::pivot_longer(cols = CTRL:RM, names_to = "Group") %>%
  mutate(Group = factor(Group, levels = c("CTRL", "RIF", "RM")))


betweenClassplot <- ggplot(ds2, aes(y= Class, x = Group, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),oob = squish, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Timepoint)

ggsave("reporting/11-2020/deseq2/Class_between_plot.png", betweenClassplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/Class_between_plot.pdf", betweenClassplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
betweenClassSig <- left_join(ds2, betweenClass, by = "Class") %>%
  mutate(Timepoint.y = factor(Timepoint.y, levels = c(1,2,3), labels = c("1_follicular", "2_ovulation", "3_luteal"))) %>%
  filter(Timepoint.x == Timepoint.y) %>%
  select(Class, Timepoint = Timepoint.x, Group, value, compRM_CTRL, compRIF_CTRL)
write.table(betweenClassSig, file = "reporting/11-2020/deseq2/class_final_between-sig.csv", sep= ";")
rm(ds, ds2, dssum)

### Order level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Order") %>%
  psmelt()  %>%
  filter(Order %in% betweenOrder$Order) %>%
  select(Order, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Order, Sample, Group, Timepoint, Abundance) %>%
  group_by(Order, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Order, Group, Timepoint, MeanAbund) %>%
  group_by(Order, Timepoint) %>%
  tidyr::pivot_wider(names_from = Group, values_from = MeanAbund) %>%
  mutate(RIF = RIF/CTRL,
         RM = RM/CTRL,
         CTRL = 1) %>%
  select(Order, Timepoint, CTRL, RIF, RM) %>%
  tidyr::pivot_longer(cols = CTRL:RM, names_to = "Group") %>%
  mutate(Group = factor(Group, levels = c("CTRL", "RIF", "RM")))


betweenOrderplot <- ggplot(ds2, aes(y= Order, x = Group, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),oob = squish, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Timepoint)

ggsave("reporting/11-2020/deseq2/Order_between_plot.png", betweenOrderplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/Order_between_plot.pdf", betweenOrderplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
betweenOrderSig <- left_join(ds2, betweenOrder, by = "Order") %>%
  mutate(Timepoint.y = factor(Timepoint.y, levels = c(1,2,3), labels = c("1_follicular", "2_ovulation", "3_luteal"))) %>%
  filter(Timepoint.x == Timepoint.y) %>%
  select(Order, Timepoint = Timepoint.x, Group, value, compRM_CTRL, compRIF_CTRL)
write.table(betweenOrderSig, file = "reporting/11-2020/deseq2/Order_final_between-sig.csv", sep= ";")
rm(ds, ds2, dssum)


### Family level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Family") %>%
  psmelt()  %>%
  filter(Family %in% betweenFamily$Family) %>%
  select(Family, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Family, Sample, Group, Timepoint, Abundance) %>%
  group_by(Family, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Family, Group, Timepoint, MeanAbund) %>%
  group_by(Family, Timepoint) %>%
  tidyr::pivot_wider(names_from = Group, values_from = MeanAbund) %>%
  mutate(RIF = RIF/CTRL,
         RM = RM/CTRL,
         CTRL = 1) %>%
  select(Family, Timepoint, CTRL, RIF, RM) %>%
  tidyr::pivot_longer(cols = CTRL:RM, names_to = "Group") %>%
  mutate(Group = factor(Group, levels = c("CTRL", "RIF", "RM")))


betweenFamilyplot <- ggplot(ds2, aes(y= Family, x = Group, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),oob = squish, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Timepoint)

ggsave("reporting/11-2020/deseq2/Family_between_plot.png", betweenFamilyplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/Family_between_plot.pdf", betweenFamilyplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
betweenFamilySig <- left_join(ds2, betweenFamily, by = "Family") %>%
  mutate(Timepoint.y = factor(Timepoint.y, levels = c(1,2,3), labels = c("1_follicular", "2_ovulation", "3_luteal"))) %>%
  filter(Timepoint.x == Timepoint.y) %>%
  select(Family, Timepoint = Timepoint.x, Group, value, compRM_CTRL, compRIF_CTRL)
write.table(betweenFamilySig, file = "reporting/11-2020/deseq2/Family_final_between-sig.csv", sep= ";")
rm(ds, ds2, dssum)


### Genus level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  tax_glom(taxrank = "Genus") %>%
  psmelt()  %>%
  filter(Genus %in% betweenGenus$Genus) %>%
  select(Genus, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(Genus, Sample, Group, Timepoint, Abundance) %>%
  group_by(Genus, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")

ds2 <- dssum %>%
  select(Genus, Group, Timepoint, MeanAbund) %>%
  group_by(Genus, Timepoint) %>%
  tidyr::pivot_wider(names_from = Group, values_from = MeanAbund) %>%
  mutate(RIF = RIF/CTRL,
         RM = RM/CTRL,
         CTRL = 1) %>%
  select(Genus, Timepoint, CTRL, RIF, RM) %>%
  tidyr::pivot_longer(cols = CTRL:RM, names_to = "Group") %>%
  mutate(Group = factor(Group, levels = c("CTRL", "RIF", "RM")))


betweenGenusPlot <- ggplot(ds2, aes(y= Genus, x = Group, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-3,3),oob = squish, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Timepoint)

ggsave("reporting/11-2020/deseq2/Genus_between_plot.png", betweenFamilyplot, width = 5, height= 3)
ggsave("reporting/11-2020/deseq2/Genus_between_plot.pdf", betweenFamilyplot, width = 5, height= 3)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
betweenGenusSig <- left_join(ds2, betweenGenus, by = "Genus") %>%
  mutate(Timepoint.y = factor(Timepoint.y, levels = c(1,2,3), labels = c("1_follicular", "2_ovulation", "3_luteal"))) %>%
  filter(Timepoint.x == Timepoint.y) %>%
  select(Genus, Timepoint = Timepoint.x, Group, value, compRM_CTRL, compRIF_CTRL)
write.table(betweenGenusSig, file = "reporting/11-2020/deseq2/Genus_final_between-sig.csv", sep= ";")
rm(ds, ds2, dssum)


### ASV level
ds <- ps.ns %>%
  transform_sample_counts(function(x) x/sum(x)*100) %>%
  psmelt() %>%
  mutate(ASV = paste(Genus, Species, OTU, sep = "|")) %>%
  filter(OTU %in% betweenASV$ASV) %>%
  select(ASV, Sample, Abundance, Group, Timepoint)

dssum <- ds %>%
  select(ASV, Sample, Group, Timepoint, Abundance) %>%
  group_by(ASV, Group, Timepoint)%>%
  summarize(MeanAbund = mean(Abundance, na.rm=T),
            sd = sd(Abundance, na.rm=T),
            n = n(),min=min(Abundance),
            max=max(Abundance), .groups = "keep")


ds2 <- dssum %>%
  select(ASV, Group, Timepoint, MeanAbund) %>%
  group_by(ASV, Timepoint) %>%
  tidyr::pivot_wider(names_from = Group, values_from = MeanAbund) %>%
  mutate(RIF = RIF/CTRL,
         RM = RM/CTRL,
         CTRL = 1) %>%
  select(ASV, Timepoint, CTRL, RIF, RM) %>%
  tidyr::pivot_longer(cols = CTRL:RM, names_to = "Group") %>%
  mutate(Group = factor(Group, levels = c("CTRL", "RIF", "RM")))

betweenASVplot <- ggplot(ds2, aes(y= ASV, x = Group, fill = value))+
  geom_tile(colour = "black") +  
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_gradient2(high = jBuPuPalette[1],
                       low = jBuPuPalette[paletteSize],
                       mid = jBuPuPalette[paletteSize/2], midpoint = 1,limits = c(-4,4),
                       oob = squish_any, 
                       name = "Rescaled Abundance", na.value = "black")+
  coord_equal(ratio = 1)+
  facet_wrap(~Timepoint)

ggsave("reporting/11-2020/deseq2/asv_between_plot.png", betweenASVplot, width = 10, height= 5)
ggsave("reporting/11-2020/deseq2/asv_between_plot.pdf", betweenASVplot, width = 10, height= 5)

# The following table merges normalized Abundance information with significance information (for subsequent plot annotation)
betweenASVsig <- ds2 %>%
  mutate(taxon = ASV,
         ASV = stringr::str_split_fixed(ASV, "\\|",n = 3)[,3]) %>%
  left_join(betweenASV, by = "ASV") %>%
  mutate(Timepoint.y = factor(Timepoint.y, levels = c(1,2,3), labels = c("1_follicular", "2_ovulation", "3_luteal"))) %>%
  filter(Timepoint.x == Timepoint.y) %>%
  select(taxon, Timepoint = Timepoint.x, Group, value, compRM_CTRL, compRIF_CTRL)
write.table(betweenASVsig, file = "reporting/11-2020/deseq2/ASV_final_between-sig.csv", sep= ";")
rm(ds, ds2, dssum)



####---- Within plot mege ----####
withinplot <- plot_grid(
  plot_grid(
    withinPhylumplot + theme(legend.position = "none"),
    withinClassplot+ theme(legend.position = "none"),
    withinOrderplot+ theme(legend.position = "none"),
    withinFamilyplot + theme(legend.position = "none"),
    withinGenusplot + theme(legend.position = "none"),
    ncol = 1,
    align = "hv"),
  plot_grid(
    get_legend(withinPhylumplot),
    get_legend(withinClassplot),
    get_legend(withinOrderplot),
    get_legend(withinFamilyplot),
    get_legend(withinGenusplot),
    ncol =1, 
    align = "hv")
  , rel_widths = c(12,3)
)
ggsave("reporting/11-2020/deseq2/withinplot-all.pdf", withinplot, width = 7, height = 14)


####---- Between plot mege ----####
betweenplot <- plot_grid(
  plot_grid(
    betweenPhylumplot + theme(legend.position = "none"),
    betweenClassplot+ theme(legend.position = "none"),
    betweenOrderplot+ theme(legend.position = "none"),
    betweenFamilyplot + theme(legend.position = "none"),
    betweenGenusPlot + theme(legend.position = "none"),
    ncol = 1,
    align = "hv"),
  plot_grid(
    get_legend(betweenPhylumplot),
    get_legend(betweenClassplot),
    get_legend(betweenOrderplot),
    get_legend(betweenFamilyplot),
    get_legend(betweenGenusPlot),
    ncol =1, 
    align = "hv")
  , rel_widths = c(12,3)
)
ggsave("reporting/11-2020/deseq2/betweenplot-all.pdf", betweenplot, width = 7, height = 14)
