library(tidyverse) ; packageVersion("tidyverse") # 1.3.2
library(phyloseq) ; packageVersion("phyloseq") # 1.42.0
library(vegan) ; packageVersion("vegan") # 2.6.4
library(dendextend) ; packageVersion("dendextend") # 1.16.0
library(viridis) ; packageVersion("viridis") # 0.6.2
library("pairwiseAdonis"); packageVersion("pairwiseAdonis") # 0.4
library(ggpubr)
library(ANCOMBC)
library(microbiome)
library(eulerr)

setwd("~/Users/Trist/OneDrive/Desktop/Microbiomes")

#uploading your taxonomy and count data
asv_tab = read.delim("ASVs_counts.tsv", sep = "\t", header=TRUE) #do for all
asv_tax = read.delim("ASVs_taxonomy.tsv", sep="\t", header=TRUE)

#need to construct around your data using your metadata table from Run Selector
sample_info_tab = read.csv("New metadata table.csv", header=TRUE)


#need to make columns the rownames all of them
#need to figure out what the row names are for all your data, likely if you followed the tutorial the tax and tab are "X")
#for your sample data file the first row needs to be the SRR that is listed in your ASV files
asv_tab = asv_tab %>% tibble::column_to_rownames("X")
asv_tax = asv_tax %>% tibble::column_to_rownames("X")
sample_df = sample_info_tab %>% tibble::column_to_rownames("Labcode")


#transforming your taxonomy data into a matrix for entry into phyloseq
asv_mat = as.matrix(asv_tab)
tax_mat = as.matrix(asv_tax)

#transform into phyloseq objects
asv = otu_table(asv_mat, taxa_are_rows=TRUE)
tax=tax_table(tax_mat)
samples=sample_data(sample_df)


#make a phyloseq object
phy =phyloseq(asv, tax, samples)

#remove chloroplast and mitochondria from your phyloseq object
phy <- phy %>% subset_taxa( family!= "mitochondria" & class!="Chloroplast") 


#number of taxa
ntaxa(phy)

#number of samples
nsamples(phy)

#number of variables
sample_variables(phy)


#plotting richness, be sure to change x and color by what you're trying to evaluate
chao = plot_richness(phy, x="Species",color = "Life_Stage", measure=c("Chao1")) +
  geom_boxplot(alpha=0.2)+
  labs(title = "Alpha Diversity Among Samples", color = "Life Stage")

chao
dev.off()

pdf("/Users/Trist/OneDrive/Desktop/Microbiomes/richnesse.pdf", height = 7, width=13)
plot(chao)
dev.off()

#if you want to show what is significant and then what you want to compare
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
comparisons_material = list( c("Dosidicus gigas", "Ommastrephes sp."), c("Dosidicus gigas", "SD Complex"), c("Dosidicus gigas", "Sthenotuethis oualanensis"), c("Dosidicus gigas", "Sthenotuethis pteropus"), c("Dosidicus gigas", "Todarodes sagittatus"), c("Sthenoteuthis oualanensis", "Ommastrephes sp."), c("Sthenotuethis oualanensis", "SD complex"), c("Sthenotuethis oualanensis", "Sthenotuethis pteropus"), c("Sthenotuethis oualanensis", "Todarodes sagittatus"), c("Ommastrephes sp.", "SD complex"), c("Ommastrephes sp.", "Sthenotuethis pteropus"), c("Ommastrephes sp.", "Todarodes sagittatus"), c("SD complex", "Sthenotuethis pteropus"), c("SD complex", "Todarodes sagittatus"), c("Sthenotuethis pteropus, Todarodes sagittatus"))
#change list( c(")) to comparisons of all variables in "source" aka species

comparisons_material2 = list( c("Early Paralarvae", "Late Paralarvae"), c("Early Paralarvae", "Subadults and adult"), c("Late Paralarvae", "Subadults and adult"))
#then rerun richness
plot_richness(phy, x="Species",color = "Life_Stage", measure=c("Chao1")) +
  geom_boxplot(alpha=0.2)+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material2, label = "p.signif", symnum.args = symnum.args)

#make a stacked barplot of the data at the phylum level
phy_phylum = tax_glom(phy, taxrank = "phylum",NArm=FALSE)
rel_abund_phylum = phyloseq::transform_sample_counts(phy_phylum,
                                                           function(x){x / sum(x)})

phy_order = tax_glom(phy, taxrank = "order",NArm=FALSE)
rel_abund_order = phyloseq::transform_sample_counts(phy_order,
                                                     function(x){x / sum(x)})

phy_class = tax_glom(phy, taxrank = "class",NArm=FALSE)
rel_abund_class = phyloseq::transform_sample_counts(phy_class,
                                                    function(x){x / sum(x)})
phy_family = tax_glom(phy, taxrank = "family",NArm=FALSE)
rel_abund_family = phyloseq::transform_sample_counts(phy_family,
                                                    function(x){x / sum(x)})

relabundance<- phyloseq::plot_bar(rel_abund_phylum, fill = "phylum")+
  geom_bar(aes(color = phylum, fill = phylum), stat = "identity", position = "stack", color = "black") +
  labs(x = "", y = "Relative Abundance\n", fill = "Phylum") +
  facet_wrap(~ Life_Stage, scales = "free") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))+
  ggtitle("Relative Abundance at Phylum Level")
relabundance

relabundanceorder<- phyloseq::plot_bar(rel_abund_order, fill = "order")+
  geom_bar(aes(color = order, fill = order), stat = "identity", position = "stack", color = "black") +
  labs(x = "", y = "Relative Abundance\n", fill = "Order") +
  facet_wrap(~ Life_Stage, scales = "free") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))+
  ggtitle("Relative Abundance at Order Level")
relabundanceorder


relabundanceclass<- phyloseq::plot_bar(rel_abund_class, fill = "class")+
  geom_bar(aes(color = class, fill = class), stat = "identity", position = "stack", color = "black") +
  labs(x = "", y = "Relative Abundance\n", fill = "Class") +
  facet_wrap(~ Life_Stage, scales = "free") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))+
  ggtitle("Relative Abundance at Class Level")
relabundanceclass

relabundancefamily<- phyloseq::plot_bar(rel_abund_family, fill = "family")+
  geom_bar(aes(color = family, fill = family), stat = "identity", position = "stack", color = "black") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Life_Stage, scales = "free") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))+
  ggtitle("Relative Abundance at Family Level")
relabundancefamily

relabundanceorder

pdf("/Users/Trist/OneDrive/Desktop/Microbiomes/relabundance.pdf", height = 7, width=13)
plot(relabundance)
dev.off()

relabundance<- phyloseq::plot_bar(rel_abund_phylum, fill = "phylum")+
  geom_bar(aes(color = phylum, fill = phylum), stat = "identity", position = "stack", color = "black") +
  labs(x = "", y = "Relative Abundance\n") +
  facet_wrap(~ Species + Life_Stage, scales = "free") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 12))
relabundance


pdf("/Users/Trist/OneDrive/Desktop/Microbiomes/relabundance.pdf", height = 7, width=13)
plot(relabundance)
dev.off()




#find the core members of the microbiome
#convert to relative abundance
pseq_rel = microbiome::transform(phy, "compositional")

#make a variable for all of the conditions you want to compare
stuff = unique(as.character(meta(phy)$Life_Stage))

stuff2 = unique(as.character(meta(phy)$Species!="Ommastrephes sp."))

sub_phy<-subset_samples(phy, phy@sam_data$Species!="Ommastrephes sp.")
stuff3 = unique(as.character(meta(sub_phy)$Species))
#make a for loop to go through each bit of "stuff" one by one and combine identified core taxa into a list
 # an empty object to store information
# for each variable n in stuff #print(paste0("Identifying Core Taxa for ", n))

list_core <- c() 
for (n in stuff){ 
  
  
  ps.sub <- subset_samples(pseq_rel, stuff== n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.09) #prevelence really matters--how much of the sample can it make upt 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

venn <- plot(venn(list_core))
venn

ep = data.table(list_core$`Early Paralarvae`)
ep

lp = data.table(list_core$`Late Paralarvae`)
lp

as = data.table(list_core$`Subadults and adult`)
as
#may not have core taxa for all samples. Can try to look at for individual species, life stage

#can also make a list of your factors and the color to give to each
#list of all your factors with the color to asign to each (replace nonCRC, CRC, H and add whatever other factors needed)
mycols = c(nonCRC="darkseagreen1", CRC="darkslategray3", H="goldenrod2") 
plot(venn(list_core),
     fills = mycols)

list_core <- c() 
for (n in stuff3){ 
  
  
  ps.sub <- subset_samples(pseq_rel, stuff3== n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.10) #prevelence really matters--how much of the sample can it make upt 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

venn <- plot(venn(list_core))
venn

#may not have core taxa for all samples. Can try to look at for individual species, life stage

#can also make a list of your factors and the color to give to each
#list of all your factors with the color to asign to each (replace nonCRC, CRC, H and add whatever other factors needed)
mycols = c(nonCRC="darkseagreen1", CRC="darkslategray3", H="goldenrod2", A="slategray4") 
plot(venn(list_core),
     fills = mycols)

#make a pcoa to evaluate each factor
phy_clr = microbiome::transform(phy, 'clr')
phy_ord = ordinate(phy_clr, "RDA", "euclidean")
plot_ordination(phy,phy_ord, type = "Labcode", color="Species", shape="Life_Stage")+
  labs(title = "PCoA Based on Species and Life Stage", shape="Life Stage", color = "Species")


pdf("/Users/Trist/OneDrive/Desktop/Microbiomes/pcoa.pdf", height = 7, width=13)
plot(xx)
dev.off()

#Test whether the seasons differ significantly from each other using the permutational ANOVA 
#Pairwise Adonis with Species
dist.uf <- phyloseq::distance(phy_clr, method = "euclidean")
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$Species, sim.method = "euclidean",
                p.adjust.m = "bonferroni")
#Pairwise Adonis with Life_Stage
pairwise.adonis(t(otu_table(phy_clr)), sample_data(phy_clr)$Life_Stage, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

#subset early paralarvae
phy_paralarvae<-subset_samples(phy, phy@sam_data$Life_Stage=='Early Paralarvae')


#find the core members of the microbiome
#convert to relative abundance
pseq_rel = microbiome::transform(phy_paralarvae, "compositional")

#make a variable for all of the conditions you want to compare
stuff = unique(as.character(meta(phy_paralarvae)$Life_Stage))


#make a for loop to go through each bit of "stuff" one by one and combine identified core taxa into a list
# an empty object to store information
# for each variable n in stuff #print(paste0("Identifying Core Taxa for ", n))

list_core <- c() 
for (n in stuff){ 
  
  
  ps.sub <- subset_samples(pseq_rel, stuff== n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.10) #prevelence really matters--how much of the sample can it make upt 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

plot(venn(list_core))

# filtered phyloseq
phy_filter = filter_taxa(phy, function(x) sum(x > 3) > (0.015*length(x)), TRUE)

phy_clr_filter = microbiome::transform(phy_filter, 'clr')
phy_ord_filter = ordinate(phy_clr_filter, "RDA", "euclidean")
plot_ordination(phy_filter,phy_ord_filter, type = "Labcode", color="Species", shape="Life_Stage")

pairwise.adonis(t(otu_table(phy_clr_filter)), sample_data(phy_clr_filter)$Life_Stage, sim.method = "euclidean",
                p.adjust.m = "bonferroni")
pairwise.adonis(t(otu_table(phy_clr_filter)), sample_data(phy_clr_filter)$Species, sim.method = "euclidean",
                p.adjust.m = "bonferroni")
#Bubble plot
library(randomcoloR)
n <- 34
palette <- distinctColorPalette(n)


type_taxa <- phy %>% 
  tax_glom(taxrank = "order") %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>% 
  arrange(order) 
summary(type_taxa$Sample)
head(type_taxa$Sample)
type_taxa$Sample

xx <- ggplot(type_taxa, aes(x = factor(Life_Stage, levels=c("Early Paralarvae", "Late Paralarvae", "Subadults and adult")), y = factor(order))) + geom_point(aes(size = Abundance, fill=order), alpha = 0.9, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,2,3,5,20)) + 
  labs( x= "Life Stage", y = "Order", size = "Relative Abundance %", fill="Phylum") + theme_bw() +  
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Prokaryotic Orders within Ommastrephid Gut") 
xx

pdf("/Users/Trist/OneDrive/Desktop/Microbiomes/orderbubble.pdf", height = 7, width=13)
plot(xx)
dev.off()


type_taxa <- phy %>% 
  tax_glom(taxrank = "phylum") %>% 
  transform_sample_counts(function(x) {x/sum(x)} ) %>% 
  psmelt() %>% 
  arrange(phylum) 
summary(type_taxa$Sample)
head(type_taxa$Sample)
type_taxa$Sample

xx <- ggplot(type_taxa, aes(x = factor(Life_Stage, levels=c("Early Paralarvae", "Late Paralarvae", "Subadults and adult")), y = factor(phylum))) + geom_point(aes(size = Abundance, fill=phylum), alpha = 0.9, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "Life Stage", y = "Order", size = "Relative Abundance %", fill="Phylum") + theme_bw() +  
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Prokaryotic Phyla within Ommastrephid Gut") 
xx

xx <- ggplot(type_taxa, aes(x = factor(Species, levels=c("Ommastrephes sp.", "Dosidicus gigas", "Todarodes sagittatus", "SD complex", "Sthenoteuthis pteropus")), y = factor(phylum))) + geom_point(aes(size = Abundance, fill=phylum), alpha = 0.9, shape = 21) + 
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,10), breaks = c(1,10,50,75)) + 
  labs( x= "Life Stage", y = "Phylum", size = "Relative Abundance %", fill="Phylum") + theme_bw() +  
  scale_fill_manual(values = palette) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ggtitle("Prokaryotic Phyla within Ommastrephid Gut") 
xx


#Ampvis2
if(!require("devtools"))
  install.packages("devtools")
#source the phyloseq_to_ampvis2() function from the gist
devtools::source_gist("8d0ca4206a66be7ff6d76fc4ab8e66c6")

ampvis2_obj <- phyloseq_to_ampvis2(phy)
ampvis2_obj


xx<- amp_heatmap(ampvis2_obj, group_by = "Life_Stage", tax_add = c("Phylum", "Order", "Genus"), tax_aggregate = "OTU",
                 tax_show = 15, color_vector = c("snow", "steelblue")) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
xx

xxSpecies<- amp_heatmap(ampvis2_obj, group_by = "Species", tax_add = c("Phylum", "Order", "Genus"), tax_aggregate = "OTU",
                 tax_show = 15, color_vector = c("snow", "steelblue")) + theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
xxSpecies

#network diagram
ig = make_network(phy, "samples", max.dist=1.0005, dist.fun = "euclidean")
plot_network(ig,phy, color="Species", shape="Life_Stage", line_weight=0.3,label=NULL)

# Data Analysis
psf_toadfish_extremes<-subset_samples(psf_toadfish, psf_toadfish@sam_data$Treatment!="35")



#Core microbiome
astr <- subset_samples(phy, Life_Stage=="Early Paralarvae") 
astr <-subset_samples(astr, Species=="Todarodes Sagittatus")

astr <- prune_taxa(taxa_sums(astr) > 0, astr)
astr <- microbiome::transform(astr, "compositional")
astr.core <- core(astr, detection = 0.1, prevalence = 0.90)
core.taxa <- taxa(astr.core)
class(core.taxa)
tax.mat <- tax_table(astr.core)
tax.df <- as.data.frame(tax.mat)
tax.df$OTU <- rownames(tax.df)
core.taxa.class <- dplyr::filter(tax.df, rownames(tax.df) %in% core.taxa)
knitr::kable(head(core.taxa.class))
# Core with compositionals:
prevalences <- seq(.07, 1, .07)
detections <- 100^seq(log10(1e-4), log10(.5), length = 100)
p1 <- plot_core(astr, 
                plot.type = "heatmap", 
                colours = c("white", "dodgerblue"),
                prevalences = prevalences, 
                detections = detections, min.prevalence = .5) +
  xlab("Detection Threshold (Relative Abundance (%))")
p1 <- p1 + theme_bw() + ylab("ASVs") + labs(title="Core Prokaryotic ASVs of Precipitates at 35 ppt") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
df <- p1$data 
list <- df$Taxa 
tax <- as.data.frame(tax_table(astr))
tax$ASV <- rownames(tax)
tax2 <- dplyr::filter(tax, rownames(tax) %in% list) 
tax.unit <- tidyr::unite(tax2, Taxa_level,c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV"), sep = "_;", remove = TRUE) #EDIT these levels to match your taxa file
tax.unit$Taxa_level <- gsub(pattern="[a-z]__",replacement="", tax.unit$Taxa_level)
df$Taxa <- tax.unit$Taxa_level
knitr::kable(head(df))
p1$data <- df
p2 <- plot(p1 + theme(axis.text.y = element_text(face="italic")))
p2





#save output of ANCOMBC
phy_sub<-subset_samples(phy,phy@sam_data$Species!="Ommastrephes sp.")
phy_sub<-subset_samples(phy_sub,phy_sub@sam_data$Species!="Todarodes sagittatus")
phy_sub<-subset_samples(phy_sub,phy_sub@sam_data$Species!="Sthenoteuthis pteropus")
phy_sub<-subset_samples(phy_sub,phy_sub@sam_data$Species!="Sthenoteuthis oualanensis")

phy_early<-subset_samples(phy, phy@sam_data$Life_Stage=="Early Paralarvae")

#adjust subset of locations for locations other than posterior
### Adjust below code "SD" and what not accordingly.
### VAR1 is the variable you are interested in. VAR2 (or VAR3, VAR4, etcs) are fixed variables you wish to account for.
output_combined_squid = ancombc2(data = phy, assay_name = "counts", tax_level = "order",
                           fix_formula = "Life_Stage", rand_formula = NULL,
                           p_adj_method = "fdr", pseudo = 0, pseudo_sens = TRUE,
                           prv_cut = 0.01, lib_cut = 0, s0_perc = 0.05,
                           group = "Life_Stage", struc_zero = TRUE, neg_lb = FALSE,
                           alpha = 0.05, n_cl = 2, verbose = TRUE,
                           global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = FALSE,
                           iter_control = list(tol = 1e-2, max_iter = 20, 
                                               verbose = TRUE),
                           em_control = list(tol = 1e-5, max_iter = 100),
                           lme_control = lme4::lmerControl(),
                           mdfdr_control = list(fwer_ctrl_method = "fdr", B = 100),
                           trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                       nrow = 2, 
                                                                       byrow = TRUE),
                                                                matrix(c(-1, 0, 1, -1),
                                                                       nrow = 2, 
                                                                       byrow = TRUE)),
                                                node = list(2, 2),
                                                solver = "ECOS",
                                                B = 100))

res_prim = output_combined_squid$res
res_prim
res_prim$taxon=as.factor(res_prim$taxon)


df_res = res_prim %>% dplyr::filter("diff_Life_StagesSubadults and adult") %>%
  dplyr::select(taxon, contains("Life_Stage")) 


df_fig_res = df_res %>% 
  arrange(desc(lfc_SpeciesSDComplex
  )) %>%
  mutate(direct = ifelse(lfc_SpeciesSDComplex
                         > 0, "Positive LFC", "Negative LFC"))


df_fig_res$taxon = factor(df_fig_res$taxon, levels = df_fig_res$taxon)
df_fig_res$direct = factor(df_fig_res$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

fig_res = df_fig_res %>%
  ggplot(aes(x = taxon, y = lfc_SpeciesSDComplex, fill = lfc_SpeciesSDComplex)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_SpeciesSDComplex - se_SpeciesSDComplex, ymax = lfc_SpeciesSDComplex + se_SpeciesSDComplex), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = "Prokaryotic Order", y = "Log fold change", 
       title = "LFC in SD Complex vs Other Species", subtitle = "Of ANCOM-BC signficant taxa corrected by Treatment") + scale_fill_viridis_c(option = "plasma") + scale_color_viridis_c(option = "plasma") +
  theme_bw() + 
  theme(legend.position = "none",
        panel.grid.minor.y = element_blank()) + coord_flip()

fig_res

pdf("/Users/Trist/OneDrive/Documents/Toadfish Microbiome/amcom3.pdf", height = 7, width=13)
plot(fig_res)
dev.off()