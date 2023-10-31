#import R programs
install.packages("phyloseq")
install.packages("vegan")
install.packages("dplyr")
install.packages("microbiome") 
library(devtools)
BiocManager::install("ANCOMBC")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"=TRUE)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages("pairwiseAdonis"); packageVersion("pairwiseAdonis")
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

library(qiime2R)
library(tidyverse)
library(vegan)
library("ANCOMBC")
library("qiime2R")
library("phyloseq")
library("vegan")
library("tidyverse")
library("microbiome") 
library("pairwiseAdonis"); packageVersion("pairwiseAdonis")

setwd("C:/Users/Adam/Desktop/Spring_2023/Tropical_Marine_Microbiomes/ProjectForRealzies")
# Read in the 16S qza files and clean them up for phyloseq
ASVtable_16S <- read_qza("table_mergeTL.qza")
ASVtable_16S <- ASVtable_16S$data # Extract the count data from list
ASVtaxa_16S <- read_qza("taxonomyTL.qza")
taxtable_16S <- ASVtaxa_16S$data %>% as_tibble() %>% separate(Taxon, sep=";",
                                                              c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
map = read.table("C:/Users/Adam/Desktop/Spring_2023/Tropical_Marine_Microbiomes/ProjectForRealzies/metadata_SCTLD_Keys16SrRNAgene.txt", row.names = 1, header = TRUE)

#enter that shit into phyloseq
physeq_16S <-phyloseq(otu_table(ASVtable_16S, taxa_are_rows= T),
                       tax_table(as.data.frame(taxtable_16S) %>% column_to_rownames("Feature.ID") %>%
                                   as.matrix()), sample_data(map))

#lets play with the data
phy = physeq_16S
phy = subset_samples(phy, Experiment=="original")
phy = filter_taxa(phy, function(x) sum(x > 5 ) > (0.15*length(x)), TRUE)

#number of taxa
ntaxa(phy)

#number of samples
nsamples(phy)

#number of variables
sample_variables(phy)

############### Modifying tutorials to make own stuff from here

#steph diplo dicho
#for the different coral types 
steph = subset_samples(phy, CoralSpecies =="Stephanocoenia_intersepta")
diplo = subset_samples(phy, CoralSpecies == "Diploria_labyrinthiformis")
dicho = subset_samples(phy, CoralSpecies == "Dichocoenia_stokesii")

steph_clr = microbiome::transform(steph, 'clr')
diplo_clr = microbiome::transform(diplo, 'clr')
dicho_clr = microbiome::transform(dicho, 'clr')

#####plot relative abundance in coral species as function of Condition!
coral_phy = subset_samples(phy, Type =="Tissue")
#coral_phy = filter_taxa(coral_phy, function(x) sum(x > 10 ) > (0.15*length(x)), TRUE)
coral_phy = tax_glom(coral_phy, taxrank = "Phylum")
ps_rel_abund = phyloseq::transform_sample_counts(coral_phy, function(x){x / sum(x)})
phyloseq::plot_bar(ps_rel_abund, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack", color = "black") +
  labs(x = "", y = "Relative Abundance\n") +
  theme(legend.title = element_text(size = 20), 
        legend.text = element_text(size = 15), 
        plot.title = element_text(size = 40)) +
  facet_wrap(~ Condition, scales = "free") +
  theme(panel.background = element_blank())

### Alpha Diversity
# Subset the phyloseq object to exclude samples where Condition is NA
phy_sub <- subset_samples(phy, !is.na(meta(phy)$Condition))

#Subset again to see individual species, or skip it to view all at once. 
phy_sub <- subset_samples(phy_sub, CoralSpecies == "Stephanocoenia_intersepta" & !is.na(Condition))
phy_sub <- subset_samples(phy_sub, CoralSpecies == "Diploria_labyrinthiformis" & !is.na(Condition))
phy_sub <- subset_samples(phy_sub, CoralSpecies == "Dichocoenia_stokesii" & !is.na(Condition))

#plotting richness, be sure to change x and color by what you're trying to evaluate
plot_richness(phy_sub, x="Condition",color = "CoralSpecies", measure=c("Chao1")) +
  geom_boxplot(alpha=0.2)

#Use this to make the figure with just S. intercepta
plot_richness(phy_sub, x="Condition", measure=c("Chao1")) +
  geom_boxplot(alpha=0.2, colour = "Blue")

#if you want to show what is significant and then what you want to compare
symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
comparisons_material = list( c("AH", "DL"), c("AH", "DU"), c("DL", "DU"))

library(ggpubr)

#then rerun richness
plot_richness(phy_sub, x="Condition",color = "CoralSpecies", measure=c("Chao1")) +
  geom_boxplot(alpha=0.2)+
  stat_compare_means(method = "wilcox.test", comparisons = comparisons_material, label = "p.signif", symnum.args = symnum.args)

##### Venn Diagrams

unique(meta(phy)$CoralSpecies)

phynNA <- subset_samples(phy, !is.na(meta(phy)$CoralSpecies))

phynNA <- subset_samples(phynNA, Condition=="AH")


#find the core members of the microbiome
#convert to relative abundance
pseq_rel = microbiome::transform(phynNA, "compositional")

pseq_rel

#make a variable for all of the conditions you want to compare
stuff <- unique(as.character(meta(phy)$Zone[!is.na(meta(phy)$Zone)]))

#make a for loop to go through each bit of "stuff" one by one and combine identified core taxa into a list
list_core <- c() # an empty object to store information
for (n in stuff){ # for each variable n in stuff
  #print(paste0("Identifying Core Taxa for ", n))
  
  ps.sub <- subset_samples(pseq_rel, stuff== n) # Choose sample from DiseaseState by n
  
  core_m <- core_members(ps.sub, # ps.sub is phyloseq selected with only samples from g 
                         detection = 0.001, # 0.001 in atleast 90% samples 
                         prevalence = 0.10) #prevelence really matters--how much of the sample can it make upt 
  print(paste0("No. of core taxa in ", n, " : ", length(core_m))) # print core taxa identified in each DiseaseState.
  list_core[[n]] <- core_m # add to a list core taxa for each group.
  #print(list_core)
}

library(VennDiagram)

# Get the core taxa for each group
Endemic_core <- list_core[["Endemic"]]
Epidemic_core <- list_core[["Epidemic"]]
Vulnerable_core <- list_core[["Vulnerable"]]

# Create a Venn Diagram
venn.diagram(
  x = list(
    Endemic = Endemic_core,
    Epidemic = Epidemic_core,
    Vulnerable = Vulnerable_core
  ),
  filename = "core_taxa_venn.png"
)

######## PCOA
###################################### AWW HELL YIS ALL OF THEM WORK

#make a pcoa to evaluate each factor
phy_clr = microbiome::transform(phy, 'clr')
phy_ord = ordinate(phy_clr, "RDA", "euclidean")
plot_ordination(phy,phy_ord, type = "samples", color="Zone", shape="Type")

#Just Sediment by Zone
ord_sed = subset_samples(phy, Type=="Sediment")
ord_sed = microbiome::transform(ord_sed, 'clr')
ord_sed = ordinate(ord_sed, "RDA", "euclidean")
plot_ordination(phy,ord_sed, type = "samples", color="ReefName", shape="Zone")  +
  geom_point(size = 5) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.title = element_text(size = 20), 
        legend.text = element_text(size = 25), 
        plot.title = element_text(size = 40)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("Sediment Microbial Community Composition by Zone")

#Just Water by Zone
ord_water = subset_samples(phy, Type=="Water")
ord_water = microbiome::transform(ord_water, 'clr')
ord_water = ordinate(ord_water, "RDA", "euclidean")
plot_ordination(phy,ord_water, type = "samples", color="ReefName", shape="Zone") +
  geom_point(size = 5) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.title = element_text(size = 20), 
        legend.text = element_text(size = 25), 
        plot.title = element_text(size = 40)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("Aqueous Microbial Community Composition by Zone")

#Apparently Healthy Tissue by Zone for each Species
ord_tis = subset_samples(phy, Type=="Tissue")

### Skip these to run all, or pick one of the following three. I don't have the patience to build a function right now
ord_tis = subset_samples(ord_tis, CoralSpecies =="Stephanocoenia_intersepta")
ord_tis = subset_samples(ord_tis, CoralSpecies == "Diploria_labyrinthiformis")
ord_tis = subset_samples(ord_tis, CoralSpecies == "Dichocoenia_stokesii")

ord_tis = subset_samples(ord_tis, Condition=="AH")
ord_tis = microbiome::transform(ord_tis, 'clr')
ord_tis = ordinate(ord_tis, "RDA", "euclidean")

#Remember to change the species name in the title... because I don't have the patience to build a function right now
plot_ordination(phy, ord_tis, type = "samples", color = "ReefName", shape = "Zone") +
  geom_point(size = 15) +
  theme(axis.line = element_line(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.title = element_text(size = 20), 
        legend.text = element_text(size = 25), 
        plot.title = element_text(size = 40)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +  
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  ggtitle("Species Name Tissue Microbial Community Composition by Zone")

#Pairwise Adonis 
  #Can't believe the dumb mistake I was making with this... You don't wanna know. I'm content to wipe it from my memory and forget it ever happened.

#trying to get pairwise adonis to work

all_clr = microbiome::transform(phy, 'clr')

all_data= as(sample_data(all_clr), "data.frame")
all_data= as(sample_data(all_clr), "data.frame")
#adonis(t(otu_table(all_clr))~ Zone,  
#       data = all_data, permutations = 999, 
#       method = "euclidean") 

#AH Tissue
env = subset_samples(phy, Conditions!="AH")
env_clr = microbiome::transform(env, 'clr')
dist.uf <- phyloseq::distance(env_clr, method = "euclidean")
pairwise.adonis(t(otu_table(env_clr)), sample_data(env_clr)$Zone, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

#By species?
dist.uf <- phyloseq::distance(dicho_clr, method = "euclidean")
pairwise.adonis(t(otu_table(dicho_clr)), sample_data(dicho_clr)$Zone, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

dist.uf <- phyloseq::distance(steph_clr, method = "euclidean")
pairwise.adonis(t(otu_table(steph_clr)), sample_data(steph_clr)$Zone, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

dist.uf <- phyloseq::distance(diplo_clr, method = "euclidean")
pairwise.adonis(t(otu_table(diplo_clr)), sample_data(diplo_clr)$Zone, sim.method = "euclidean",
                p.adjust.m = "bonferroni")


#Significant differences as expected

#Let's compare all the 5 sample types even though my PCOAs only look at Apparently Healthy tissue. Maybe we'll see something interesting
#env = subset_samples(phy, Type!="Tissue")
#env_clr = microbiome::transform(env, 'clr')

#Not sure why this subset isn't working properly so I'm just gonna use all and truncate the table. It's not ideal but I gotta have some values to talk about!

dist.uf <- phyloseq::distance(all_clr, method = "euclidean")
pairwise.adonis(t(otu_table(all_clr)), sample_data(all_clr)$Conditions, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

#Hey look something interesting. No difference between EpAH and EpDU

#Let's look at sediment
env = subset_samples(phy, Type!="Sediment")
env_clr = microbiome::transform(env, 'clr')
dist.uf <- phyloseq::distance(env_clr, method = "euclidean")
pairwise.adonis(t(otu_table(env_clr)), sample_data(env_clr)$Zone, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

#And water
env = subset_samples(phy, Type!="Water")
env_clr = microbiome::transform(env, 'clr')
dist.uf <- phyloseq::distance(env_clr, method = "euclidean")
pairwise.adonis(t(otu_table(env_clr)), sample_data(env_clr)$Zone, sim.method = "euclidean",
                p.adjust.m = "bonferroni")


#########################################################################################
########################################################
################### Differential Abundance

#Fighting with R to do the damn thing

install.packages("ancombc")
install.packages("ancombc", repos = "http://cran.us.r-project.org")
install.packages("remotes")
remotes::install_version("ancombc", version = "1.1.3")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ancombc")

library(ANCOMBC)
library(dplyr)

if (!requireNamespace("phyloseq", quietly = TRUE))
  install.packages("phyloseq")
library(phyloseq)

# Subset for AH, sediment, and water
AH <- subset_samples(phy, Condition =="AH")
sed_phy = subset_samples(phy, Type =="Sediment")
wat_phy = subset_samples(phy, Type =="Water")



# Subset the phyloseq object to the Family level
sed_family <- tax_glom(sed_phy, taxrank = "Family")
wat_family <- tax_glom(wat_phy, taxrank = "Family")

sample_variables(sed_family)

# Run ancombc with prevalence cutoff of 90% for sediment
sedout <- ancombc(tax_level = "Family",
               assay_name = "counts",
               phyloseq = sed_phy,
               formula = "Zone",
               p_adj_method = "holm",
               group = "Zone",
               lib_cut = 1000,
               struc_zero = TRUE,
               neg_lb = TRUE,
               conserve = FALSE,
               alpha = 0.05,
               global = TRUE,
               prv_cut = 0.9)

#And for water, because I'm too lazy to write a function
watout <- ancombc(tax_level = "Family",
               assay_name = "counts",
               phyloseq = wat_phy,
               formula = "Zone",
               p_adj_method = "holm",
               group = "Zone",
               lib_cut = 1000,
               struc_zero = TRUE,
               neg_lb = TRUE,
               conserve = FALSE,
               alpha = 0.05,
               global = TRUE,
               prv_cut = 0.9)

#And for apparently healthy tissue because I'm too lazy to write a function
AHout <- ancombc(tax_level = "Family",
               assay_name = "counts",
               phyloseq = AH,
               formula = "Zone",
               p_adj_method = "holm",
               group = "Zone",
               lib_cut = 1000,
               struc_zero = TRUE,
               neg_lb = TRUE,
               conserve = FALSE,
               alpha = 0.05,
               global = TRUE,
               prv_cut = 0.9)

#Also for tissue types for final cross-comparison... because I'm too lazy to write a function even though it probably would have taken less time than writing all these comments about it
Tout <- ancombc(tax_level = "Family",
                assay_name = "counts",
                phyloseq = phy,
                formula = "Condition",
                p_adj_method = "holm",
                group = "Condition",
                lib_cut = 1000,
                struc_zero = TRUE,
                neg_lb = TRUE,
                conserve = FALSE,
                alpha = 0.05,
                global = TRUE,
                prv_cut = 0.9)

#Let's look at tissue types by zone? Oh wait never mind, I can only do this with the Epidemic samples because there's no lesioned tissue samples from corals without SCTLD! Derp
Epi = subset_samples(phy, Zone=="Epidemic")
End = subset_samples(phy, Zone=="Endemic")
Vul = subset_samples(phy, Zone=="Vulnerable")

Epiout <- ancombc(tax_level = "Family",
                assay_name = "counts",
                phyloseq = Epi,
                formula = "Condition",
                p_adj_method = "holm",
                group = "Condition",
                lib_cut = 1000,
                struc_zero = TRUE,
                neg_lb = TRUE,
                conserve = FALSE,
                alpha = 0.05,
                global = TRUE,
                prv_cut = 0.9)


library(tibble)

#Let's look at the tables! Well... they're not actually tables now... They're lists. Wish they were tables so I could flip them.
#out$res$diff_abn

sedout$res$diff_abn
watout$res$diff_abn
AHout$res$diff_abn
Tout$res$diff_abn

Epiout$res$diff_abn #Welp so much for that

#Clean them up. InTeRcEpT *insert mocking SpongeBob GIF

colnames(sedout$res$diff_abn)[2] <- "ZoneEndemic"
colnames(watout$res$diff_abn)[2] <- "ZoneEndemic"
colnames(AHout$res$diff_abn)[2] <- "ZoneEndemic"
colnames(Tout$res$diff_abn)[2] <- "ConditionAH"
colnames(Epiout$res$diff_abn)[2] <- "ConditionAH"

###Couldn't get this stuff to work the way I wanted so I'm cross referencing the tables above to get at the taxa maybe hopefully.

#Find out which taxa are associated with AH tissue in Endemic and Epidemic but not in Vulnerable
AHEE <- AHout$res$diff_abn[ AHout$res$diff_abn$ZoneEndemic == TRUE & AHout$res$diff_abn$ZoneEpidemic == TRUE & AHout$res$diff_abn$ZoneVulnerable == FALSE, ]

#Find out which taxa are associated with sediment in Endemic and Epidemic but not in Vulnerable
sedEE <- sedout$res$diff_abn[ sedout$res$diff_abn$ZoneEndemic == TRUE & sedout$res$diff_abn$ZoneEpidemic == TRUE & sedout$res$diff_abn$ZoneVulnerable == FALSE,]

#Find out which taxa are associated with water in Endemic and Epidemic but not in Vulnerable
watEE <- watout$res$diff_abn[ watout$res$diff_abn$ZoneEndemic == TRUE & watout$res$diff_abn$ZoneEpidemic == TRUE & watout$res$diff_abn$ZoneVulnerable == FALSE, ]

#Prevalence differences in AH, associated with disease and not vulnerable
AHEE
#Differential abundance in AH
AHout

#Prevalence differences in sediment, associated with disease and not vulnerable
sedEE
#Differential abundance in sediment
sedout

#Prevalence differences in water, associated with disease and not vulnerable
watEE
#Differential abundance in water
watout

Tout$res$diff_abn
Tout

Epiout$res$diff_abn
Epiout

######Disregard everything beyond this point except if you want to figure out what I was doing wrong with the bits that I didn't nix, and then laugh at whatever dumb mistake I was making.
# Convert tax_table to a data frame and add tax_id column
tax_table_df <- as.data.frame(tax_table(sed_family))
tax_table_df$OTU <- rownames(tax_table(sed_family))

head(tax_table_df)

# Join tax_table_df to ancom_taxa
ancom_taxa <- out$res$diff_abn %>%
  left_join(tax_table_df, by = "OTU") %>%
  .$Family

# Extract the differentially abundant taxa at the Family level
ancom_taxa <- out$res$diff_abn %>%
  as.data.frame() %>%
  tibble::rownames_to_column("taxid") %>%
  left_join(tax_table(sed_family), by = "taxid") %>%
  .$Family

head(tax_table(sed_family))

# Extract the differentially abundant taxa at the Family level
ancom_taxa <- out$res$diff_abn %>%
  as.data.frame() %>%
  tibble::rownames_to_column("taxid") %>%
  left_join(tax_table(sed_family), by = "taxid", copy = TRUE) %>%
  .$Family

ancom_taxa


out <- ancombc2(data = sed_family,
                tax_level = "Family",
                fix_formula = "~Zone",
                prv_cut = 0.9)


#####

# Subset the phyloseq object to the Family level
sed_family <- tax_glom(sed_phy, taxrank = "Family")

# Run ancombc2 with prevalence cutoff of 90%
out <- ancombc2(data = sed_family,
                tax_level = "Family",
                fix_formula = "~ Zone",
                prv_cut = 0.9)

ancom_taxa <- out$res$diff_abn %>%
  as.data.frame() %>%
  tibble::rownames_to_column("taxid") %>%
  left_join(tax_table(sed_family), by = "taxid") %>%
  filter(Zone) %>%
  .$Family

####