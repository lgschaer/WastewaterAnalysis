#packages used
library(phyloseq)
library(tidyverse)
library(vegan)
library(RColorBrewer)

#set theme
theme_set(theme_bw())

#load data into phyloseq

##sample data
sdata <- as.csv("/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Unk/sample_info.csv", row.names = 1, header = TRUE, sep = ",", check.names = TRUE, stringsAsFactors = TRUE)
head(sdata)

sdata2 <- sdata %>% 
  rownames_to_column(var = "SampleID") %>%
  as_tibble()
head(sdata2)

sdata3 <- sdata2 %>%
  mutate(
    SampleID,
    Sample_ID = SampleID
  )%>%
  column_to_rownames(var = "Sample_ID")
head(sdata3)

##sequence table
sequence_table <- readRDS("/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Unk/seqtab.rds")
colnames(sequence_table) <- NULL
sequence_table <- as.data.frame(sequence_table)
sequence_table[1:4,1:4]
sequence_table <- rownames_to_column(sequence_table, var = "SampleID")
sequence_table[1:4,1:4]
class(sequence_table)

sequence_table <- sequence_table %>%
  as.data.frame()
sequence_table[1:13,1:4]

#join metadata to sequence table so we can filter samples
seq_table <- sdata2 %>% 
  full_join(sequence_table, by = "SampleID")
seq_table[1:10,1:10]


#get rid of "SampleID" header, this is where we would subset if needed
subset_all <- seq_table %>%
  column_to_rownames(var = "SampleID")
subset_all[1:10,1:10]

subset_all$volume <- NULL
subset_all$sampletype <- NULL
subset_all$letter <- NULL
subset_all$depth.ft. <- NULL
subset_all$distance_ww <- NULL

subset_all[1:10,1:10]
class(subset_all)

subset_all <- as.matrix(subset_all)
m <- (colSums(subset_all, na.rm=TRUE) != 0)                   #T if colSum is not 0, F otherwise
subset_all[1:3,1:3]
nonzero <- subset_all[, m]                                    #all the non-zero columns
zero <- subset_all[, !m]                                      #all the zero columns
s_table <- nonzero[,colSums((nonzero)>0)>1, drop=FALSE]
nonzero[1:5,1:5]

##taxa table
taxa <- readRDS("/home/lgschaer/old/boatbiome/April-Boat-Microbiome/Unk/taxa.rds")
taxa_t<-t(taxa)
colnames(taxa_t) <- NULL
taxa_table<-t(taxa_t)
taxa_table <- as.matrix(taxa_table)
taxa_table[1:5,1:5]

#make phyloseq object
samdata = sample_data(sdata3)

colnames(nonzero) <- NULL
seqtab = otu_table(nonzero, taxa_are_rows = FALSE)

taxtab = tax_table(taxa_table)
rownames(taxtab) <- NULL

phyloseq_object_all = phyloseq(otu_table(seqtab), tax_table(taxtab), sample_data(samdata))
phyloseq_object_all

#normalize data
#Delete samples with a mean of less than 1000
samplesover1000_all <- subset_samples(phyloseq_object_all, sample_sums(phyloseq_object_all) > 1000)

#Check if there are OTUs with no counts, if so how many?
any(taxa_sums(samplesover1000_all) == 0)
sum(taxa_sums(samplesover1000_all) == 0)

#Prune OTUs with no counts 
prune_samplesover1000_all <- prune_taxa(taxa_sums(samplesover1000_all) > 0, samplesover1000_all)
any(taxa_sums(prune_samplesover1000_all) == 0)

#make sure seed is set the same each time, set to 81 here
rarefy_samplesover1000_all <- rarefy_even_depth(prune_samplesover1000_all, rngseed= 81, sample.size = min(sample_sums(prune_samplesover1000_all)))

ww <- rarefy_samplesover1000_all
ww

head(sample_data(ww))

ww2 <- subset_samples(ww, letter == "W" | letter == "C")
head(sample_data(ww2))
ww2

saveRDS(ww2, "/home/lgschaer/old/wastewater_analysis/ww2.RDS")

#count samples
sample_counts <- sample_data(ww) %>%
  group_by(sampletype) %>%
  count()
sample_counts

#assign colors to samples
sample_colors <- c("w-1-75mL" = "#FEEDDE", 
                   "w-2-300mL" = "#FDBE85", 
                   "w-3-200mL" = "#FD8D3C", 
                   "w-4-250mL" = "#D94701", 
                   "A1" = "#DEEBF7", 
                   "A2" = "#9ECAE1", 
                   "A3" = "#3182BD", 
                   "B1" = "#EFEDF5", 
                   "B2" = "#BCBDDC", 
                   "B3" = "#756BB1",
                   "C1" = "#EDF8E9",
                   "C2" = "#74C476", 
                   "C3" = "#006D2C")

letter_colors <- colorRampPalette(brewer.pal(5, "Greens"))(3)
letter_colors

sample_types <- c("lakewater", "wastewater")

sample_labels <- c("wastewater" = "Waste\nWater", "lakewater" = "Lake\nWater")


#dot plot of alpha diversity Observed OTUs and Shannon Diversity

ww2 %>%                                                                                #phyloseq object
  plot_richness(
    x = "letter",                                                                      #compare diversity of letter
    measures = c("Observed", "Shannon")) +                                             #choose diversity measures
  geom_point(aes(fill = SampleID),size = 6, shape = 21, show.legend = TRUE)+           #make violin plot, set fill aes to sampletype
 # geom_boxplot(width=0.1) +                                                           #add boxplot, set width
  theme_classic()+                                                                     #change theme to classic
  xlab(NULL)+                                                                          #no label on x-axis
  scale_fill_manual(values = sample_colors)+                                           #set fill colors
  ggtitle("Alpha Diversity") +                                                         #add title
  theme(strip.text = element_text(face = "bold", size = 25))+                          #adjust headings
  theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5),    #adjust x-axis text
        axis.title.y = element_text(size = 20),                                        #adjust y-axis title
        axis.text.y = element_text(size = 15),                                         #adjust y-axis text
        legend.text = element_text(size = 15),                                         #adjust legend text size
        legend.title = element_text(size = 22),                                        #adjust legend title size
        plot.title=element_text(size = 25, face = "bold", hjust = 0.5))                #change title size, face and position

#finding means and medians
richness <- ww2 %>%
  estimate_richness(measures = c("Observed", "Shannon")) %>%           #specify which measures
  rownames_to_column(var = "SampleID") %>%                             #add column name to SampleID column
  as_tibble() %>%
  mutate(letter = ifelse(SampleID == "C1"|SampleID == "C2"|SampleID == "C3", "C", "W"))      #add letter column to differentiate sample type
head(richness)

rich_summary <- richness %>%
  group_by(letter) %>%                                                 #group by sample type
  summarize(                                                           #add columns for stat summaries
    maxObs = max(Observed),
    maxShan = max(Shannon),
    minObs = min(Observed),
    minShan = min(Shannon),
    rangeObs = maxObs - minObs,
    meanObserved = mean(Observed),
    medianObserved = median(Observed),
    rangeShan = maxShan - minShan,
    meanShannon = mean(Shannon),
    medianShannon = median(Shannon)
  )
View(rich_summary)                                  #check that we have the columns we want

#ANOVA
set.seed(81)

#calculate alpha diversity measures
alpha_data_all <- ww2 %>%
  estimate_richness(measures = c("Observed", "Shannon")) %>%           #specify which measures
  rownames_to_column(var = "SampleID")                                #add column name to SampleID column
head(alpha_data_all)

data_all <- alpha_data_all %>% 
  mutate(
  SampleID = gsub(".", "-", SampleID, fixed = TRUE))%>%
  left_join(sdata3, by = "SampleID") 
head(data_all)

#Observed ANOVA
anova_obs <- aov(Observed ~ sampletype, data_all)
summary(anova_obs)

#Shannon ANOVA
anova_shan <- aov(Shannon ~ sampletype, data_all)
summary(anova_shan)

#tukeyHSD, post-hoc
TukeyHSD(anova_obs)
TukeyHSD(anova_shan)

#PERMANOVA
set.seed(81)

#subset phyloseq object, all samples by datatype
wl <- subset_samples(ww2, sampletype %in% c("wastewater", "lakewater"))



# Calculate bray curtis distance matrix, all samples
wl_bray <- phyloseq::distance(wl, method = "bray")


# make a data frame from the sample_data, all samples
wl.sam <- data.frame(sample_data(wl))


# Adonis test, all samples
adonis(wl_bray ~ sampletype, data = wl.sam)


#PCoA Plot

#ordination
all_pcoa <- ordinate(
  physeq = ww2, 
  method = "PCoA", 
  distance = "bray"
)

head(sample_data(ww2))

#plot
plot_ordination(
  physeq = ww2,                                                          #phyloseq object
  ordination = all_pcoa)+                                                #ordination
  geom_point(aes(fill = SampleID), shape = 21, size = 3) +                         #sets fill color to sampletype
  scale_fill_manual(values = sample_colors) +
  theme_classic() +                                                      #changes theme, removes grey background
  theme(                             
    legend.text = element_text(size = 20),                               #changes legend size
    legend.title = element_blank(),                                      #removes legend title
    legend.background = element_rect(fill = "white", color = "black"))+  #adds black boarder around legend
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20))+
  guides(fill = guide_legend(override.aes = list(shape = 21)))           #fills legend points based on the fill command



#EXPLORING TAXA

#filter out eukaryotes and mitochondria
justbacteria <- ww2 %>%
  subset_taxa(
    Kingdom == "Bacteria" &                   #only bacteria
      Family  != "mitochondria"               #filter out mitochondria
  )
justbacteria

#-------------------------------------------------#
#DETERMINE IF THERE ARE COMMON WASTEWATER ORGANISMS

#Aglomerate at genus level, melt to long format
genusabundance <- justbacteria %>%
  tax_glom(taxrank = "Genus") %>%                      # agglomerate at order level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Genus) %>%
  as.tibble()
head(genusabundance)


#Select necessary variables
sulf <- genusabundance %>%
  select(Phylum, Class, Order, Family, Genus, SampleID, Abundance) %>%
  group_by(Phylum, Class, Order, Family, Genus, SampleID) 
head(sulf)
#dim(sulf)


#Filter out core wastewater organisms
sulf2 <- sulf %>%
  ungroup() %>%
  mutate(Genus = as.character(Genus)) %>%
  mutate(Temp = str_replace(Genus, "Desulfo|desulfo|Thio|thio|Proteus|proteus|
                            Desulf|desulf|Sulf|Sulf|Byssovorax|Pedomicrobium|Rhizomicrobium|
                            Chloro|chloro|Epsilon|epsilon|Sediminibacterium|Flavobacterium|
                            Roseomonas|Trichococcus|Alkanidiges|Sphingomonas|Anaerolinea|
                            Ignavibacterium|Verrucomicrobium|Blastopirellula|Aquabacterium|Gp7|Gp17|incertae|
                            Syntrophorhabus|Luteimonas|Thauera|Lewinella|Longilinea|Levilinea|Gp4|Gp6|
                            Ohtaekwangia|Caldilinea|Zavarzinella|Pasteuria|Luteibacter|Hyphomicrobium|
                            Terrimonas|Haliscomenobacter|Zoogloea|Leptolinea|Prostheocobacter|Ferruginibacter|
                            Lactococcus|Azoarcus|Ferribacterium|Clostridium|Gemmatimonas|Nitrospra|Arcobacter|
                            Paludibacter|Leeia|Macellibacteroides|
                            Chromatium|Beggiatoa|Paracoccus|Achromatium|Campylobacter|Salmonella", "__xSulfurx__")) %>%
  separate(Temp, into = c("Garbage1", "Temp"), sep = "__x") %>%
  separate(Temp, into = c("Sulfur", "Garbage2"), sep = "x__") %>%
  filter(Sulfur == "Sulfur") %>%
  select(-c("Garbage1", "Garbage2", "Sulfur")) %>%
  filter(Abundance > 0.01)
head(sulf2)


#make vector for breaks and labels
b <- c("w-1-75mL", "w-2-300mL", "w-3-200mL", "w-4-250mL", "C1", "C2", "C3")

li <- c("w-1-75mL", "w-2-300mL", "w-3-200mL", "w-4-250mL", "C1", "C2", "C3")

la <- c("w-1-75mL" ="W1", "w-2-300mL" = "W2", "w-3-200mL" = "W3", "w-4-250mL" = "W4", "C1" = "C1", "C2" = "C2", "C3" = "C3")

#make color palette
colors1 <- c(
  "purple4", "darkcyan",     "orchid1",   "green",       "coral1",        "blue", "yellow", "red"
 # "grey47",  "cyan",         "yellow",    "darkgreen",   "palegoldenrod",     "tan4",
 # "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1",
 # "white",   "dodgerblue",   "firebrick", "yellowgreen", "magenta", "black", "lightpink"
)  

#Plot abundance of sulfur degrading organisms

#Bars reflect proportion of overall abundance represented
ggplot(sulf2)+
  geom_col(mapping = aes(x = SampleID, y = Abundance, fill = Genus), color = "black", position = "stack", show.legend = TRUE)+
  ylab("Proportion of Community") +
  ylim(c(0,1))+
  scale_fill_manual(values = colors1) +
  scale_x_discrete(
    breaks = b,
    labels = la,
    limits = li)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))

#Bars reflect proportion of this subset
ggplot(sulf2)+
  geom_col(mapping = aes(x = SampleID, y = Abundance, fill = Genus), color = "black", position = "fill", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors1) +
  scale_x_discrete(
    breaks = b,
    labels = la,
    limits = li)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))


#save table of wastewater organisms present in data
wworg <-sulf2 %>%
  select(c("SampleID", "Phylum", "Class", "Order","Family", "Genus", "Abundance"))%>%
  filter(Abundance > 0.01)
head(wworg)
dim(wworg)

as.csv(wworg, "/home/lgschaer/old/wastewater_analysis/wworg.csv")

#--------------------------------------------------------#
#EXPLORE VARIATION IN TAXA THROUGHOUT WW TREATMENT PROCESS

#Summarize abundance of each class
classabundance <- justbacteria %>%
  tax_glom(taxrank = "Class") %>%                      # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format
  arrange(Class) 
head(classabundance)

#Select and summarize necessary variables
all <- classabundance %>%
  select(Phylum, Class, SampleID, Abundance) %>%
  mutate(
    Phylum = as.character(Phylum),   #change Phylum and Class cols to character vectors
    Class = as.character(Class),
    Taxa = ifelse(Phylum == "Proteobacteria", Class, Phylum),     #Create "Taxa" column which shows Class for Proteobacteria, Phylum for all other phyla
    Taxa.1p = ifelse(Abundance < 0.01, "<1%", Taxa),              #Label taxa present at low abundance "< 1%" for both Taxa and Class
    Class.1p = ifelse(Abundance < 0.01, "<1%", Class)
  ) %>%
  arrange(desc(Abundance, SampleID))                              #Arrange with descending abundances
head(all)

#MAKING A TAXA PLOT BY CLASS - this is what we used!

all2 <- all %>%
  mutate(
    Class.1p = ifelse(Abundance < 0.01, "<1%", Class),
    Class = Class.1p
    )

#make vector for breaks and labels
b <- c("w-1-75mL", "w-2-300mL", "w-3-200mL", "w-4-250mL", "C1", "C2", "C3")

li <- c("w-1-75mL", "w-2-300mL", "w-3-200mL", "w-4-250mL", "C1", "C2", "C3")

la <- c("w-1-75mL" ="W1", "w-2-300mL" = "W2", "w-3-200mL" = "W3", "w-4-250mL" = "W4", "C1" = "C1", "C2" = "C2", "C3" = "C3")

#make color vector... harder than you'd think!

#colors8 <- c("black", "darkgray", "blue", "chocolate", "green","gold","purple","cyan","black","firebrick","lightgoldenrod","orange","blue","orangered","magenta","orange","olivedrab3","navy","seagreen","royalblue","red","black","white")

#colors9 <- c(
#  "black","#CBD588", "darkgrey", "#5F7FC7", "orange","#DA5724", "purple", "#508578", "cyan", "#CD9BCD",
#  "#AD6F3B", "green", "orchid", "yellow", "#652926", "#34eaff", 
#  "#8569D5", "#5E738F","#D1A33D", "#8A7C64","blue","darkred", "white"
#)

colors10 <- c(
  "black",   "darkcyan",     "orchid1",   "green",       "coral1",        "blue",
  "grey47",  "cyan",         "yellow",    "darkgreen",   "palegoldenrod",     "tan4",
  "grey77",  "darkblue",     "orange",    "red",         "mediumpurple1",    "purple4",
  "white",   "dodgerblue",   "firebrick", "yellowgreen", "magenta"
)  


#plot
ggplot(all2)+
  geom_col(mapping = aes(x = SampleID, y = Abundance, fill = Class), color = "black", position = "fill", show.legend = TRUE)+
  ylab("Proportion of Community") +
  scale_fill_manual(values = colors10) +
  scale_x_discrete(
    breaks = b,
    labels = la,
    limits = li)+
  xlab(NULL)+
  theme_minimal()+
  theme(axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))

#save CSV of abundance of taxa throughout ww treatment process
head(all2)

csv <- all2 %>%
  select(SampleID, Phylum, Class, Abundance) %>%
  filter(Abundance > 0.01)%>%
  arrange(desc(Abundance))
head(csv)

write_csv(csv, "/home/lgschaer/old/wastewater_analysis/taxa_abundance.csv")

csv2 <- all2 %>%
  select(SampleID, Phylum, Abundance) %>%
  filter(Abundance > 0.01) %>%
  group_by(SampleID, Phylum) %>%
  summarize(
    SumAb = sum(Abundance)
  )
head(csv2)

View(csv2)

write_csv(csv2, "/home/lgschaer/old/wastewater_analysis/taxa_abundance2.csv")


#make table comparing abundance to previous sample, comparable by Phylum and Taxa variables

pSamp <- all %>%
  select(SampleID, Phylum, Taxa, Abundance) %>%                     #select rows needed
  mutate(                                                           #rename SampleID and Abundance
    First = SampleID,
    FirstAb = Abundance,
    NextA = ifelse(SampleID == "w-1-75mL", "w-2-300mL", SampleID),  #Next 6 lines: create "SampleID" column representing the next sample in sequence using ifelse arguments
    NextB = ifelse(SampleID == "w-2-300mL", "w-3-200mL", NextA),
    NextC = ifelse(SampleID == "w-3-200mL", "w-4-250mL", NextB),
    NextD = ifelse(SampleID == "w-4-250mL", "C1", NextC),
    NextE = ifelse(SampleID == "C1", "C2", NextD),
    SampleID = ifelse(SampleID == "C2", "C3", NextE)
  ) %>%
  select(First, FirstAb, Phylum, Taxa, SampleID) %>%               #select needed variables
  left_join(all, by = c("SampleID", "Phylum", "Taxa")) %>%         #join "all" data frame by colums
  mutate(
    Next = SampleID,                                               #rename joined variables
    NextAb = Abundance
  ) %>%
  select(First, FirstAb, Next, NextAb, Phylum, Taxa)%>%            #select necessary variables
  filter(First != "C3") %>%                                        #remove samples that have C3 as "First" since C3 is last in sequence and not needed for comparison
  unite(Samples, First, Next, sep = "_") %>%                       #merge first and next into a column separated with "_"
  group_by(Samples, Taxa) %>%                                      #group by variables we will be using to plot data
  mutate(
    DifAb = mean(FirstAb) - mean(NextAb)  )                        #calculate difference in abundances between compared samples
head(pSamp)
#View(pSamp)

#preparing to plot:

#make vector for breaks and labels
breaks <- c("w-1-75mL_w-2-300mL", "w-2-300mL_w-3-200mL", "w-3-200mL_w-4-250mL", "w-4-250mL_C1", "C1_C2", "C2_C3")

limits <- c("w-1-75mL_w-2-300mL", "w-2-300mL_w-3-200mL", "w-3-200mL_w-4-250mL", "w-4-250mL_C1", "C1_C2", "C2_C3")

labels <- c("w-1-75mL_w-2-300mL" ="W1-W2", "w-2-300mL_w-3-200mL" = "W2-W3", "w-3-200mL_w-4-250mL" = "W3-W4", "w-4-250mL_C1" = "W4-C1", "C1_C2" = "C1-C2", "C2_C3" = "C2-C3")

#make colorpalette
colors2 <- colorRampPalette(brewer.pal(8, "Set1"))(23)

#by taxa (Phylum, Proteobacteria as Class)
ggplot(pSamp, mapping = aes(x = Samples, y = DifAb))+
  geom_col(aes(fill = Taxa), position = "dodge", show.legend = TRUE, na.rm = TRUE)+
  ylab("Difference in Abundance") +
  scale_fill_manual(values = colors2) +
  scale_x_discrete(
    breaks = breaks,
    labels = labels,
    limits = limits)+
  xlab(NULL)+
  theme_minimal()+
  theme(legend.position = "bottom",
        axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 20),
        title = element_text(size = 25))

#Make table to compare differences in abundance by phylum, taxa and class
head(all)

pSamp2 <- all %>%
  select(SampleID, Phylum, Taxa.1p, Class.1p, Abundance) %>%
  mutate(
    First = SampleID,
    FirstAb = Abundance,
    NextA = ifelse(SampleID == "w-1-75mL", "w-2-300mL", SampleID),
    NextB = ifelse(SampleID == "w-2-300mL", "w-3-200mL", NextA),
    NextC = ifelse(SampleID == "w-3-200mL", "w-4-250mL", NextB),
    NextD = ifelse(SampleID == "w-4-250mL", "C1", NextC),
    NextE = ifelse(SampleID == "C1", "C2", NextD),
    SampleID = ifelse(SampleID == "C2", "C3", NextE)
  ) %>%
  select(First, FirstAb, Phylum, Taxa.1p, Class.1p, SampleID) %>%
  left_join(all, by = c("SampleID", "Phylum", "Taxa.1p", "Class.1p")) %>%
  mutate(
    Next = SampleID,
    NextAb = Abundance
  ) %>%
  select(First, FirstAb, Next, NextAb, Phylum, Taxa.1p, Class.1p)%>%
  filter(First != "C3") %>%
  unite(Samples, First, Next, sep = "_") %>%
  group_by(Samples, Phylum, Class.1p) %>%
  mutate(
    DifAb = mean(FirstAb) - mean(NextAb)) #%>%
  #filter(DifAb > 0.02 | DifAb < -0.02)
head(pSamp2)
View(pSamp2)

#plotting taxa vs change in abundances!

#make color palette
colors3 <- colorRampPalette(brewer.pal(8, "Set1"))(23)

#PLOT 1

#add column with "Samples" as a factor
pSamp3 <- pSamp2 %>%
  mutate(Samples2 = factor(Samples, levels=limits)) %>%   #using limits vector from earlier
  filter(DifAb > 0.0001 | DifAb < -0.0001)
head(pSamp3)

#plot differences with stacked bar chart
ggplot(pSamp3, mapping = aes(x = Phylum, y = DifAb))+
  geom_col(aes(fill = Class.1p), position = "stack", show.legend = TRUE, na.rm = TRUE)+
  facet_wrap( ~ Samples2, nrow = 3, labeller = labeller(Samples2 = labels))+
  geom_hline(yintercept = 0, color = "black")+
  ylab("Difference in Abundance") +
  scale_fill_manual(values = colors3) +
  xlab(NULL)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 20, face = "bold"),
        legend.position = "right",
        axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 25))

#PLOT 2

#make color palette
colors4 <- colorRampPalette(brewer.pal(8, "Set1"))(17)

#add column with "Samples" as a factor
pSamp4 <- pSamp2 %>%
  mutate(Samples2 = factor(Samples, levels=limits)) %>%   #using limits vector from earlier
  filter(DifAb > 0.005 | DifAb < -0.005)
head(pSamp4)

#plot differences with stacked bar chart
ggplot(pSamp4, mapping = aes(x = Phylum, y = DifAb))+
  geom_col(aes(fill = Class.1p), position = "stack", show.legend = TRUE, na.rm = TRUE)+
  facet_wrap( ~ Samples2, nrow = 3, labeller = labeller(Samples2 = labels))+
  geom_hline(yintercept = 0, color = "black")+
  ylab("Difference in Abundance") +
  scale_fill_manual(values = colors4) +
  xlab(NULL)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 20, face = "bold"),
        legend.position = "right",
        axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 25))

#PLOT 3

#filtering a different amount
#add column with "Samples" as a factor
pSamp5 <- pSamp2 %>%
  mutate(Samples2 = factor(Samples, levels=limits)) %>%   #using limits vector from earlier
  filter(DifAb > 0.01 | DifAb < -0.01)
head(pSamp3)

#make color palette
colors5 <- colorRampPalette(brewer.pal(8, "Set1"))(17)

#plot differences with stacked bar chart
ggplot(pSamp5, mapping = aes(x = Phylum, y = DifAb))+
  geom_col(aes(fill = Class.1p), color = "black", position = "stack", show.legend = TRUE, na.rm = TRUE)+
  facet_wrap( ~ Samples2, nrow = 3, labeller = labeller(Samples2 = labels))+
  geom_hline(yintercept = 0, color = "black")+
  ylab("Difference in Abundance") +
  scale_fill_manual(values = colors8) +
  xlab(NULL)+
  theme_bw()+
  theme(strip.text.x = element_text(size = 20, face = "bold"),
        legend.position = "right",
        axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1),
        axis.title.y = element_text(size = 20),
        legend.text = element_text(size = 15),
        title = element_text(size = 25))


#experimenting..
pSamp6 <- pSamp3 %>%
  mutate(
    XEND = ifelse(abs(DifAb) >= 0.1, FirstAb, NA),
    Label = ifelse(abs(DifAb) >= 0.1, abs(DifAb), NA),
    YEND = ifelse(XEND >= 0, NextAb+DifAb, NA),
    Class = Class.1p
  )
head(pSamp5)
#View(pSamp5)

#make color palette
#colors6 <- colorRampPalette(brewer.pal(8, "Set1"))(17)

#colors7 <- c("black","dodgerblue2", "chocolate", "green", "lightgray", "yellow","purple","lightgoldenrod","darkgreen","orange","blue","orangered","magenta","dimgrey","olivedrab3","navy","seagreen","royalblue","red","black","white")

#making sure all phyla in this figure are the same color as the other figure
colors11 <- c(
  "black",      "darkcyan",   "orchid1",      "green",         "blue",
  "grey47",     "cyan",       "darkgreen",    "palegoldenrod", "tan4",
  "darkblue",   "orange",     "mediumpurple1","purple4",
  "dodgerblue", "firebrick",  "magenta"
)  

#install.packages("ggrepel")
#library(ggrepel)

#plot differences with regression chart
ggplot(pSamp6, mapping = aes(x = FirstAb, y = NextAb))+
  geom_abline(color = "black")+
  geom_point(aes(fill = Class), shape = 21, size = 6, show.legend = TRUE, na.rm = TRUE)+
  facet_wrap( ~ Samples2, nrow = 3, labeller = labeller(Samples2 = labels))+
  geom_segment(aes(xend = XEND, yend = YEND), linetype = "dashed", na.rm = TRUE) +
  geom_text(aes(label = round(Label, digits = 2)), size = 4, hjust = 1.5, vjust = .5, na.rm = TRUE)+
  ylab("Next Sample Abundance") +
  scale_fill_manual(values = colors11) +
  xlab("First Sample Abundance")+
  theme_bw()+
  theme(strip.text.x = element_text(size = 20, face = "bold"),
        legend.position = "right",
        axis.text.y.left = element_text(size = 20),
        axis.text.x = element_text(size = 15, angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 15),
        title = element_text(size = 25))
