library(yingtools2)
library(tidyverse)
library(phyloseq)
rm(list=ls())
phy.combined.all <- readRDS("data/phy.combined.rds")



# amplicon 16S
# -PCR-amplify (make copies) of the 16S region
# -needs ~ 50K seqs
# -sequence all the 16S
# -denoising (get rid of errors and group together into ASVs)
# -classify seqs
# -levels: Domain, Kingom, ..., Genus, Species, ASV.

# shotgun
# -'bag of genes' - various regions from various bacteria
# -infer amounts and classification in a fuzzy logic way
# -needs ~ 1 million+
# -unlike amplicon, can get functional traits, abx-resistance, etc.
# -every little bit is classified (sort of like blast)
# -levels: Domain, Kingom, ..., Genus, Species, taxids

# amplicon 16S: to make comparable,
# change levels: Domain, Kingom, ..., Genus, Species, taxid


all.levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "full_taxonomy", "taxid")
taxid.levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "taxid")
bact.levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phy.combined <- phy.combined.all %>%
  phy.collapse(taxranks=taxid.levels) %>%
  mutate(otu=taxid) %>%
  select(-taxid)

s <- phy.combined %>% get.samp()
s %>% count(source)
s %>% group_by(Sample_ID) %>%
  dt()
s %>% filter(MRN=="35540296")  %>% count(source)


# old way
# s <- phy.combined %>% get.samp() %>% 
#   filter(MRN=="35540296")
# phy.onept <- prune_samples(s$sample,phy.combined)

# new way
phy.pt <- phy.combined %>%
  filter(MRN=="35540296")
otu <- phy.pt %>%
  get.otu.melt()
s <- phy.pt %>% get.samp()



# show all
ggplot() +
  geom_taxonomy(data=otu,aes(x=sample,y=pctseqs,fill=otu,label=Species)) +
  # geom_text(data=s,aes(x=sample,y=1,label=source),color="blue",angle=90,size=3) +  
  theme(axis.text.x=element_text(angle=90))

# separate by facet
ggplot(otu,aes(x=Sample_ID,y=pctseqs,fill=otu,label=Species)) +
  geom_taxonomy(position="fill") + 
  facet_grid(source ~ .) +
  theme(axis.text.x=element_text(angle=90))

# show number of seqs
ggplot() +
  expand_limits(y=1.2) + 
  geom_taxonomy(data=otu,aes(x=Sample_ID,y=pctseqs,fill=otu,label=Species),
                position="fill") +
  geom_text(data=s,aes(x=Sample_ID,y=1,label=source),color="blue",angle=90,size=3) +
  facet_grid(source ~ .) +
  theme(axis.text.x=element_text(angle=90))







