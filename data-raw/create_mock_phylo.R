
}library(tidyverse)
library(yingtools2)
library(ytdata)
library(phyloseq)

setwd("C:/Users/Ying/Desktop/ephraim")
rm(list=ls())
load("tyler.phylo.RData")
load("phy.bmt.all.RData")

# get rid of duplication from BMT.
pt.single.bmt <- phy.bmt.all %>%
  get.samp() %>%
  group_by(MRN,BMT) %>%
  summarize(median.bmtday=median(bmtday),
            .groups="drop") %>%
  group_by(MRN) %>%
  slice(which.min(abs(median.bmtday))) %>%
  ungroup()

# find patients with any duplicates
s <- phy.bmt.all %>% get.samp() %>%
  filter(nseqs>5000) %>%
  filter(SampleType=="Stool") %>%
  inner_join(pt.single.bmt,by=c("MRN","BMT")) %>%
  group_by(MRN,BMT) %>%
  mutate(nsamps.pt.total=n()) %>%
  ungroup() %>%
  group_by(Sample_ID) %>%
  mutate(has.reseq=n()>1) %>%
  ungroup() %>%
  group_by(MRN,DateCollection) %>%
  mutate(has.sameday=n()>1 & !has.reseq) %>%
  ungroup()
s.duplicates <- s %>% group_by(MRN,BMT) %>%
  filter(sum(has.sameday+has.reseq)>1) %>% ungroup() %>%
  select(MRN,BMT,sample,Sample_ID,DateCollection,bmtday,has.sameday,has.reseq,nseqs)

bmt.phy <- prune_samples(s.duplicates$sample,phy.bmt.all) %>%
  prune_unused_taxa()
sample_data(bmt.phy) <- s.duplicates %>% set.samp()

#combine tyler and bmt dups

taxa_names(tyler.phy) <- refseq(tyler.phy) %>% as.character()
taxa_names(bmt.phy) <- refseq(bmt.phy) %>% as.character()

s.bmt <- get.samp(bmt.phy) %>% mutate(source="bmt")
s.tyler <- get.samp(tyler.phy) %>%
  mutate(source="tyler",MRN="tyler",DateCollection=as.Date("2021-07-15"),
         has.sameday=TRUE,has.reseq=TRUE)

s.mock <- bind_rows(s.bmt,s.tyler) %>%
  transmute(sample,MRN,BMT,DateCollection,
            Sample_ID=ifelse(source=="bmt",Sample_ID,letter),
            day=ifelse(source=="bmt",bmtday,1),
            has.sameday,has.reseq,source) %>%
  mutate(pt=MRN,
         pt.day=str_glue("{pt}_d{day}"),
         pt.day.samp=str_glue("{pt.day}_{Sample_ID}"))

mock.phy <- merge_phyloseq(
  otu_table(tyler.phy),
  otu_table(bmt.phy),
  tax_table(tyler.phy),
  tax_table(bmt.phy),
  refseq(tyler.phy),
  refseq(bmt.phy))

tax_table(mock.phy) <- mock.phy %>% get.tax() %>%
  select(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species) %>%
  set.tax()
sample_data(mock.phy) <- s.mock %>% set.samp()

#just to make the names simpler
sample_names(mock.phy) <- str_replace(sample_names(mock.phy),"\\.\\.BMT.*$","")

save(mock.phy,file="mock_phylo.RData")


