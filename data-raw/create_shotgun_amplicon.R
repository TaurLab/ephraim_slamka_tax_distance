

library(yingtools2)
library(tidyverse)
library(ytpipeline)
library(phyloseq)
rm(list=ls())

load("data-raw/tax.RData",envir=shotgundata<-new.env())
load("data-raw/human_microbiota_dada2.RData",envir=bact16sdata<-new.env())

taxid.shotgun <- shotgundata$t$taxon_id %>% unique()


tax.bact16s <- bact16sdata$tax.blast.full %>% 
  group_by(otu) %>%
  slice(1) %>%
  ungroup()
taxid.16s <- tax.bact16s %>% pull(staxid)  
all.taxid <- c(taxid.shotgun,taxid.16s) %>% unique() %>% sort()
bact.ranks <- c("Domain","Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax <- tibble(taxid=all.taxid) %>%
  arrange(as.numeric(taxid)) %>%
  mutate(full_taxonomy=get_full_taxonomy(taxid)) %>% 
  get_ranks_from_full_taxonomy(taxranks=bact.ranks) %>%
  mutate(full_taxonomy=full_taxonomy_list_to_char(full_taxonomy))


otu.shotgun <- shotgundata$t %>% 
  left_join(transmute(shotgundata$s,Sample_ID,sample=basename(tsvfile)),by="Sample_ID") %>%
  rename(taxid=taxon_id,numseqs=reads) %>%
  mutate(source="shotgun",
         otu=paste0("shotgun.otu=",taxid))

otu.dict <- tibble(refseq=as.character(refseq(bact16sdata$phy.dada2)),
                   otu=taxa_names(bact16sdata$phy.dada2))
otu.16s <- bact16sdata$phy.dada2 %>% 
  # filter(sample %in% c("2435D..pool1036.2", "1252T..pool661", "1373B..pool764")) %>%
  select(sample,Sample_ID) %>%
  get.otu.melt(sample_data=TRUE,tax_data=FALSE) %>%
  left_join(transmute(tax.bact16s,otu,taxid=as.character(staxid)),by="otu") %>%
  select(-pctseqs) %>%
  mutate(source="amplicon.16s") %>%
  left_join(otu.dict,by="otu") %>%
  select(-otu) %>% rename(otu=refseq) %>%
  mutate(otu=paste0("amplicon.16s.otu=",otu))


phy.allsamps <- bind_rows(otu.shotgun,otu.16s) %>%
  get.phyloseq.from.melt(tax_ranks="taxid",
                         sample_vars=c("Sample_ID","source"),
                         sample_id="sample",
                         abundance_var="numseqs",
                         taxa_id="otu")

s.allsamps <- get.samp(phy.allsamps) %>%
  left_join(bact16sdata$samp,by="Sample_ID") %>%
  select(sample,source,Sample_ID,MRN,DateCollection,Consistency,SampleType)
sample_data(phy.allsamps) <- s.allsamps %>% set.samp()

t.allsamps <- get.tax(phy.allsamps) %>% left_join(tax,by="taxid") %>%
  select(Domain,Kingdom,Phylum,Class,Order,Family,Genus,Species,full_taxonomy,taxid)

tax_table(phy.allsamps) <- t.allsamps %>% set.tax()

save(phy.allsamps,file="data-raw/phy.allsamps.RData")

library(yingtools2)
library(tidyverse)
library(ytpipeline)
library(phyloseq)
rm(list=ls())
load("data-raw/phy.allsamps.RData")

mrn.subset <- c("35454047", "35469632", "35540296", 
                "00787607", "00396816", "35516265", 
                "35519740", "35425665", "35383160", 
                "35515917", "35537584", "35367876", 
                "38053990", "38084143", "00354766")

phy.combined <- phy.allsamps %>%
  filter(SampleType=="Stool",!is.na(MRN),
         MRN %in% mrn.subset) %>%
  filter(Domain=="Bacteria") %>%
  prune_samples(sample_sums(.)>50,.)

s.combined <- phy.combined %>% get.samp()
s.allsamps <- phy.allsamps %>% get.samp()

# save(phy.combined,file="data/phy.combined.RData")
saveRDS(phy.combined,file="data/phy.combined.rds")


if (FALSE) {
  xx <- s.combined %>%
    group_by(MRN) %>%
    mutate(has.same.pt=n()>1) %>%
    ungroup() %>%
    group_by(MRN,DateCollection,source) %>%
    mutate(has.same.day=n()>1) %>%
    ungroup() %>%
    group_by(MRN,DateCollection,Sample_ID,source) %>%
    mutate(has.same.sample.source=n()>1) %>%
    ungroup() %>%
    group_by(MRN,DateCollection,Sample_ID) %>%
    mutate(has.cross.platform=n_distinct(source)>1) %>%
    ungroup() %>%
    mutate(has.same.day.only=has.same.day & !has.same.sample.source)
  
  xx %>% count(has.same.day.only,has.same.sample.source)
  xx %>% count(has.same.day,has.same.sample.source)
  xx %>% count(has.same.pt,
               has.same.day,
               has.same.sample.source,
               has.cross.platform)
  
  
  yy <- xx %>% 
    group_by(MRN) %>%
    summarize(n.samps=n(),
              n.amplicon=sum(source=="amplicon.16s"),
              n.shotgun=sum(source=="shotgun"),
              across(.cols=c(has.same.pt,has.same.day,has.same.day.only,has.same.sample.source,has.cross.platform),
                     .fns=sum)) %>%
    ungroup() %>%
    arrange(desc(n.samps))
  
  yy %>% select(-has.same.day) %>%
    # filter(has.cross.platform>0,has.same.sample.source>0) %>%
    arrange(
      !(has.cross.platform>=6 & has.same.sample.source>0),
      !(has.same.day.only==10 & has.cross.platform==76),
      !(has.same.sample.source>=10),
      desc(n.samps)) %>%
    mutate(row=row_number()) %>% 
    # slice(1:15) %>% pull(MRN) %>% copy.as.Rcode()
    dt()  
  
  
  look <- function(mrn) {
    s <- phy.combined %>%
      filter(MRN==mrn) %>%
      get.samp()
    otu <- phy.combined %>%
      filter(MRN==mrn) %>%
      phy.collapse.bins() %>%
      get.otu.melt() %>%
      mutate(sample=fct_reordern(sample,MRN,DateCollection))
    ggplot() +
      geom_taxonomy(data=otu,aes(x=sample,y=pctseqs,fill=Species,label=Species)) +
      geom_text(data=s,aes(x=sample,y=1,label=nseqs),angle=-90,color="blue") +
      theme(axis.text.x=element_text(angle=-90))
  }  

  gg <- mrn.subset %>% map(look,.progress = TRUE)
  
  mrn <- "00396816"
  look(mrn)
  gg[[2]]
  
  pdf("asdf.pdf",width=14,height=8)
  gg[1:5]
  dev.off()
  pdf("asdf2.pdf",width=14,height=8)
  gg[6:10]
  dev.off()
}













