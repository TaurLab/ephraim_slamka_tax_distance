
library(yingtools2)
library(tidyverse)
library(phyloseq)
rm(list=ls())
phy.combined.all <- readRDS("data/phy.combined.rds")



all.levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "full_taxonomy", "taxid")
taxid.levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "taxid")
bact.levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phy.combined <- phy.combined.all %>%
  phy.collapse(taxranks=taxid.levels) %>%
  mutate(otu=taxid) %>%
  select(-taxid)

s.combined <- phy.combined %>% get.samp() %>% select(sample,source,MRN,DateCollection,Sample_ID)
s1 <- s.combined %>% rename_with(~paste0(.,"1"))
s2 <- s.combined %>% rename_with(~paste0(.,"2"))
pairs.combined <- sample_names(phy.combined) %>% 
  combn(2) %>% t() %>%
  as_tibble() %>%
  rename(sample1=V1,sample2=V2) %>%
  left_join(s1,by="sample1") %>%
  left_join(s2,by="sample2") %>%
  mutate(status=case_when(
    Sample_ID1==Sample_ID2 ~ "(4) Re-seq",
    MRN1==MRN2 & DateCollection1==DateCollection2 ~ "(3) Same subj/day",
    MRN1==MRN2 ~ "(2) Same subj",
    TRUE ~ "(1) Diff subj"
  ),
  platform=case_when(
    source1!=source2 ~ "cross",
    source1=="shotgun" & source2=="shotgun" ~ "shotgun",
    source1=="amplicon.16s" & source2=="amplicon.16s" ~ "16s"
  ),
  status2=str_glue("{status}\n[{platform}]"))
pairs.combined %>% count(status2)




view_violin <- function(dist,title="") {
  pairdata <- dist %>% get.pairwise() %>%
    inner_join(pairs.combined,by=c("sample1","sample2"))
  ggplot(pairdata,aes(x=status2,y=dist,fill=status)) + 
    geom_violin() +
    geom_boxplot(alpha=0.2) + 
    ggtitle(title) + 
    expand_limits(y=c(0,1)) + theme(legend.position="none") +
    scale_y_continuous(str_glue("{title}\ndistance")) +
    scale_x_discrete("")
}

bray.dist <- calc.distance(phy.combined,method="bray")
view_violin(bray.dist,"bray")
mean.horn.dist <- calc.distance(phy.combined,method="mean.horn")
view_violin(mean.horn.dist,"mean.horn")
unfold.horn.dist <- calc.distance(phy.combined,method="unfold.horn")
view_violin(unfold.horn.dist,"unfold.horn")























