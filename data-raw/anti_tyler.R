
##create mock sequences that are mod of original 


library(yingtools2)
library(phyloseq)
library(tidyverse)
library(ggtree)
library(scales)
library(grid)
library(gridExtra)

rm(list=ls())
load("C:\\Users\\taury\\Desktop\\tyler.phylo.RData")

otu <- get.otu.melt(tyler.phy)
pal <- get.yt.palette2(otu)

ggplot(otu,aes(x=sample,y=pctseqs,fill=Species)) +
  geom_col() +
  theme(legend.position="none") +
  scale_fill_manual(values=pal)


# dist <- distance(tyler.phy,"bray")
# dist
# 
# mat <- as.matrix(dist)
# xy <- colnames(mat) %>% combn(2) %>% t()
# pwdist <- data.frame(xy, dist=mat[xy]) %>%
#   rename(sample1=X1,sample2=X2,dist=dist)



# asdf --------------------------------------------------------------------

tyler1 <- tyler.phy
tyler2 <- tyler.phy

taxa_names(tyler1) <- refseq(tyler.phy) %>% as.character() %>% unname() %>% seq_along() %>% paste0("asv",.,"_a")
taxa_names(tyler2) <- refseq(tyler.phy) %>% as.character() %>% unname() %>% seq_along() %>% paste0("asv",.,"_b")
sample_names(tyler1) <- paste0(sample_names(tyler1),"_a")
sample_names(tyler2) <- paste0(sample_names(tyler2),"_b")

tyler.combined <- merge_phyloseq(otu_table(tyler1),otu_table(tyler2),
                                 tax_table(tyler1),tax_table(tyler2),
                                 refseq(tyler1),refseq(tyler2),
                                 sample_data(tyler1),sample_data(tyler2))

s.combined <- tyler.combined %>% get.samp() %>%
  mutate(group=ifelse(grepl("_b$",sample),"A","B"))


plot.pt2(tyler.combined,"morisita")

plot.pt2(tyler.combined,"bray")







