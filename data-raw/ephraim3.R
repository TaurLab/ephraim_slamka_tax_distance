
library(yingtools2)
library(ytrecipes)
library(tidyverse)
library(phyloseq)
library(ape)
library(ggtree)
library(scales)
library(grid)
library(gridExtra)
rm(list=ls())


load("tyler.phylo.RData")

plot.pt2 <- function(phy,method="bray") {
  ssub <- phy %>% get.samp() %>%
    arrange(sample) %>%
    mutate(sample=factor(sample,levels=sample),
           x=as.numeric(sample))
  sample_data(phy) <- ssub %>% set.samp()
  
  otusub <- phy %>%
    get.otu.melt() %>%
    tax.plot(data=TRUE)
  pal <- get.yt.palette2(otusub)
  
  dist <- distance(phy,method=method)
  mat <- as.matrix(dist)

  adjpairs <- ssub %>%
    transmute(sample1=sample,
              sample2=lead(sample)) %>%
    filter(!is.na(sample2)) %>%
    mutate(dist=map2_dbl(sample1,sample2,~mat[.x,.y]))
  otusub <- otusub %>% mutate(sample=as.numeric(sample))
  adjpairs <- adjpairs %>% mutate(sample1=as.numeric(sample1),
                                  sample2=as.numeric(sample2),
                                  x.mid=(sample1+sample2)/2)

  g.tax <- ggplot() +
    geom_col(data=otusub,aes(x=x,y=pctseqs,fill=Species),width=0.95) +
    scale_fill_manual(values=pal) +
    geom_text(data=otusub,aes(x=x,y=y.text,label=tax.label),angle=-90) +
    # geom_text(data=ssub,aes(x=sample,y=0,label=sample),vjust=1) +
    geom_col(data=ssub,aes(x=x,y=1),width=0.95,color="black",fill=NA) +
    scale_y_continuous("Relative Abundance",label=percent) +
    # scale_x_continuous("Transplant Day") +
    theme(legend.position="none") +
    geom_point(data=adjpairs,aes(x=x.mid,y=1.1,size=dist),color="steelblue") +
    geom_text(data=adjpairs,aes(x=x.mid,y=1.1,label=short_number(dist))) +
    # geom_segment(data=adjpairs,aes(x=sample1+0.05,xend=bmtday2-0.05,y=1.07,yend=1.07)) +
    expand_limits(size=0) +
    # scale_size_continuous(range=c(1,10)) +
    ggtitle(method)
  
  
  g.tax
  
}


methods <- c("morisita","canberra","gower","altGower","horn","bray","euclidean","kulczynski","jaccard","manhattan",
             "unifrac","wunifrac","dpcoa")

plot.pt2(tyler.phy,"bray")
plot.pt2(tyler.phy,"unifrac")

glist <- methods %>% map(~plot.pt2(tyler.phy,.))

pdf("tyler_dists.pdf",width=15,height=8)
glist
dev.off()


shell.exec("tyler_dists.pdf")
















