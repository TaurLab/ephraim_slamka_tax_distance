
# load data and functions.R (run this first) ---------------------------------------------------------------

library(phyloseq) #phyloseq objects house 16S data
library(grid) #for manipulating plots
library(gridExtra) #for manipulating plots
library(vegan) #calculates known distance metrics like bray-curtis
library(ggtree) #for drawing trees
library(ape) #used to manipulate trees
library(tidyverse)
library(yingtools2) #ying's suite of data tools

rm(list=ls())
phy <- readRDS("data/mock_phylo_compact.rds")
source("R/functions.R")
s <- get.samp(phy)

# assess all pairwise combinations and classify as one of 4 groups.
s1 <- get.samp(phy) %>% select(sample,pt,pt.day,pt.day.samp) %>% rename_with(~paste0(.,"1"))
s2 <- get.samp(phy) %>% select(sample,pt,pt.day,pt.day.samp) %>% rename_with(~paste0(.,"2"))
pairs <- sample_names(phy) %>%
  combn(2) %>% t() %>%
  as_tibble() %>%
  rename(sample1=V1,sample2=V2) %>%
  left_join(s1,by="sample1") %>%
  left_join(s2,by="sample2") %>%
  mutate(status=case_when(
    pt.day.samp1==pt.day.samp2 ~ "same pt \r\n day \r\n sample",
    pt.day1==pt.day2 ~ "  same pt \r\n day \r\n (diff sample)",
    pt1==pt2 ~ "  same pt \r\n (diff day/sample) ",
    TRUE ~"  diff pt  ",
  ))
rm(s1,s2)

# most samples are from different pts
pairs %>% count(status)

# create various distances, including the customized ------------------------------------------------

 
 
weighted.mean <- function(x) {
  weights <- length(x):1
  sum( x*weights ) / sum(weights)
}

# usual distances 

dist.manhattan <- distance(phy,"manhattan")
dist.bray <- distance(phy,"bray")
dist.euclidean <- distance(phy,"euclidean")
dist.horn <- distance(phy,"horn")
dist.unifrac <- distance(phy,"unifrac")
dist.wunifrac <- distance(phy,"wunifrac")


#ignore what follows, just an attempt at learning, any mention of these functions was just an attempt at learning-----------
#this isn't finished, I'm just trying to figure out how all of this works.
athing<- function(list){sqrt(sum(list)) }
aathing<- function(list) {n<-length(list)
pwr<- 1/n
(sum(list))^pwr}

aabthing<-function(x) {n<-length(x)
pwr<- 1/n
this<-((sum(x))^pwr)
this/n}

aacthing<-function(x) {n<-length(x)
pwr<- 1/n
(sum(x))^pwr}
geomean<- function(set){n<-length(set)
pwr<-1/n
(prod(set))^pwr}
dist.taxhorn.aabthing<- get.taxdist(phy, fn=aabthing, method="horn")
g.hc.aabthing<-view.hclust(dist.taxhorn.aabthing, title= "aabthing")
grid.draw(g.hc.aabthing)

#working here
abthing<-function(x){
  n<-(sum(x))*(prod(x))
  m<-(n^(1/length(x)))/x
  mean(m)/length(m)
  #pwr<- 1/(length(n))
  # (prod(n))^pwr
}
acthing<-function(x){
  n<-(sum(x))*(prod(x))
  m<-(n^(1/length(x)))/x
  mean(m)/mean(x)
  #pwr<- 1/(length(n))
  # (prod(n))^pwr
}
#test these two ^^vv
geomean2<- function(x){
  
  n<-length(x)
  pwr<- (1/n)
  (prod(x))^pwr}

bthing<- function(x){
  n<-x^x
  m<-sum(n)
  m/length(n)
}
bathing<-function(x){
  n<-x^(1/(x^2))
  m<-sum(n)
  m/length(n)
}
bbthing<-function(x){
  n<-x^(sin(x))
  m<-sum(n)
  m/length(n)
}
bcthing<-function(x){
  n<-x^tan(x)
  m<-sum(n)
  m/length(n)
}
#this is the area I'm working on
#note for future ephraim, if elses are doomed to fail in 
#nonspecific cases. try using a median or mean or something, idk.
#be clever future E. past E believes in you.
bdthing<-function(x){
  
  d<-mean(x)
  n<-(x)^x
  o<-(x)^(1/x)
  w<-mean(n)
  s<-mean(o)
  ifelse(x>d, return(w), return(s))
}
frog<-c(0,2,44,57,0.47,9)
mean(frog[-1])

nsk<-function(x){
  mean(x[-1])
}
wnsk1<-function(x){
  set<-x[-1]
  weights<-length(set):1
  sum(set*weights)/sum(weights)
}
wnsk2<-function(x){
  set<-x[-1]
  #weights<-length(set):1
  weights<-c(2,3,4,5,4,2,1)
  sum(set*weights)/sum(weights)
}
weighted.mean <- function(x) {
  weights <- length(x):1
  sum( x*weights ) / sum(weights)
}


#custom distances
dist.taxhorn.mean <- get.taxdist(phy,fn=mean,method="horn")
dist.taxhorn.weightedmean <- get.taxdist(phy,fn=weighted.mean,method="horn")
dist.taxhorn.aabthing<- get.taxdist(phy,fn=aabthing,method="horn")
dist.taxhorn.abthing<- get.taxdist(phy,fn=abthing, method="horn")
dist.taxhorn.acthing<-get.taxdist(phy,fn=acthing,method="horn")
dist.taxhorn.bthing<-get.taxdist(phy,fn=bthing,method="horn")
dist.taxhorn.bdthing<-get.taxdist(phy,fn=bdthing,method="horn")
dist.taxhorn.bbthing<-get.taxdist(phy,fn=bbthing,method="horn")
dist.taxhorn.nsk<-get.taxdist(phy,fn=nsk,method="horn")
dist.taxhorn.wnsk2<-get.taxdist(phy,fn=wnsk2,method="horn")

dist.taxhorn.wnsk1<-get.taxdist(phy,fn=wnsk1,method="horn")
dist.taxhorn.bcthing<-get.taxdist(phy,fn=bcthing,method="horn")
# dist.taxhorn.scratch <- get.taxdist(phy,fn=mean,method="horn",show.work = TRUE)
# x <- dist.taxhorn.scratch$dist.list[[10]]
# mean(x)
# weighted.mean(x)



# examine each pt's samples at different tax levels -------------------------------------------------

# get a subset phyloseq and run plot.dist. 
physub <- subset_samples(phy,pt=="HV")
g.tyler.bray <- plot.dist(physub,"bray")
g.tyler.bray
g.tyler.horn <- plot.dist(physub,"horn")
g.tyler.horn

# to plot this for all subjects, create list of ggplots:
glist <- get.samp(phy) %>% group_by(pt) %>%
  group_map(~{
    pt <- .x$pt[1]
    metric <- "horn"
    physub <- prune_samples(.x$sample,phy) %>% prune_unused_taxa()
    plot.dist(physub, method=metric, sortby=c("pt","pt.day","pt.day.samp")) +
      ggtitle(str_glue("subject: {pt} / metric: {metric}"))
  },.keep=TRUE)

pdf("plots/mock_metrics.pdf",width=34,height=17)
glist
dev.off()
shell.exec(normalizePath("plots/mock_metrics.pdf"))




# view some of the pairs to get an idea ----------------------------------------------------------


#select 10 pairs from each comparison and show
otu <- phy %>% get.otu.melt() %>%
  tax.plot(data=TRUE)
pal <- get.yt.palette2(otu)

pairs.select <- pairs %>% 
  group_by(status) %>% slice(1:10) %>% 
  transmute(pair.id=row_number(),status,sample1,sample2) %>% 
  ungroup() %>% 
  pivot_longer(cols=c(sample1,sample2),values_to="sample") %>%
  left_join(otu,by="sample")

g.pairs.select <- 
  ggplot(pairs.select,aes(x=name,y=pctseqs,fill=Species)) + geom_col(show.legend=FALSE) +
  # geom_text(aes(y=y.text,label=tax.label),angle=-90) +
  scale_fill_manual(values=pal) +
  facet_grid(status~pair.id,scales="free_x",space="free_x")
g.pairs.select

pdf("plots/pairs.grouping.taxview.pdf",width=20,height=12)
g.pairs.select
dev.off()

shell.exec(normalizePath("plots/pairs.grouping.taxview.pdf"))


# view all resequenced samples ----------------------------------------------------------

s.reseq <- s %>% group_by(pt.day.samp) %>%
  filter(n()>1)

otu.reseq <- phy %>% prune_samples(s.reseq$sample,.) %>%
  get.otu.melt() %>%
  tax.plot(data=TRUE)
pal <- get.yt.palette2(otu.sameday)

g.reseq.sample <- ggplot(otu.reseq,aes(x=sample,y=pctseqs,fill=Species)) +
  geom_col(show.legend = FALSE) + 
  geom_text(aes(y=y.text,label=tax.label),angle=-90) +
  scale_fill_manual(values=pal) +
  facet_grid(.~pt.day.samp,scales="free_x",space="free_x")

g.reseq.sample

ggsave("plots/g.all.reseq.pdf",g.reseq.sample,width=20,height=10)
shell.exec(normalizePath("plots/g.all.reseq.pdf"))



# view all same day samples (not same sample)-----------------------------------------------

s.sameday <- s %>% 
  group_by(pt.day.samp) %>%
  slice(1) %>%
  group_by(pt.day) %>%
  filter(n()>1) %>%
  ungroup()

otu.sameday <- phy %>% prune_samples(s.sameday$sample,.) %>%
  get.otu.melt() %>%
  tax.plot(data=TRUE)
pal <- get.yt.palette2(otu.sameday)

g.sameday.sample <- ggplot(otu.sameday,aes(x=sample,y=pctseqs,fill=Species)) +
  geom_col(show.legend = FALSE) + 
  geom_text(aes(y=y.text,label=tax.label),angle=-90) +
  scale_fill_manual(values=pal) +
  facet_grid(.~pt.day,scales="free_x",space="free_x")
g.sameday.sample
ggsave("plots/g.same.sample.pdf",g.sameday.sample,width=15,height=8)
shell.exec("plots/g.same.sample.pdf")


# violin of distances -----------------------------------------------------

#compare by group:


view_violin <- function(dist,title) {
  pairdata <- dist %>% get.pairwise() %>%
    inner_join(pairs,by=c("sample1","sample2"))
  ggplot(pairdata,aes(x=status,y=dist,fill=status)) + geom_violin(position = position_dodge(10)) +
    geom_boxplot(alpha=0.1) + ggtitle(title) + expand_limits(y=c(0,1)) +
    theme(legend.position= "none") 
}

g.v.manhattan <- view_violin(dist.manhattan,title="manhattan")
g.v.bray <- view_violin(dist.bray,title="bray")
g.v.euclidean <- view_violin(dist.euclidean,title="euclidean")
g.v.horn <- view_violin(dist.horn,title="horn")
g.v.unifrac <- view_violin(dist.unifrac,title="unifrac")
g.v.wunifrac <- view_violin(dist.wunifrac,title="wunifrac")
g.v.taxhorn.mean <- view_violin(dist.taxhorn.mean,title="taxhorn.mean")
g.v.taxhorn.weightedmean <- view_violin(dist.taxhorn.weightedmean,title="taxhorn.weightedmean")
#g.v.taxhorn.aabthing <- view_violin(dist.taxhorn.aabthing, title="taxhorn.aabthing")
#g.v.taxhorn.abthing<-view_violin(dist.taxhorn.abthing,title="taxhorn.abthing")
#g.v.taxhorn.acthing<-view_violin(dist.taxhorn.acthing,title="taxhorn.acthing")
#g.v.taxhorn.bthing<-view_violin(dist.taxhorn.bthing,title="taxhorn.bthing")
#g.v.taxhorn.bdthing<-view_violin(dist.taxhorn.bdthing,title="taxhorn.bdthing")
#g.v.taxhorn.bbthing<-view_violin(dist.taxhorn.bbthing,title="taxhorn.bbthing")
#g.v.taxhorn.wnsk2<-view_violin(dist.taxhorn.wnsk2,title="taxhorn.wnsk2")

g.v.taxhorn.wnsk1<-view_violin(dist.taxhorn.wnsk1,title="taxhorn.wnsk")
#g.v.taxhorn.nsk<-view_violin(dist.taxhorn.nsk,title="taxhorn.nsk")
#g.v.taxhorn.bcthing<-view_violin(dist.taxhorn.bcthing,title="taxhorn.bcthing")
vlist <- list( 
  g.v.bray,
  g.v.manhattan,
  g.v.euclidean,
  g.v.horn,
  g.v.unifrac,
  g.v.wunifrac,
  g.v.taxhorn.mean,
  #g.v.taxhorn.weightedmean,
  #g.v.taxhorn.aabthing,
  #g.v.taxhorn.abthing,
  #g.v.taxhorn.acthing,
  #g.v.taxhorn.bthing,
  #g.v.taxhorn.bdthing,
  #g.v.taxhorn.bbthing,
  #g.v.taxhorn.bcthing),
  #g.v.taxhorn.nsk,
  g.v.taxhorn.wnsk1)
  #g.v.taxhorn.wnsk2)
  

do.call(grid.arrange,vlist)
#find way to seperate violin plot collumns


mg <- marrangeGrob(vlist,ncol=3,nrow=3)
ggsave("plots/violin_compare_groups.pdf",mg,width=20,height=12)
shell.exec(normalizePath("plots/violin_compare_groups.pdf"))


# pca of samples ---------------------------------------------------------------------

view.pca <- function(metric) {
  # metric <- "euclidean"
  pcadata <- ordinate(phy, method = "NMDS", distance = metric) %>%
    scores(display = "sites") %>%
    as.data.frame() %>%
    rownames_to_column("sample") %>%
    left_join(s, by = "sample")
  s.dups <- pcadata %>% filter(has.sameday)
  g.pca <- ggplot(pcadata, aes(x = NMDS1, y = NMDS2)) +
    geom_path(data=s.dups,aes(group=pt.day),color="dark gray",size=1) +
    geom_point(aes(color=pt),alpha = 0.75) +
    # geom_text(data=s.dups,aes(label=pt.day)) +
    theme(aspect.ratio = 1) + ggtitle(metric)
  g.pca
}

view.pca("bray")
view.pca("horn")

metrics <- c("euclidean","bray","unifrac","wunifrac","horn")
g.pca.list <- metrics %>% map(view.pca)
do.call(grid.arrange,c(g.pca.list,list(nrow=2)))


# hclust of samples -------------------------------------------------------


g.hc.wnsk<-view.hclust(dist.taxhorn.wnsk,title="taxhorn.wnsk")
g.hc.manhattan <- view.hclust(dist.manhattan,title="manhattan")
g.hc.bray <- view.hclust(dist.bray,title="bray")
g.hc.euclidean <- view.hclust(dist.euclidean,title="euclidean")
g.hc.horn <- view.hclust(dist.horn,title="horn")
g.hc.unifrac <- view.hclust(dist.unifrac,title="unifrac")
g.hc.wunifrac <- view.hclust(dist.wunifrac,title="wunifrac")
g.hc.taxhorn.mean <- view.hclust(dist.taxhorn.mean,title="taxhorn.mean")
g.hc.taxhorn.weightedmean <- view.hclust(dist.taxhorn.weightedmean,title="taxhorn.weightedmean")
  
grid.newpage()
grid.draw(g.hc.manhattan)

# g.hc.taxhorn.weightedmean <- view.hclust(dist.taxhorn.weightedmean,title="taxhorn.weightedmean",label.pct.cutoff = 0.1)
# grid.draw(g.hc.taxhorn.weightedmean)

mg <- marrangeGrob(list(
  g.hc.wnsk,
  g.hc.bray,
  g.hc.manhattan,
  g.hc.euclidean,
  g.hc.horn,
  g.hc.aabthing,
  g.hc.unifrac,
  g.hc.wunifrac,
  g.hc.taxhorn.mean,
  g.hc.taxhorn.weightedmean),
  ncol=1,nrow=1)



ggsave("plots/hclust_compare.pdf",mg,width=20,height=12)
shell.exec(normalizePath("plots/hclust_compare.pdf"))
