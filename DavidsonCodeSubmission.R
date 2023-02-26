#Written in R 4.2.2
#github:https://github.com/TaurLab/ephraim_slamka_tax_distance 
#note: the packages 'yingtools2' 'ytrecipes' and 'vegan' must be installed through their github pages.
rm(list=ls())
#install recquired packages
library(ggplot2)
library(base)
library(phyloseq) #phyloseq objects house 16S data
library(grid) #for manipulating plots
library(gridExtra) #for manipulating plots
library(vegan) #calculates known distance metrics like bray-curtis
library(ggtree) #for drawing trees
library(ape) #used to manipulate trees
library(tidyverse)
library(yingtools2) #ying's suite of data tools

#note: To properly run the rest of file, simply pull the code from the github, this will allow you to access the dataset.
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

#Functions:
  # given a distance matrix, return a table of pairwise distances.
  # e.g.: dist(mtcars) %>% get.pairwise()
  get.pairwise <- function(dist,diag=TRUE) {
  mat <- as.matrix(dist,diag=diag)
  # xy <- t(combn(colnames(mat), 2))
  # data.frame(xy, dist=mat[xy]) %>% as_tibble() %>% rename(sample1=X1,sample2=X2)
  rows <- rownames(mat)
  xy <- mat %>% as.data.frame() %>%
    rownames_to_column("sample1") %>%
    pivot_longer(cols=-sample1,names_to="sample2",values_to="dist") %>%
    mutate(sample1=factor(sample1,levels=rows),
           sample2=factor(sample2,levels=rows)) %>%
    filter(as.numeric(sample1)<=as.numeric(sample2))
  xy
}


  # given a table of pairwise distances, return the distance matrix.
  # dist(mtcars) %>% get.pairwise() %>% get.dist()
  get.dist <- function(pw) {
  mat <- pw %>% pivot_wider(id_cols=sample2,names_from=sample1,values_from=dist) %>%
    column_to_rownames("sample2") %>%
    as.matrix()
  as.dist(mat,diag=TRUE)
}


  # visually display a distance matrix.
  # dist(mtcars) %>% view.dist()
  # dist(mtcars) %>% view.dist(show.numbers=TRUE)
  view.dist <- function(dist,show.numbers=FALSE) {
  pairwise <- get.pairwise(dist) %>%
    mutate(sample2=fct_rev(sample2))
  ggplot(pairwise,aes(x=sample1,y=sample2,label=pretty_number(dist),fill=dist)) +
    geom_tile(color="black") +
    {if (show.numbers) geom_text(color="yellow") else NULL} +
    expand_limits(fill=c(0,1)) +
    scale_x_discrete(position = "top") +
    theme(axis.text.x=element_text(angle=90,hjust=0))
}


  # given a phyloseq object, show stack plot and overlay with distance metric in adjacent samples.
  plot.dist <- function(phy,method="bray",sortby="sample",
                      levels=c("Superkingdom","Phylum","Class",
                               "Order","Family","Genus","Species")) {
  barwidth <- 0.5
  s <- get.samp(phy,stats=TRUE,measures="InvSimpson") %>%
    arrange(!!!syms(sortby)) %>%
    mutate(sample=fct_inorder(sample),
           xvar=as.numeric(sample))
  sample_data(phy) <- s %>% set.samp()
  otu <- phy %>%
    phy.collapse() %>%
    get.otu.melt(filter.zero=FALSE) %>%
    tax.plot(data=TRUE,label.pct.cutoff = 0.1) %>%
    arrange(Species,otu) %>%
    mutate(otu=fct_inorder(otu),
           tax.label=sub(" ","\n",tax.label))
  pal <- get.yt.palette2(otu)
  lvl.list <- map(seq_along(levels),~levels[1:.x])
  names(lvl.list) <- map(lvl.list,last)
  phy.list <- map(lvl.list,~phy.collapse(phy,taxranks=.x))
  phy.list <- c(phy.list,list(asv=phy))
  
  dist.list <- phy.list %>% map(~distance(.x,method=method))
  mat.list <- dist.list %>% map(as.matrix)
  adjsamples <- s %>%
    transmute(sample1=sample,
              sample2=lead(sample),
              xvar1=xvar,
              xvar2=lead(xvar)) %>%
    filter(!is.na(sample2)) %>%
    mutate(xmid=(xvar1+xvar2)/2,
           dist.values=map2(sample1,sample2,function(sample1,sample2) {
             lst <- map_dbl(mat.list,function(mat) {
               mat[as.character(sample1),as.character(sample2)]
             }) %>% {tibble(taxlevel=names(.),ranklevel=rev(seq_along(.)),dist=unname(.))}
             lst
           })) %>%
    unnest(cols=dist.values) %>%
    mutate(y=scales::rescale(ranklevel,to=c(1.1,2)))
  
  taxlbls <- adjsamples %>% select(taxlevel,ranklevel,y) %>%
    distinct() %>%
    mutate(x=0.5)
  otu2 <- otu %>% mutate(xvar=map(xvar,~.+c(-1,1)*barwidth/2)) %>% unnest(cols=xvar)
  ggplot() +
    geom_col(data=otu,aes(x=xvar,y=pctseqs,fill=Species),color="black",size=0.25,width=barwidth) +
    geom_area(data=otu2,aes(x=xvar,y=numseqs,fill=Species),alpha=0.35,position="fill") +
    geom_text(data=s,aes(x=xvar,y=0,label=str_glue("{sample}\nInvSimp: {pretty_number(InvSimpson)}\nnseqs={short_number(nseqs,sig.digits=1)}")),
              lineheight=0.75,hjust=1.1,angle=90) +
    geom_text(data=otu,aes(x=xvar,y=y.text,label=tax.label),angle=0,lineheight=0.75) +
    scale_fill_manual(values=pal) +
    theme(legend.position = "none") +
    geom_point(data=adjsamples,aes(x=xmid,y=y,size=dist),color="steelblue") +
    geom_text(data=adjsamples,aes(x=xmid,y=y,label=pretty_number(dist,3))) +
    geom_text(data=taxlbls,aes(x=x,y=y,label=taxlevel)) +
    scale_size_continuous(range=c(0.5,10)) +
    expand_limits(size=c(0,1),y=-0.3) +
    ggtitle(method)
}


  # this function will calculate horn (or other) distance, and then 
  # apply fn to the values.
  get.taxdist <- function(phy,fn=mean,method="horn",show.work=FALSE) {
  ranks <- rank_names(phy)
  samples <- sample_names(phy)
  # create multiple phyloseq objects, collapsed at the 
  # Superkingdom, Phylum, .... , Species level. 
  phy.levels <- ranks %>% seq_along() %>%
    map(~ranks[1:.x]) %>% map(~phy.collapse(phy,taxranks=.x)) %>%
    setNames(ranks)
  phy.levels <- c(phy.levels,list("asv"=phy))
  all.levels <- names(phy.levels)
  # calculate the distance matrix (metric=method) for each level.
  # this is a list of distance matrices.
  dist.levels <- phy.levels %>% map(~distance(.x,method=method))
  # run get.pairwise() to get a list of pairwise distances.
  pairwise.levels <- dist.levels %>% imap(~{
    newname <- str_glue("dist.{.y}")
    get.pairwise(.x) %>% rename(!!sym(newname):=dist)
  })
  pairwise.all <- pairwise.levels[[1]]
  
  for (i in seq_along(all.levels)[-1]) {
    pairwise.all <- pairwise.all %>% full_join(pairwise.levels[[i]],by=c("sample1","sample2"))
  }
  pairwise.melt <- pairwise.all %>% pivot_longer(cols=-c(sample1,sample2),
                                                 names_to="dist.type",values_to="dist")
  pairwise.calcdist <- pairwise.melt %>% group_by(sample1,sample2) %>%
    summarize(dist.list=list(setNames(dist,dist.type)),
              dist=map_dbl(dist.list,fn),
              .groups = "drop")
  if (show.work) {
    return(pairwise.calcdist)
  }
  taxdist <- get.dist(pairwise.calcdist)
  taxdist
}



  view.hclust <- function(dist,.phy=phy,title="",label.pct.cutoff=0.3) {
  # by.group <- "pt.day"
  
  hc <- hclust(dist)
  tr <- as.phylo(hc)
  
  droptips <- setdiff(tr$tip.label,sample_names(.phy))
  tr <- drop.tip(tr,droptips)
  gt <- ggtree(tr) %<+% get.samp(.phy)
  gd <- gt$data
  pad <- (max(gd$x)-min(gd$x)) * 0.025
  ylim <- range(gd$x) + c(0,pad)
  s.dups.hc <- gd %>% filter(has.sameday) %>%
    arrange(y) %>%
    mutate(lbl=ifelse(has.sameday,pt.day,NA_character_),
           row=get.row(y,y,row=pt.day,min.gap=1)) 
  map <- gd %>% filter(isTip) %>%
    select(sample=label,x,y)
  xlim <- range(map$y) + c(-0.5,0.5)
  otu <- .phy %>% get.otu.melt() %>%
    left_join(map,by="sample") %>%
    tax.plot(data=TRUE,label.pct.cutoff = label.pct.cutoff)
  pal <- get.yt.palette2(otu)
  
  g.tax <- ggplot(otu,aes(x=y,y=pctseqs,fill=Species)) + 
    geom_col(show.legend=FALSE) +
    geom_text(aes(y=y.text,label=tax.label),angle=-90) +
    scale_fill_manual(values=pal) +
    coord_cartesian(xlim=xlim,expand=FALSE) +
    theme(axis.text = element_blank(),
          axis.title.y=element_blank()) +
    xlab("Sample") 
  
  row.range <- range(s.dups.hc$row) + c(-0.5,0.5)
  g.groups <- ggplot() + 
    geom_path(data=s.dups.hc,aes(x=y,y=row,group=pt.day),color="dark gray",linetype="longdash") + 
    geom_point(data=s.dups.hc,aes(x=y,y=row)) +
    geom_segment(data=map,aes(x=y,xend=y,y=row.range[1],yend=row.range[2]),color="light gray",alpha=0.35) +
    coord_cartesian(xlim=xlim,expand=FALSE) + expand_limits(y=row.range) + theme_void() 
  
  g.hclust <- gt + 
    expand_limits(x=ylim) +
    geom_point2(data=gd,aes(subset=isTip,color=pt,size=has.sameday),alpha=0.75) +
    geom_point2(data=gd,aes(subset=isTip & has.reseq,size=has.sameday),shape=1) +
    # geom_text2(data=gd,aes(subset=isTip & has.sameday, label=pt.day.samp, size=has.sameday),angle=-90,hjust=-0.05,size=3) +
    scale_x_reverse() +
    coord_flip(ylim=xlim,expand=FALSE) +
    ggtitle(title)
  
  
  # gg.stack(g.hclust,g.groups,g.tax,heights=c(1.5,1,4),align.xlim=FALSE,as.gtable = TRUE)
  gg.stack(g.hclust,g.tax,heights=c(1.5,4),align.xlim=FALSE,as.gtable = TRUE)
}




  # most samples are from different pts
  pairs %>% count(status)

# create various distances, including the customized ------------------------------------------------



weighted.mean <- function(x) {
  weights <- length(x):1
  sum( x*weights ) / sum(weights)
}
  
# function for weighted mean without superkingdom
  
  wnsk1<-function(x){
    set<-x[-1]
    weights<-length(set):1
    sum(set*weights)/sum(weights)
  }
  
# usual distances 

dist.manhattan <- distance(phy,"manhattan")
dist.bray <- distance(phy,"bray")
dist.euclidean <- distance(phy,"euclidean")
dist.horn <- distance(phy,"horn")
dist.unifrac <- distance(phy,"unifrac")
dist.wunifrac <- distance(phy,"wunifrac")


#custom distances
  dist.taxhorn.mean <- get.taxdist(phy,fn=mean,method="horn")
  dist.taxhorn.weightedmean <- get.taxdist(phy,fn=weighted.mean,method="horn")

  #dist.taxhorn.nsk<-get.taxdist(phy,fn=nsk,method="horn")
  dist.taxhorn.wnsk1<-get.taxdist(phy,fn=wnsk1,method="horn")


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
# I elected to use violin plots at they show the biases of a metric towards certain extremes.
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
#g.v.taxhorn.nsk<-view_violin(dist.taxhorn.nsk,title="taxhorn.nsk")

#the one below is the final taxHorn
g.v.taxhorn.wnsk1<-view_violin(dist.taxhorn.wnsk1,title="taxhorn.wnsk")

vlist <- list( 
  g.v.bray,
  g.v.manhattan,
  g.v.euclidean,
  g.v.horn,
  g.v.unifrac,
  g.v.wunifrac,
  g.v.taxhorn.mean,
  g.v.taxhorn.weightedmean,
  #g.v.taxhorn.nsk,
  g.v.taxhorn.wnsk1)


do.call(grid.arrange,vlist)

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
  g.hc.unifrac,
  g.hc.wunifrac,
  g.hc.taxhorn.mean,
  g.hc.taxhorn.weightedmean),
  ncol=1,nrow=1)



ggsave("plots/hclust_compare.pdf",mg,width=20,height=12)
shell.exec(normalizePath("plots/hclust_compare.pdf"))



