
# load data (run this) ---------------------------------------------------------------

library(tidyverse)
library(yingtools2)
library(phyloseq)
library(grid)
library(gridExtra)
library(vegan)
library(ggtree)
library(ape)
rm(list=ls())

phy <- readRDS("data/mock_phylo_compact.RData")

# functions (run this) ---------------------------------------------------------------

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





# practice distance manipulation ----------------------------------------------------------------

d1 <- dist(mtcars)
view.dist(d1)
pw <- get.pairwise(d1)
d2 <- get.dist(pw)
view.dist(d2)

# categorize dists between samples ----------------------------------------------------------------

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
    pt.day.samp1==pt.day.samp2 ~ "same pt/day/sample",
    pt.day1==pt.day2 ~ "same pt/day (different sample)",
    pt1==pt2 ~ "same pt (different day/sample)",
    TRUE ~ "different pt",
  ))
rm(s1,s2)



# most samples are from different pts
pairs %>% count(status)




# levels ------------------------------------------------------------------


# example of plot.dist
# tyler's samples, with bray or horn.
physub <- subset_samples(phy,pt=="tyler")
plot.dist(physub,"horn")
plot.dist(physub,"bray")

# create list of ggplots, one for each patient.
glist <- get.samp(phy) %>% group_by(pt) %>%
  group_map(~{
    pt <- .x$pt[1]
    metric <- "horn"
    physub <- prune_samples(.x$sample,phy) %>% prune_unused_taxa()
    plot.dist(physub, method=metric, sortby=c("pt","pt.day","pt.day.samp")) +
      ggtitle(str_glue("subject: {pt} / metric: {metric}"))
  },.keep=TRUE)

#create view of each patient

pdf("mock_metrics.pdf",width=34,height=17)
glist
dev.off()
shell.exec("mock_metrics.pdf")



# pca/hclust ---------------------------------------------------------------------


s.dups <- get.samp(phy) %>%
  group_by(pt.day) %>%
  mutate(dups=n()>1) %>%
  ungroup()

phy.dups <- phy
sample_data(phy.dups) <- s.dups %>% set.samp()

metric <- "bray"
pcadata <- ordinate(phy.dups, method = "NMDS", distance = metric) %>%
  scores(display = "sites") %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(s.dups, by = "sample")
g.pca <- ggplot(pcadata, aes(x = NMDS1, y = NMDS2, color=pt.day, size=dups)) +
  geom_point(alpha = 0.75,show.legend = FALSE) +
  theme(aspect.ratio = 1) + ggtitle(metric)
g.pca



dist <- distance(phy,metric)
hc <- hclust(dist)
tr <- as.phylo(hc)
gt <- ggtree(tr) %<+% get.samp(phy)
gd <- gt$data
gt + geom_point2(data=gd,aes(subset=isTip,color=pt),size=4) +
  geom_text2(data=gd,aes(subset=isTip,label=label),hjust=0)
g.hclust <- gt + geom_point2(data=gd,aes(subset=isTip,color=pt.day),size=4,show.legend=FALSE) +
  geom_text2(data=gd,aes(subset=isTip,label=label),hjust=0,show.legend=FALSE)
g.hclust


# hclust ------------------------------------------------------------------

# create a plot with hierarchical clustering and stacked taxa.
phy.hc <- phy
dist <- distance(phy.hc,"horn")
hc <- hclust(dist)
tr <- as.phylo(hc)
gt <- ggtree(tr)

# make sure plot aesthetically aligns.
xlim <- gt$data$y %>% range() + c(-0.5,0.5)
g.hc <- gt + scale_x_reverse() +
  coord_flip(ylim=xlim,expand=FALSE)
g.hc

map <- g.hc$data %>% filter(isTip) %>% select(x=y,sample=label)
s <- get.samp(phy.hc) %>%
  left_join(map,by="sample")
otu <- phy.hc %>% get.otu.melt() %>%
  tax.plot(data=TRUE) %>%
  left_join(map,by="sample")
pal <- get.yt.palette2(otu)
g.tax <- ggplot() + expand_limits(y=-0.06) +
  geom_col(data=otu,aes(x=x,y=pctseqs,fill=Species),show.legend = FALSE) +
  geom_text(data=otu,aes(x=x,y=y.text,label=tax.label),angle=-90) +
  geom_point(data=s,aes(x=x,y=-0.01,color=pt),size=4) +
  scale_fill_manual("Bacteria phylotype",values=pal) +
  coord_cartesian(xlim=xlim,expand=FALSE) +
  theme(axis.text=element_blank(),
        axis.title=element_blank(),
        axis.ticks=element_blank())

gg.stack(g.hc,g.tax,heights=c(2,5),align.xlim=FALSE)

pdf("hclust.pdf",width=20,height=10)
gg.stack(g.hc,g.tax,align.xlim=FALSE,heights=c(2,5),newpage = FALSE)
dev.off()

shell.exec('hclust.pdf')


# hclust 2 ----------------------------------------------------------------

# function to view various distances
view.hclust <- function(dist,phy,title="hclust") {
  hc <- hclust(dist)
  tr <- as.phylo(hc)
  gt <- ggtree(tr)
  
  xlim <- gt$data$y %>% range() + c(-0.5,0.5)
  g.hc <- gt + scale_x_reverse() + coord_flip(ylim=xlim,expand=FALSE) +
    ggtitle(title)
  map <- g.hc$data %>% filter(isTip) %>% select(x=y,sample=label)
  s <- get.samp(phy) %>%
    left_join(map,by="sample")
  otu <- phy %>% get.otu.melt() %>%
    tax.plot(data=TRUE) %>%
    left_join(map,by="sample")
  pal <- get.yt.palette2(otu)
  g.tax <- ggplot() + expand_limits(y=-0.06) +
    geom_col(data=otu,aes(x=x,y=pctseqs,fill=Species),show.legend = FALSE) +
    geom_text(data=otu,aes(x=x,y=y.text,label=tax.label),angle=-90) +
    geom_point(data=s,aes(x=x,y=-0.01,color=pt),size=4) +
    scale_fill_manual("Bacteria phylotype",values=pal) +
    coord_cartesian(xlim=xlim,expand=FALSE) +
    theme(axis.text=element_blank(),
          axis.title=element_blank(),
          axis.ticks=element_blank())
  gg.stack(g.hc,g.tax,heights=c(2,5),align.xlim=FALSE)
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
  # pairwise.calcdist <- pairwise.melt %>% group_by(sample1,sample2) %>%
  #     summarize(dist=fn(dist),.groups="drop")
  taxdist <- get.dist(pairwise.calcdist)
  taxdist
  
}



# to make a custom function, use show.work=TRUE to get a table with dist values.
work <- get.taxdist(phy,fn=mean,show.work = TRUE)
# pick a row to view the calculation
row <- work %>% slice(100)
row$dist.list
row$dist
x <- row$dist.list[[1]]



# say we want to define a function
custom_fun <- function(x) {
  weights <- length(x):1
  sum( x*weights ) / sum(weights)
}
custom_fun(x)


# generate a few distance matrices:


dist.taxhorn.mean <- get.taxdist(phy,method="horn",fn=mean)

dist.taxhorn.custom <- get.taxdist(phy,method="horn",fn=custom_fun)

dist.euclidean <- distance(phy,"euclidean")

dist.bray <- distance(phy,"bray")



#view these and see which is best
dist.taxhorn.mean %>% view.hclust(phy,title="mean across levels")

dist.taxhorn.custom %>% view.hclust(phy,title="weighted mean across levels")

dist.euclidean %>% view.hclust(phy,title="euclidean")

dist.bray %>% view.hclust(phy,title="bray")

pdf("compare_distances.pdf",width=18,height=10)
dist.taxhorn.mean %>% view.hclust(phy,title="mean across levels")
dist.taxhorn.custom %>% view.hclust(phy,title="weighted mean across levels")
dist.euclidean %>% view.hclust(phy,title="euclidean")
dist.bray %>% view.hclust(phy,title="bray")
dev.off()
shell.exec("compare_distances.pdf")


#compare by group:

view_violin <- function(dist,title) {
  pairdata <- dist %>% get.pairwise() %>% 
    inner_join(pairs,by=c("sample1","sample2"))
  ggplot(pairdata,aes(x=status,y=dist,fill=status)) + geom_violin() +
    geom_boxplot(alpha=0.1) + ggtitle(title)
}


g1 <- view_violin(dist.taxhorn.mean,"mean across levels")
g2 <- view_violin(dist.taxhorn.custom,"weighted mean across levels")
g3 <- view_violin(dist.euclidean,"euclidean")
g4 <- view_violin(dist.bray,"bray")


pdf("violin_compare_groups.pdf",width=20,height=12)
grid.arrange(g1,g2,g3,g4)
dev.off()

shell.exec("violin_compare_groups.pdf")

