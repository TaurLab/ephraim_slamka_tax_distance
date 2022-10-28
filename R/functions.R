
library(tidyverse)
library(yingtools2) 

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
  g.hclust
  
  gg.stack(g.hclust,g.groups,g.tax,heights=c(3,1,4),align.xlim=FALSE,as.gtable = TRUE)
}
