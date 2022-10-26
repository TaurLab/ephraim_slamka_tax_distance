Slamka taxdist project
================
myself

## title

Quarto enables you to weave together content and executable code into a
finished document. To learn more about Quarto see <https://quarto.org>.

## Running Code

When you click the **Render** button a document will be generated that
includes both content and the output of embedded code. You can embed
code like this:

``` r
1 + 1
```

    [1] 2

You can add options to executable code like this

    [1] 4

The `echo: false` option disables the printing of code (only output is
displayed).

``` r
  library(tidyverse)
  library(yingtools2)

  # given a distance matrix, return a table of pairwise       distances.
  # e.g.: dist(mtcars) %>% get.pairwise()
  get.pairwise <- function(dist,diag=TRUE) {
    mat <- as.matrix(dist,diag=diag)
    # xy <- t(combn(colnames(mat), 2))
    # data.frame(xy, dist=mat[xy]) %>% as_tibble() %>%   rename(sample1=X1,sample2=X2)
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
```

``` r
d<-dist(mtcars)
view.dist(d)
```

![](README_files/figure-gfm/unnamed-chunk-4-1.png)
