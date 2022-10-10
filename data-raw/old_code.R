
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

