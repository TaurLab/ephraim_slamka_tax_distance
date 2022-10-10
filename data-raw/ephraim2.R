
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

pts <- tibble(MRN=c("00036481", "00034315", "00368257"),
              BMT=as.Date(c("2018-05-09", "2018-05-31", "2017-03-17")))


pt1.bmt %>% group_by(MRN,BMT) %>%
  filter(n()>1) %>% dt()

p <- tibble(MRN="00034315",BMT=as.Date("2018-05-31"))
ssub <- phy1 %>%
  get.samp() %>%
  inner_join(p,by=c("MRN","BMT"))

physub <- prune_samples(ssub$sample,phy1) %>%
  phy.collapse()
otusub <- physub %>%
  get.otu.melt() %>%
  tax.plot(data=TRUE)
pal <- get.yt.palette2(otusub)

g.tax <- ggplot() +
  geom_col(data=otusub,aes(x=bmtday,y=pctseqs,fill=Species),width=0.95) +
  scale_fill_manual(values=pal) +
  geom_text(data=otusub,aes(x=bmtday,y=y.text,label=tax.label),angle=-90) +
  geom_text(data=ssub,aes(x=bmtday,y=0,label=bmtday),vjust=1) +
  geom_col(data=ssub,aes(x=bmtday,y=1),width=0.95,color="black",fill=NA) +
  scale_y_continuous("Relative Abundance",label=percent) +
  scale_x_continuous("Transplant Day") +
  coord_cartesian(xlim=c(-7,15)) +
  theme(legend.position="none")
g.tax

ord <- ordinate(physub,distance="bray")
g1 <- plot_ordination(physub,ord,label = "bmtday") + theme(aspect.ratio=1)


grid.arrange(g.tax,g1)
ordinate(physeq, method = "DCA", distance = "bray", formula = NULL,
         ...)


dist <- distance(physub,method="bray")
pca <- prcomp(dist)

pcadata <- pca$x[,1:2] %>% data.frame() %>%
  rownames_to_column("sample") %>%
  left_join(ssub,by="sample")
loadings <- summary(pca)$importance["Proportion of Variance",][1:2]
axis.labels <- paste0(names(loadings)," (",percent(loadings),")")

g.pca1 <- ggplot(pcadata,aes(x=PC1,y=PC2,label=bmtday)) +
  geom_point() +
  geom_text(vjust=1.1) +
  scale_color_manual(values=pal) +
  xlab(axis.labels[1]) + ylab(axis.labels[2]) +
  theme(aspect.ratio=1,legend.position="none")

pca$rotation
results$rotation


g.pca1
results <- prcomp(USArrests, scale = TRUE)


grid.arrange(g.tax,g.pca1)


mat <- as.matrix(dist)
xy <- colnames(mat) %>% combn(2) %>% t()
pwdist <- data.frame(xy, dist=mat[xy]) %>%
  rename(sample1=X1,sample2=X2,dist=dist)

adjpairs <- ssub %>%
  arrange(bmtday) %>%
  transmute(sample1=sample,
            bmtday1=bmtday,
            sample2=lead(sample),
            bmtday2=lead(bmtday)) %>%
  filter(!is.na(sample2)) %>%
  mutate(bmtday.mid=midpoint(bmtday1,bmtday2),
         dist=map2_dbl(sample1,sample2,~mat[.x,.y]))

g.tax +
  geom_text(data=adjpairs,aes(x=bmtday.mid,y=1.1,label=pretty_number(dist))) +
  geom_segment(data=adjpairs,aes(x=bmtday1+0.05,xend=bmtday2-0.05,y=1.07,yend=1.07))






library(yingtools2)
library(tidyverse)
library(ytrecipes)
library(phyloseq)
library(scales)
library(grid)
library(gridExtra)

p <- pts %>% slice(1) 

plot.pt <- function(p,method="bray") {
  ssub <- phy1 %>%
    get.samp() %>%
    inner_join(p,by=c("MRN","BMT")) %>%
    mutate(bmtday=rank(bmtday))

  physub <- prune_samples(ssub$sample,phy1) %>%
    phy.collapse()
  sample_data(physub) <- ssub %>% set.samp()
  otusub <- physub %>%
    get.otu.melt() %>%
    tax.plot(data=TRUE)
  pal <- get.yt.palette2(otusub)

  dist <- distance(physub,method=method)
  # pca <- prcomp(dist)
  # pcadata <- pca$x[,1:2] %>% data.frame() %>%
  #   rownames_to_column("sample") %>%
  #   left_join(ssub,by="sample")
  # loadings <- summary(pca)$importance["Proportion of Variance",][1:2]
  # axis.labels <- paste0(names(loadings)," (",percent(loadings),")")
  mat <- as.matrix(dist)
  adjpairs <- ssub %>%
    arrange(bmtday) %>%
    transmute(sample1=sample,
              bmtday1=bmtday,
              sample2=lead(sample),
              bmtday2=lead(bmtday)) %>%
    filter(!is.na(sample2)) %>%
    mutate(bmtday.mid=midpoint(bmtday1,bmtday2),
           dist=map2_dbl(sample1,sample2,~mat[.x,.y]))
  
  
  
  g.tax <- ggplot() +
    geom_col(data=otusub,aes(x=bmtday,y=pctseqs,fill=Species),width=0.95) +
    scale_fill_manual(values=pal) +
    geom_text(data=otusub,aes(x=bmtday,y=y.text,label=tax.label),angle=-90) +
    geom_text(data=ssub,aes(x=bmtday,y=0,label=bmtday),vjust=1) +
    geom_col(data=ssub,aes(x=bmtday,y=1),width=0.95,color="black",fill=NA) +
    scale_y_continuous("Relative Abundance",label=percent) +
    scale_x_continuous("Transplant Day") +
    theme(legend.position="none") +
    geom_point(data=adjpairs,aes(x=bmtday.mid,y=1.1,size=dist),color="steelblue") +
    geom_text(data=adjpairs,aes(x=bmtday.mid,y=1.1,label=short_number(dist))) +
    geom_segment(data=adjpairs,aes(x=bmtday1+0.05,xend=bmtday2-0.05,y=1.07,yend=1.07)) +
    expand_limits(size=0) +
    scale_size_continuous(range=c(1,10)) +
    ggtitle(method)
  g.tax

}

plot.pt.multi <- function(p,methods=c("bray","euclidean")) {
  glist <- methods %>% map(~{
    message(.x)
    p %>% plot.pt(method=.x)
  })
  do.call(grid.arrange,glist)
}

load("tyler.phylo.RData")



p <- tibble(MRN="00034315",BMT=as.Date("2018-05-31"))
# plot.pt(p,"bray")
methods <- c("bray","euclidean","jsd","manhattan","canberra")
pt1 %>% filter(!is.na(BMT)) %>% slice(7) %>% plot.pt.multi(methods)

pt1 %>% filter(!is.na(BMT)) %>% slice(7) %>% select(MRN,BMT) %>% copy.as.Rcode()



distance
p <- tibble(MRN="00034315",BMT=as.Date("2018-05-31"))
methods <- c("bray","euclidean","morisita","horn","jsd","manhattan","canberra")


p %>% plot.pt.multi(methods)


pdf("three.pdf",width=20,height=15)
pts %>% group_by(1:n()) %>%
  group_map(~plot.pt.multi(.x,methods))
dev.off()
shell.exec("three.pdf")



p <- tyler.phy

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
















