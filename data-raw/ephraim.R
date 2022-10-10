
library(yingtools2)
library(ytrecipes)
library(tidyverse)
library(phyloseq)
library(ape)
library(ggtree)
library(scales)
library(grid)
library(gridExtra)
library(ggbiplot)
rm(list=ls())



# 
# pca <- prcomp(mtcars)
# pcadata <- pca$x[,1:2] %>% data.frame() %>% 
#   rownames_to_column("car")
# loadings <- summary(pca)$importance["Proportion of Variance",][1:2]
# axis.labels <- paste0(names(loadings)," (",percent(loadings),")")
# g1 <- ggplot(pcadata,aes(x=PC1,y=PC2,label=car)) + 
#   geom_point() + geom_text() +
#   scale_color_manual(values=pal) +
#   xlab(axis.labels[1]) + ylab(axis.labels[2]) + 
#   theme(aspect.ratio=1,legend.position="none")
# 
# 
# 
# dist <- dist(mtcars)
# pca2 <- prcomp(dist)
# pcadata2 <- pca2$x[,1:2] %>% data.frame() %>% 
#   rownames_to_column("car")
# loadings2 <- summary(pca2)$importance["Proportion of Variance",][1:2]
# axis.labels2 <- paste0(names(loadings2)," (",percent(loadings2),")")
# g2 <- ggplot(pcadata2,aes(x=PC1,y=PC2,label=car)) + 
#   geom_point() + geom_text() +
#   scale_color_manual(values=pal) +
#   xlab(axis.labels2[1]) + ylab(axis.labels2[2]) + 
#   theme(aspect.ratio=1,legend.position="none")
# 
# grid.arrange(g1,g2,nrow=1)
# 
# 
# ordination






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
g1 <- plot_ordination(physub,ord) + theme(aspect.ratio=1)
g.pca1
grid.arrange(g1,g.pca1)


pcaCars <- princomp(mtcars, cor = TRUE)
carsDf <- data.frame(pcaCars$scores, "cluster" = factor(carsClusters))

library(ggplot2)
library(ggrepel)
ggplot(carsDf,aes(x=Comp.1, y=Comp.2)) +
  geom_text_repel(aes(label = rownames(carsDf))) +
  theme_classic() +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  geom_point(aes(color = cluster), alpha = 0.55, size = 3) +
  xlab("PC1") +
  ylab("PC2") + 
  xlim(-5, 6) + 
  ggtitle("PCA plot of Cars")




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

plot.pt <- function(p,method="bray") {
  ssub <- phy1 %>% 
    get.samp() %>% 
    inner_join(p,by=c("MRN","BMT"))
  physub <- prune_samples(ssub$sample,phy1) %>% 
    phy.collapse()
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




p <- tibble(MRN="00034315",BMT=as.Date("2018-05-31"))
# plot.pt(p,"bray")
methods <- c("bray","euclidean","jsd","manhattan","canberra")
pt1 %>% filter(!is.na(BMT)) %>% slice(4) %>% plot.pt.multi(methods)







