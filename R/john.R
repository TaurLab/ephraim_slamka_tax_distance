
##### basically the code of analyze, without output

# load data and functions.R (run this first) ---------------------------------------------------------------

library(phyloseq) #phyloseq objects house 16S data
library(grid) #for manipulating plots
library(gridExtra) #for manipulating plots
library(vegan) #calculates known distance metrics like bray-curtis
library(ggtree) #for drawing trees
library(ape) #used to manipulate trees
library(tidyverse)
library(yingtools2) #ying's suite of data tools
library(patchwork)
library(rlang)
setwd(here::here())
rm(list=ls())
load("data/mock_phylo.RData")
# phy <- readRDS("data/mock_phylo.rds")
# source("R/functions.R")
s <- get.samp(mock.phy)

# assess all pairwise combinations and classify as one of 4 groups.
s1 <- get.samp(mock.phy) %>% select(sample,pt,pt.day,pt.day.samp) %>% rename_with(~paste0(.,"1"))
s2 <- get.samp(mock.phy) %>% select(sample,pt,pt.day,pt.day.samp) %>% rename_with(~paste0(.,"2"))
pairs <- sample_names(mock.phy) %>%
  combn(2) %>% t() %>%
  as_tibble() %>%
  rename(sample1=V1,sample2=V2) %>%
  left_join(s1,by="sample1") %>%
  left_join(s2,by="sample2") %>%
  mutate(status=case_when(
    pt.day.samp1==pt.day.samp2 ~ "same pt/day/sample",
    pt.day1==pt.day2 ~ "same pt/day (diff sample)",
    pt1==pt2 ~ "same pt (diff day/sample)",
    TRUE ~ "diff pt",
  ))
rm(s1,s2)


prs <- pairs %>% 
  group_by(status) %>%
  sample_n(5) %>%
  ungroup()



dists <- c("bray","euclidean","chao","horn","jsd","taxhorn")



data <- dists %>% map(~{
  prs %>% mutate(dname=.x,
                 dist=calc.pairwise(sample1,sample2,mock.phy,method=.x))
}) %>% list_rbind()

test <- data %>% 
  group_by(dname) %>%
  summarize(test1.median.dist=median(dist[status=="same pt/day/sample"]),
            test1.max.dist=max(dist[status=="same pt/day/sample"]),
            test2.max.dist=max(dist),
            .groups="drop")

test


gg.list <- dists %>% map(~{
  prs %>% mutate(dist=calc.pairwise(sample1,sample2,mock.phy,method=.x)) %>%
    ggplot(aes(x=status,y=dist,fill=status)) + geom_violin() + ggtitle(.x)
})

inject(wrap_plots(!!!gg.list,guides="collect"))






# create various distances, including the customized ------------------------------------------------





# # usual distances 
# 
# dist.manhattan <- distance(phy,"manhattan")
# dist.bray <- distance(phy,"bray")
# dist.euclidean <- distance(phy,"euclidean")
# dist.horn <- distance(phy,"horn")
# dist.unifrac <- distance(phy,"unifrac")
# dist.wunifrac <- distance(phy,"wunifrac")
# 
# 
# 
# 
# #ignore what follows, just an attempt at learning, any mention of these functions was just an attempt at learning-----------
# #this isn't finished, I'm just trying to figure out how all of this works.
# athing<- function(list){sqrt(sum(list)) }
# aathing<- function(list) {n<-length(list)
# pwr<- 1/n
# (sum(list))^pwr}
# 
# aabthing<-function(x) {n<-length(x)
# pwr<- 1/n
# this<-((sum(x))^pwr)
# this/n}
# 
# aacthing<-function(x) {n<-length(x)
# pwr<- 1/n
# (sum(x))^pwr}
# geomean<- function(set){n<-length(set)
# pwr<-1/n
# (prod(set))^pwr}
# 
# #working here
# abthing<-function(x){
#   n<-(sum(x))*(prod(x))
#   m<-(n^(1/length(x)))/x
#   mean(m)/length(m)
#   #pwr<- 1/(length(n))
#   # (prod(n))^pwr
# }
# acthing<-function(x){
#   n<-(sum(x))*(prod(x))
#   m<-(n^(1/length(x)))/x
#   mean(m)/mean(x)
#   #pwr<- 1/(length(n))
#   # (prod(n))^pwr
# }
# #test these two ^^vv
# geomean2<- function(x){
#   
#   n<-length(x)
#   pwr<- (1/n)
#   (prod(x))^pwr}
# 
# bthing<- function(x){
#   n<-x^x
#   m<-sum(n)
#   m/length(n)
# }
# bathing<-function(x){
#   n<-x^(1/(x^2))
#   m<-sum(n)
#   m/length(n)
# }
# bbthing<-function(x){
#   n<-x^(sin(x))
#   m<-sum(n)
#   m/length(n)
# }
# bcthing<-function(x){
#   n<-x^tan(x)
#   m<-sum(n)
#   m/length(n)
# }
# #this is the area I'm working on
# #note for future ephraim, if elses are doomed to fail in 
# #nonspecific cases. try using a median or mean or something, idk.
# #be clever future E. past E believes in you.
# bdthing<-function(x){
#   
#   d<-mean(x)
#   n<-(x)^x
#   o<-(x)^(1/x)
#   w<-mean(n)
#   s<-mean(o)
#   ifelse(x>d, return(w), return(s))
# }
# 
# nsk<-function(x){
#   mean(x[-1])
# }
# wnsk1<-function(x){
#   set<-x[-1]
#   weights<-length(set):1
#   sum(set*weights)/sum(weights)
# }
# wnsk2<-function(x){
#   set<-x[-1]
#   #weights<-length(set):1
#   weights<-c(2,3,4,5,4,2,1)
#   sum(set*weights)/sum(weights)
# }
# weighted.mean <- function(x) {
#   weights <- length(x):1
#   sum( x*weights ) / sum(weights)
# }
# 
# 
# 
# 
# 
# #custom distances
# dist.taxhorn.mean <- get.taxdist(phy,fn=mean,method="horn")
# dist.taxhorn.weightedmean <- get.taxdist(phy,fn=weighted.mean,method="horn")
# dist.taxhorn.aabthing<- get.taxdist(phy,fn=aabthing,method="horn")
# dist.taxhorn.abthing<- get.taxdist(phy,fn=abthing, method="horn")
# dist.taxhorn.acthing<-get.taxdist(phy,fn=acthing,method="horn")
# dist.taxhorn.bthing<-get.taxdist(phy,fn=bthing,method="horn")
# dist.taxhorn.bdthing<-get.taxdist(phy,fn=bdthing,method="horn")
# dist.taxhorn.bbthing<-get.taxdist(phy,fn=bbthing,method="horn")
# dist.taxhorn.nsk<-get.taxdist(phy,fn=nsk,method="horn")
# dist.taxhorn.wnsk2<-get.taxdist(phy,fn=wnsk2,method="horn")
# 
# dist.taxhorn.wnsk1<-get.taxdist(phy,fn=wnsk1,method="horn")
# dist.taxhorn.bcthing<-get.taxdist(phy,fn=bcthing,method="horn")
# # dist.taxhorn.scratch <- get.taxdist(phy,fn=mean,method="horn",show.work = TRUE)
# # x <- dist.taxhorn.scratch$dist.list[[10]]
# # mean(x)
# # weighted.mean(x)
# # g.hc.aabthing<-view.hclust(dist.taxhorn.aabthing, title= "aabthing")
# # grid.draw(g.hc.aabthing)
# 
# 
# 
# # examine each pt's samples at different tax levels -------------------------------------------------
# 
# # get a subset phyloseq and run plot.dist. 
# physub <- subset_samples(phy,pt=="HV")
# g.tyler.bray <- plot.dist(physub,"bray")
# g.tyler.horn <- plot.dist(physub,"horn")
# 
# # to plot this for all subjects, create list of ggplots:
# glist <- get.samp(phy) %>% group_by(pt) %>%
#   group_map(~{
#     pt <- .x$pt[1]
#     metric <- "horn"
#     physub <- prune_samples(.x$sample,phy) %>% prune_unused_taxa()
#     plot.dist(physub, method=metric, sortby=c("pt","pt.day","pt.day.samp")) +
#       ggtitle(str_glue("subject: {pt} / metric: {metric}"))
#   },.keep=TRUE)
# 
# glist[[1]]
# 
# # pdf("plots/mock_metrics.pdf",width=34,height=17)
# # glist
# # dev.off()
# # shell.exec(normalizePath("plots/mock_metrics.pdf"))
# 
# 
# # view some of the pairs to get an idea ----------------------------------------------------------
# 
# 
# #select 10 pairs from each comparison and show
# otu <- phy %>% get.otu.melt() %>%
#   tax.plot(data=TRUE)
# pal <- get.yt.palette2(otu)
# 
# pairs.select <- pairs %>% 
#   group_by(status) %>% slice(1:10) %>% 
#   transmute(pair.id=row_number(),status,sample1,sample2) %>% 
#   ungroup() %>% 
#   pivot_longer(cols=c(sample1,sample2),values_to="sample") %>%
#   left_join(otu,by="sample")
# 
# g.pairs.select <- 
#   ggplot(pairs.select,aes(x=name,y=pctseqs,fill=Species)) + geom_col(show.legend=FALSE) +
#   # geom_text(aes(y=y.text,label=tax.label),angle=-90) +
#   scale_fill_manual(values=pal) +
#   facet_grid(status~pair.id,scales="free_x",space="free_x")
# # g.pairs.select
# 
# # pdf("plots/pairs.grouping.taxview.pdf",width=20,height=12)
# # g.pairs.select
# # dev.off()
# # 
# # shell.exec(normalizePath("plots/pairs.grouping.taxview.pdf"))
# 
# # view all same day samples (not same sample)-----------------------------------------------
# 
# s.sameday <- s %>% 
#   group_by(pt.day.samp) %>%
#   slice(1) %>%
#   group_by(pt.day) %>%
#   filter(n()>1) %>%
#   ungroup()
# 
# otu.sameday <- phy %>% prune_samples(s.sameday$sample,.) %>%
#   get.otu.melt() %>%
#   tax.plot(data=TRUE)
# pal <- get.yt.palette2(otu.sameday)
# 
# g.sameday.sample <- ggplot(otu.sameday,aes(x=sample,y=pctseqs,fill=Species)) +
#   geom_col(show.legend = FALSE) + 
#   geom_text(aes(y=y.text,label=tax.label),angle=-90) +
#   scale_fill_manual(values=pal) +
#   facet_grid(.~pt.day,scales="free_x",space="free_x")
# # g.sameday.sample
# # ggsave("plots/g.same.sample.pdf",g.sameday.sample,width=15,height=8)
# # shell.exec("plots/g.same.sample.pdf")
# 
# 
# # view all resequenced samples ----------------------------------------------------------
# 
# s.reseq <- s %>% group_by(pt.day.samp) %>%
#   filter(n()>1)
# 
# otu.reseq <- phy %>% prune_samples(s.reseq$sample,.) %>%
#   get.otu.melt() %>%
#   tax.plot(data=TRUE)
# pal <- get.yt.palette2(otu.sameday)
# 
# g.reseq.sample <- ggplot(otu.reseq,aes(x=sample,y=pctseqs,fill=Species)) +
#   geom_col(show.legend = FALSE) + 
#   geom_text(aes(y=y.text,label=tax.label),angle=-90) +
#   scale_fill_manual(values=pal) +
#   facet_grid(.~pt.day.samp,scales="free_x",space="free_x")
# 
# # g.reseq.sample
# # 
# # ggsave("plots/g.all.reseq.pdf",g.reseq.sample,width=20,height=10)
# # shell.exec(normalizePath("plots/g.all.reseq.pdf"))
# 
# 
# 
# # violin of distances -----------------------------------------------------
# 
# #compare by group:
# 
# pairs <- pairs %>% 
#   mutate(status2=recode2(status,
#                          recodes=c("diff pt"="(1) Different\nsubjects", 
#                                    "same pt (diff day/sample)"="(2) Same subj,\ndiff days/samps",
#                                    "same pt/day (diff sample)"="(3) Same subj/day,\ndiff samps", 
#                                    "same pt/day/sample"="(4) Re-seq of\nsame samp"),
#                          as.factor=TRUE))
# 
# 
# view_violin <- function(dist,title) {
#   pairdata <- dist %>% get.pairwise() %>%
#     inner_join(pairs,by=c("sample1","sample2"))
#   ggplot(pairdata,aes(x=status2,y=dist,fill=status2)) + geom_violin() +
#     geom_boxplot(alpha=0.2) + ggtitle(title) + expand_limits(y=c(0,1)) + theme(legend.position="none") +
#     scale_y_continuous(str_glue("{title}\ndistance")) +
#     scale_x_discrete("")
# }
# 
# 
# 
# g.v.manhattan <- view_violin(dist.manhattan,title="Manhattan")
# g.v.bray <- view_violin(dist.bray,title="Bray-Curtis")
# g.v.euclidean <- view_violin(dist.euclidean,title="Euclidean") + scale_y_continuous("Euclidean\ndistance",labels = pretty_number)
# g.v.horn <- view_violin(dist.horn,title="Horn-Morisita")
# g.v.unifrac <- view_violin(dist.unifrac,title="Unifrac (unweighted)")
# g.v.wunifrac <- view_violin(dist.wunifrac,title="wunifrac")
# g.v.taxhorn.mean <- view_violin(dist.taxhorn.mean,title="taxhorn.mean")
# g.v.taxhorn.weightedmean <- view_violin(dist.taxhorn.weightedmean,title="taxhorn.weightedmean")
# g.v.taxhorn.aabthing <- view_violin(dist.taxhorn.aabthing, title="taxhorn.aabthing")
# g.v.taxhorn.abthing<-view_violin(dist.taxhorn.abthing,title="taxhorn.abthing")
# g.v.taxhorn.acthing<-view_violin(dist.taxhorn.acthing,title="taxhorn.acthing")
# g.v.taxhorn.bthing<-view_violin(dist.taxhorn.bthing,title="taxhorn.bthing")
# g.v.taxhorn.bdthing<-view_violin(dist.taxhorn.bdthing,title="taxhorn.bdthing")
# g.v.taxhorn.bbthing<-view_violin(dist.taxhorn.bbthing,title="taxhorn.bbthing")
# g.v.taxhorn.wnsk2<-view_violin(dist.taxhorn.wnsk2,title="taxhorn.wnsk2")
# 
# g.v.taxhorn.wnsk1<-view_violin(dist.taxhorn.wnsk1,title="taxHorn")
# g.v.taxhorn.nsk<-view_violin(dist.taxhorn.nsk,title="taxhorn.nsk")
# g.v.taxhorn.bcthing<-view_violin(dist.taxhorn.bcthing,title="taxhorn.bcthing")
# vlist <- list( 
#   #g.v.bray,
#   #g.v.manhattan,
#   #g.v.euclidean,
#   g.v.horn,
#   g.v.unifrac,
#   #g.v.wunifrac,
#   g.v.taxhorn.mean,
#   #g.v.taxhorn.weightedmean,
#   #g.v.taxhorn.aabthing,
#   #g.v.taxhorn.abthing,
#   #g.v.taxhorn.acthing,
#   #g.v.taxhorn.bthing,
#   #g.v.taxhorn.bdthing,
#   #g.v.taxhorn.bbthing,
#   #g.v.taxhorn.bcthing),
#   g.v.taxhorn.nsk,
#   g.v.taxhorn.wnsk1,
#   g.v.taxhorn.wnsk2)
# 
# 
# 
# 
# 
# # do.call(grid.arrange,vlist)
# 
# # pdf("violin_compare.pdf",width=6,height=16)
# # gg.stack(
# #   g.v.bray,
# #   g.v.euclidean,
# #   g.v.horn,
# #   g.v.unifrac,
# #   g.v.taxhorn.wnsk1,
# #   adjust.themes = FALSE,
# #   newpage = FALSE
# # )
# # dev.off()
# # shell.exec("violin_compare.pdf")
# 
# # mg <- marrangeGrob(vlist,ncol=3,nrow=3)
# # ggsave("plots/violin_compare_groups.pdf",mg,width=20,height=12)
# # shell.exec(normalizePath("plots/violin_compare_groups.pdf"))
# 
# # pca of samples ---------------------------------------------------------------------
# 
# view.pca <- function(metric) {
#   # metric <- "euclidean"
#   pcadata <- ordinate(phy, method = "NMDS", distance = metric) %>%
#     scores(display = "sites") %>%
#     as.data.frame() %>%
#     rownames_to_column("sample") %>%
#     left_join(s, by = "sample")
#   s.dups <- pcadata %>% filter(has.sameday)
#   g.pca <- ggplot(pcadata, aes(x = NMDS1, y = NMDS2)) +
#     geom_path(data=s.dups,aes(group=pt.day),color="dark gray",size=1) +
#     geom_point(aes(color=pt),alpha = 0.75) +
#     # geom_text(data=s.dups,aes(label=pt.day)) +
#     theme(aspect.ratio = 1) + ggtitle(metric)
#   g.pca
# }
# 
# 
# # view.pca("bray")
# # view.pca("horn")
# # 
# # metrics <- c("euclidean","bray","unifrac","wunifrac","horn")
# # g.pca.list <- metrics %>% map(view.pca)
# # do.call(grid.arrange,c(g.pca.list,list(nrow=2)))
# 
# 
# # hclust of samples -------------------------------------------------------
# 
# 
# physub <- phy %>% subset_samples(pt!="PT3" &  day<100)
# 
# g.hc.wnsk1<-view.hclust(dist.taxhorn.wnsk1,.phy=physub,title="Taxhorn")
# g.hc.manhattan <- view.hclust(dist.manhattan,.phy=physub,title="manhattan")
# g.hc.bray <- view.hclust(dist.bray,.phy=physub,title="bray")
# g.hc.euclidean <- view.hclust(dist.euclidean,.phy=physub,title="Euclidean")
# g.hc.horn <- view.hclust(dist.horn,.phy=physub,title="Horn-Morisita")
# g.hc.unifrac <- view.hclust(dist.unifrac,.phy=physub,title="unifrac")
# g.hc.wunifrac <- view.hclust(dist.wunifrac,.phy=physub,title="wunifrac")
# g.hc.taxhorn.mean <- view.hclust(dist.taxhorn.mean,.phy=physub,title="taxhorn.mean")
# g.hc.taxhorn.weightedmean <- view.hclust(dist.taxhorn.weightedmean,.phy=physub,title="taxhorn.weightedmean")
# 
# 
# 
# 
# 
# # pdf("compare.pdf",height=16,width=17)
# # grid.arrange(g.hc.euclidean,
# #              g.hc.horn,
# #              g.hc.wnsk1,
# #              g.hc.taxhorn.mean,
# #              ncol=1)
# # dev.off()
# # 
# # shell.exec("compare.pdf")
# 
# 
# 
# 
# # grid.newpage()
# # grid.draw(g.hc.manhattan)
# 
# # g.hc.taxhorn.weightedmean <- view.hclust(dist.taxhorn.weightedmean,title="taxhorn.weightedmean",label.pct.cutoff = 0.1)
# # grid.draw(g.hc.taxhorn.weightedmean)
# 
# # mg <- marrangeGrob(list(
# #   g.hc.wnsk,
# #   g.hc.bray,
# #   g.hc.manhattan,
# #   g.hc.euclidean,
# #   g.hc.horn,
# #   g.hc.aabthing,
# #   g.hc.unifrac,
# #   g.hc.wunifrac,
# #   g.hc.taxhorn.mean,
# #   g.hc.taxhorn.weightedmean),
# #   ncol=1,nrow=1)
# 
# 
# 
# # ggsave("plots/hclust_compare.pdf",mg,width=20,height=12)
# # shell.exec(normalizePath("plots/hclust_compare.pdf"))
# 
# 
# 
# 
# 
# # view stacked samples (slide) --------------------------------------------
# 
# 
# 
# s <- get.samp(phy) %>%
#   group_by(pt) %>% 
#   mutate(pos=row_number(day),
#          status=coalesce_indicators(has.reseq,has.sameday,else.value="none",first.hit.only = TRUE),
#          pt.label=recode2(pt,c("HV"="healthy volunteer", "PT1"="Patient 1", "PT2"="Patient 2", "PT3"="Patient 3"))) %>%
#   ungroup()
# 
# s.dups <- s %>% 
#   filter(status!="none")
# 
# 
# otu <- phy %>% phy.collapse() %>% 
#   get.otu.melt() %>% 
#   tax.plot(data=TRUE,label.pct.cutoff = 0.5) %>%
#   mutate(tax.label=sub(" ","\n",tax.label)) %>% 
#   left_join_replace(select(s,sample,pos,status,pt.label),by="sample")
# pal <- get.tax.palette(otu)
# legend <- get.tax.legend(fontsize=4) %>%  
#   gtable::gtable_add_rows(unit(9,"null")) %>%  
#   gtable::gtable_add_rows(unit(9,"null"),pos=0)
# 
# 
# g.tax <- ggplot() +
#   expand_limits(y=-0.05) +
#   geom_col(data=otu,aes(x=pos,y=pctseqs,fill=Species),show.legend = FALSE) +
#   geom_text(data=otu,aes(x=pos,y=y.text,label=tax.label),size=2.5,angle=-90,lineheight=0.75) +
#   geom_text(data=s,aes(x=pos,y=0,label=day),vjust=1.08,size=2.5) +
#   scale_fill_manual(values=pal)  +
#   scale_x_continuous("Time (day)") +
#   scale_y_continuous("Bacterial Relative Abundance",labels = scales::percent,breaks=c(0,0.5,1)) +
#   facet_grid(pt.label ~ .,) +
#   theme(legend.position = c(0.85,0.93),
#         axis.text.x=element_blank(),
#         axis.ticks = element_blank())
# 
# g.dups <- list(
#   geom_col(data=s.dups,aes(x=pos,y=1,color=status),fill=NA,size=1.5,show.legend=FALSE),
#   geom_path(data=s,aes(x=pos,y=1.05,group=pt.day,color=status),size=1.5),
#   scale_color_manual("Similar Samples",values=c("has.reseq"="blue", "has.sameday"="red"),
#                      labels=c("has.reseq"="Same sample, resequenced", "has.sameday"="Distinct samples on same day"))
# )
# 
# 
# png("plots/taxdata.png",height=7,width=14,units="in",res=150)
# grid.arrange(
#   g.tax,
#   legend,
#   nrow=1,
#   widths=c(3,1))
# dev.off()
# 
# shell.exec("plots/taxdata.png")
# 
# 
# 
# png("plots/taxdata2.png",height=8,width=16,units="in",res=150)
# grid.arrange(
#   g.tax+g.dups,
#   legend,
#   nrow=1,
#   widths=c(3,1))
# dev.off()
# 
# # shell.exec("plots/taxdata2.png")
# 
# 
# 
# # adjacent dists ----------------------------------------------------------
# 
# 
# 
# physub <- phy %>% subset_samples(pt %in% c("HV","PT1"))
# 
# s <- get.samp(physub) %>%
#   group_by(pt) %>% 
#   mutate(pos=row_number(day),
#          status=coalesce_indicators(has.reseq,has.sameday,else.value="none",first.hit.only = TRUE),
#          pt.label=recode2(pt,c("HV"="healthy volunteer", "PT1"="Patient 1", "PT2"="Patient 2", "PT3"="Patient 3"))) %>%
#   ungroup()
# 
# s.dups <- s %>% 
#   filter(status!="none")
# 
# 
# otu <- physub %>% phy.collapse() %>% 
#   get.otu.melt() %>% 
#   tax.plot(data=TRUE,label.pct.cutoff = 0.5) %>%
#   mutate(tax.label=sub(" ","\n",tax.label)) %>% 
#   left_join_replace(select(s,sample,pos,status,pt.label),by="sample")
# pal <- get.tax.palette(otu)
# legend <- get.tax.legend(fontsize=4) %>%  
#   gtable::gtable_add_rows(unit(9,"null")) %>%  
#   gtable::gtable_add_rows(unit(9,"null"),pos=0)
# 
# 
# g.tax <- ggplot() +
#   expand_limits(y=-0.05) +
#   geom_col(data=otu,aes(x=pos,y=pctseqs,fill=Species),show.legend = FALSE) +
#   geom_text(data=otu,aes(x=pos,y=y.text,label=tax.label),size=2.5,angle=-90,lineheight=0.75) +
#   geom_text(data=s,aes(x=pos,y=0,label=day),vjust=1.08,size=2.5) +
#   scale_fill_manual(values=pal)  +
#   scale_x_continuous("Time (day)") +
#   scale_y_continuous("Bacterial Relative Abundance",labels = scales::percent,breaks=c(0,0.5,1)) +
#   facet_grid(pt.label ~ .,) +
#   theme(legend.position = c(0.85,0.93),
#         axis.text.x=element_blank(),
#         axis.ticks = element_blank())
# 
# g.dups <- list(
#   geom_col(data=s.dups,aes(x=pos,y=1,color=status),fill=NA,size=1.5,show.legend=FALSE),
#   geom_path(data=s,aes(x=pos,y=1.05,group=pt.day,color=status),size=1.5),
#   scale_color_manual("Similar samples",values=c("has.reseq"="blue", "has.sameday"="red"),
#                      labels=c("has.reseq"="Same sample, resequenced", "has.sameday"="Distinct samples on same day"))
# )
# 
# 
# # adjacent only
# adj <- s %>% group_by(pt) %>%
#   arrange(pt,pos) %>%
#   mutate(sample1=sample,
#          sample2=lead(sample),
#          pos1=pos,
#          pos2=lead(pos)) %>%
#   ungroup() %>% 
#   select(pt,pt.label,sample1,sample2,pos1,pos2) %>% 
#   filter(!is.na(sample2)) %>%
#   mutate(adjacent=TRUE,
#          mid.pos=midpoint(pos1,pos2))
# 
# pairs.euclidean <- dist.euclidean %>%
#   get.pairwise() %>%
#   rename(euclidean.dist=dist)
# pairs.bray <- dist.bray %>%
#   get.pairwise() %>%
#   rename(bray.dist=dist)
# pairs.horn <- dist.horn %>%
#   get.pairwise() %>%
#   rename(horn.dist=dist)
# pairs.unifrac <- dist.unifrac %>%
#   get.pairwise() %>%
#   rename(unifrac.dist=dist)
# pairs.taxhorn <- dist.taxhorn.wnsk1 %>%
#   get.pairwise() %>%
#   rename(taxhorn.dist=dist)
# 
# pairs2  <- pairs %>%
#   inner_join(pairs.euclidean,by = c("sample1", "sample2")) %>% 
#   inner_join(pairs.bray,by = c("sample1", "sample2")) %>% 
#   inner_join(pairs.horn,by = c("sample1", "sample2")) %>% 
#   inner_join(pairs.unifrac,by = c("sample1", "sample2")) %>% 
#   inner_join(pairs.taxhorn,by = c("sample1", "sample2")) %>%
#   mutate(pos1=s$pos[match(sample1,s$sample)],
#          pos2=s$pos[match(sample2,s$sample)],
#          adjacent=(pt1==pt2) & (abs(pos1-pos2)==1))
# 
# 
# adj.pairs2 <- pairs2 %>%
#   filter(adjacent) %>%
#   mutate(mid.pos=(pos1+pos2)/2,
#          pt.label=recode2(pt1,c("HV"="healthy volunteer", "PT1"="Patient 1", "PT2"="Patient 2", "PT3"="Patient 3")))
# 
# adj.plot <- function(col=bray.dist,y=1.1,color="steelblue",lbl=NULL,size=3) {
#   col <- enquo(col)
#   if (is.null(lbl)){
#     lbl <- as_label(col)  
#   }
#   list(
#     expand_limits(size=c(0,1)),
#     geom_point(data=adj.pairs2,aes(x=mid.pos,y=y,size=!!col),color=color,show.legend = FALSE),
#     geom_text(data=adj.pairs2,aes(x=mid.pos,y=y,label=pretty_number(!!col)),size=size,check_overlap = TRUE),
#     annotate("text",x=-1,y=y,label=lbl)
#   )
# }
# 
# png("plots/adj_dist_compare.png",width=14,height=8,units = "in",res=150)
# g.tax + g.dups + 
#   expand_limits(y=1.6) +
#   adj.plot(bray.dist,y=1.1,color="steelblue",lbl="Bray-Curtis") +
#   adj.plot(unifrac.dist,y=1.3,color="darkolivegreen3",lbl="Unifrac") +
#   # adj.plot(taxhorn.dist,y=1.5,color="red",lbl="taxHorn") +
#   scale_size_continuous(range=c(0.5,10))
# dev.off()
# shell.exec("plots/adj_dist_compare.png")
# 
# 
# png("plots/adj_dist_compare2.png",width=14,height=8,units = "in",res=150)
# g.tax + g.dups + 
#   expand_limits(y=1.6) +
#   adj.plot(bray.dist,y=1.1,color="steelblue",lbl="Bray-Curtis") +
#   adj.plot(unifrac.dist,y=1.3,color="darkolivegreen3",lbl="Unifrac") +
#   adj.plot(taxhorn.dist,y=1.5,color="red",lbl="taxHorn") +
#   scale_size_continuous(range=c(0.5,10))
# dev.off()
# 
# 
# shell.exec("plots/adj_dist_compare2.png")
# 
# 
# 
# # hclust slide ------------------------------------------------------------
# 
# 
# view.hclust2 <- function(dist,.phy=phy,title="",label.pct.cutoff=0.3) {
#   # by.group <- "pt.day"
#   
#   hc <- hclust(dist)
#   tr <- as.phylo(hc)
#   
#   droptips <- setdiff(tr$tip.label,sample_names(.phy))
#   tr <- drop.tip(tr,droptips)
#   gt <- ggtree(tr) %<+% get.samp(.phy)
#   gd <- gt$data
#   pad <- (max(gd$x)-min(gd$x)) * 0.025
#   ylim <- range(gd$x) + c(0,pad)
#   s.dups.hc <- gd %>% filter(has.sameday) %>%
#     arrange(y) %>%
#     mutate(lbl=ifelse(has.sameday,pt.day,NA_character_),
#            row=get.row(y,y,row=pt.day,min.gap=1)) 
# 
#   map <- gd %>% filter(isTip) %>%
#     select(sample=label,x,y)
#   xlim <- range(map$y) + c(-0.5,0.5)
#   otu <- .phy %>% get.otu.melt() %>%
#     left_join(map,by="sample") %>%
#     tax.plot(data=TRUE,label.pct.cutoff = label.pct.cutoff) %>%
#     mutate(tax.label=sub(" ","\n",tax.label))
#   pal <- get.yt.palette2(otu)
#   
#   g.tax <- ggplot() + 
#     geom_col(data=otu,aes(x=y,y=pctseqs,fill=Species),show.legend=FALSE) +
#     geom_text(data=otu,aes(x=y,y=y.text,label=tax.label),angle=-90,lineheight=0.75) +
#     scale_fill_manual(values=pal) +
#     coord_cartesian(xlim=xlim,expand=FALSE) +
#     theme(axis.text = element_blank(),
#           axis.title.y=element_blank()) +
#     xlab("Sample")
#   # +  
#   #   geom_path(data=s.dups.hc,aes(x=y,y=1.02,group=pt.day,color=has.reseq),size=1) +
#   #   geom_point(data=s.dups.hc,aes(x=y,y=1.02,group=pt.day,color=has.reseq)) +
#   #   # geom_path(data=s.dups.hc,aes(x=y,y=-0.02,group=pt.day,color=has.reseq),size=1) +
#   #   scale_color_manual("Similar Samples",values=c("TRUE"="red","FALSE"="blue"),
#   #                      label=c("TRUE"="Same sample, resequenced","FALSE"="Distinct samples on same day"))
#   
#   row.range <- range(s.dups.hc$row) + c(-0.5,0.5)
#   g.groups <- ggplot() + 
#     geom_path(data=s.dups.hc,aes(x=y,y=row,group=pt.day,color=has.reseq)) + 
#     geom_point(data=s.dups.hc,aes(x=y,y=row,color=has.reseq)) +
#     geom_segment(data=map,aes(x=y,xend=y,y=row.range[1],yend=row.range[2]),alpha=0.35) +
#     scale_color_manual("Similar Samples",values=c("TRUE"="red","FALSE"="blue"),
#                        label=c("TRUE"="Sample sample, resequenced","FALSE"="Distinct samples on same day")) +
#     coord_cartesian(xlim=xlim,expand=FALSE) + expand_limits(y=row.range) + theme_void()
#   
#   g.hclust <- gt + 
#     expand_limits(x=ylim) +
#     geom_point2(data=gd,aes(subset=isTip,color=pt),alpha=0.75) +
#     # geom_point2(data=gd,aes(subset=isTip & has.reseq,size=has.sameday),shape=1) +
#     # geom_text2(data=gd,aes(subset=isTip & has.sameday, label=pt.day.samp, size=has.sameday),angle=-90,hjust=-0.05,size=3) +
#     scale_x_reverse() +
#     coord_flip(ylim=xlim,expand=FALSE) +
#     ggtitle(title)
#   
#   
#   gg.stack(g.hclust,g.groups,g.tax,heights=c(4,1,8),align.xlim=FALSE,as.gtable = TRUE)
#   # gg.stack(g.hclust,g.tax,heights=c(1.5,4),align.xlim=FALSE,as.gtable = TRUE)
# }
# 
# 
# 
# physub <- phy %>% subset_samples(pt %in% c("HV","PT1"))
# 
# g.hc.wnsk1 <- view.hclust2(dist.taxhorn.wnsk1,.phy=physub,title="taxHorn")
# g.hc.bray <- view.hclust2(dist.bray,.phy=physub,title="Bray-Curtis")
# g.hc.manhattan <- view.hclust2(dist.manhattan,.phy=physub,title="manhattan")
# g.hc.euclidean <- view.hclust2(dist.euclidean,.phy=physub,title="Euclidean")
# g.hc.horn <- view.hclust2(dist.horn,.phy=physub,title="Horn-Morisita")
# g.hc.unifrac <- view.hclust2(dist.unifrac,.phy=physub,title="unifrac")
# g.hc.wunifrac <- view.hclust2(dist.wunifrac,.phy=physub,title="wunifrac")
# 
# grid.arrange(
#   # g.hc.bray,
#   g.hc.euclidean,
#   # g.hc.bray,
#   g.hc.wnsk1,ncol=1
# )
# 
# 
# 
# pdf("plots/hclust_slide.pdf",width=17,height=10)
# grid.arrange(
#   g.hc.euclidean,
#   g.hc.wnsk1,ncol=1
# )
# dev.off()
# shell.exec("plots/hclust_slide.pdf")
# 
# legend  <- get.tax.legend()
# pdf("plots/legend.pdf",width=4,height=5)
# grid.arrange(legend)
# dev.off()
# shell.exec("plots/legend.pdf")
# 
# 
# # qrcode ------------------------------------------------------------------
# 
# 
# library(qrcode)
# code <- qr_code("https://github.com/TaurLab/ephraim_slamka_tax_distance")
# 
# png("plots/qrcode.png")
# plot(code)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # distribution testing,  chisq test ----------------------------------------------------------
# 
# 
# xx <- pairs2 %>%
#   mutate(euclidean.dist=scales::rescale(euclidean.dist,to=c(0,1))) %>% 
#   pivot_longer(cols=ends_with("dist"),names_to="dist.metric",values_to="dist.value") %>%
#   mutate(dist.value.bin=cut2(dist.value,lower=c(0.25,0.5,0.75)))
# 
# dists <- c("euclidean.dist", "bray.dist", "horn.dist", "unifrac.dist", "taxhorn.dist")
# 
# yy <- combn(dists,2) %>% t() %>% as_tibble() %>%
#   mutate(res=map2(V1,V2,~{
#     # .x="euclidean.dist"
#     # .y="bray.dist"
#     sub <- xx %>% filter(dist.metric %in% c(.x,.y))
#     zz <- with(sub,chisq.test(dist.metric,dist.value.bin))
#     tibble(chisq.stat=zz$statistic,pval=zz$p.value)
#   })) %>% unnest(res)
# 
# 
# yy %>% dt()
# 
# 
# 
# 
# 
# 
# 
