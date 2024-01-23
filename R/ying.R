
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
library(rlang)
library(patchwork)

rm(list=ls())
phy <- readRDS("data/mock_phylo_compact.rds")
# source("R/functions.R")
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
    pt.day.samp1==pt.day.samp2 ~ "(4) Re-seq of\nsame samp",
    pt.day1==pt.day2 ~ "(3) Same subj/day,\ndiff samps",
    pt1==pt2 ~ "(2) Same subj,\ndiff days/samps",
    TRUE ~ "(1) Different\nsubjects",
  ),
  status=fct_reordern(status,status))

rm(s1,s2)


# calc dists (run) --------------------------------------------------------------


calc.tax.dpcoa.distance <- function(phy) {
  phy2 <- phy
  t <- get.tax(phy2)
  form <- as.formula(paste0("~",paste(c(rank_names(phy),"otu"),collapse="/")))
  tr <- as.phylo.formula2(x=form,data=t)
  phy_tree(phy2) <- tr
  distance(phy2,"dpcoa")
}

calc.distance2 <- function(.phy, method, rarefy=FALSE, use.pct=FALSE, use.tax.tree=FALSE, ...) {
  
  if (rarefy) {
    .phy <- rarefy_even_depth(.phy)
  }
  if (use.pct) {
    .phy <- transform_sample_counts(.phy, function(x) x / sum(x) )
  }
  if (use.tax.tree) {
    form <- as.formula(paste("~",paste(c(rank_names(.phy),"otu"),collapse="/")))
    tax <- get.tax(.phy)
    tax.tree <- as.phylo.formula2(form,data=tax)
    phy_tree(.phy) <- tax.tree
  }

  if (rlang::is_function(method)) {
    dist <- method(phy)
  } else if (is.character(method)) {
    if (method=="taxhorn") {
      dist <- calc.taxhorn.distance(phy)
    } else if (method=="tax.dpcoa") {
      dist <- calc.tax.dpcoa.distance(phy)
    } else {
      dist <- distance(physeq=phy, method=method, ...)
    }
  } else {
    stop("unknown method")
  }
  return(dist)
}


check_same <- function(method) {
  phy.rel <- transform_sample_counts(.phy, function(x) x / sum(x) )
  d1 <- distance(phy,method)
  d2 <- distance(phy.rel,method)
  p1 <- d1 %>% get.pairwise()
  p2 <- d2 %>% get.pairwise()
  compare <- inner_join(p1,p2,by=c("sample1","sample2")) %>%
    mutate(diff=dist.x-dist.y)
  
  g <- ggplot(compare,aes(x=dist.x,y=dist.y)) + geom_point()
  print(g)
  all(abs(compare$diff) < 1e-10)
}

# check_same("euclidean")
# check_same("bray")
# check_same("chao")


phy.unfold.taxranks <- function(phy) {
  # samp <- sample_data(phy)
  # phy <- phyloseq(otu_table(phy),tax_table(phy))
  # phy <- yingtools2:::make_taxonomy_distinct.phyloseq(phy,add.rank=TRUE)
  ranks <- rank_names(phy)
  # list of phyloseqs
  phy.levels <- ranks %>% seq_along() %>%
    map(~ranks[1:.x]) %>% map(~{
      phy_ <- phy.collapse(phy,taxranks=.x)
      t <- get.tax(phy_) %>%
        mutate(newotu=paste(!!!syms(.x),sep="|"))
      taxa_names(phy_) <- t$newotu
      return(phy_)
    }) %>%
    setNames(ranks) %>% c(list("asv"=phy))
  # all.levels <- names(phy.levels) # all ranks including 'asv'
  all.tax <- phy.levels %>% map(get.tax) %>% list_rbind()
  all.otu <- phy.levels %>% map(get.otu,as.matrix=FALSE) %>% list_rbind()
  if (anyDuplicated(all.tax$otu)!=0) {
    stop("YTError: There were duplicated taxa names when unfolding!")
  }
  phy.unfold <- phyloseq(set.otu(all.otu),set.tax(all.tax),sample_data(phy))
  n.taxa <- phy.levels %>% map_int(ntaxa) %>% rev()
  n.taxa.text <- n.taxa %>% paste(collapse=" + ") %>% paste0(" = ",nrow(all.tax)," taxa")
  message(str_glue("Created new unfolded phyloseq:\n{n.taxa.text}"))
  return(phy.unfold)
}

calc.mean.distance <- function(phy,method,weights=c(0,7,6,5,4,3,2,1)) {
  # method="horn"
  # weights <- c(0,7,6,5,4,3,2,1)
  phy <- phyloseq(otu_table(phy),tax_table(phy))
  ranks <- rank_names(phy)
  # list of phyloseqs
  phy.levels <- ranks %>% seq_along() %>%
    map(~ranks[1:.x]) %>% map(~phy.collapse(phy,taxranks=.x)) %>%
    setNames(ranks) %>% c(list("asv"=phy))
  all.levels <- names(phy.levels) # all ranks including 'asv'
  # calculate the distance matrix (metric=method) for each level.
  # this is a list of distance matrices.
  dist.levels <- phy.levels %>% map(~distance(.x,method=method))
  # run get.pairwise() to get a list of pairwise distances.
  wts <- tibble(taxlevel=all.levels,weight=weights)
  pairwise.all <- dist.levels %>% 
    map(get.pairwise) %>%
    list_rbind(names_to = "taxlevel") %>%
    left_join(wts,by="taxlevel")
  pairwise.calcdist <- pairwise.all %>% group_by(sample1,sample2) %>%
    summarize(dist=sum(dist*weight)/sum(weight),
              .groups="drop")
  taxdist <- get.dist(pairwise.calcdist)
  taxdist
}


calc.unfold.distance <- function(phy,method) {
  phy.unfold <- phy.unfold.taxranks(phy)
  distance(phy.unfold,method=method)
}





dists <- list(euclidean = calc.distance2(phy,"euclidean"), 
              euclidean.pct = calc.distance2(phy,"euclidean",use.pct=TRUE), # should be normalized
              manhattan.pct = calc.distance2(phy,"manhattan"),
              manhattan = calc.distance2(phy,"manhattan",use.pct=TRUE), # should be normalized
              bray = calc.distance2(phy,"bray"),
              bray.pct = calc.distance2(phy,"bray",use.pct=TRUE),
              chao = calc.distance2(phy,"chao"),
              chao.pct = calc.distance2(phy,"chao",use.pct=TRUE),
              cao = calc.distance2(phy,"cao"),
              cao.pct = calc.distance2(phy,"cao",use.pct=TRUE),
              altGower = calc.distance2(phy,"altGower"),
              altGower.pct = calc.distance2(phy,"altGower",use.pct=TRUE),
              canberra = calc.distance2(phy,"canberra"),
              canberra.pct = calc.distance2(phy,"canberra",use.pct=TRUE),
              kulczynski = calc.distance2(phy,"kulczynski"),
              kulczynski.pct = calc.distance2(phy,"kulczynski",use.pct=TRUE),
              gower = calc.distance2(phy,"gower"),
              gower.pct = calc.distance2(phy,"gower",use.pct=TRUE),
              unifrac = calc.distance2(phy,"unifrac",rarefy=TRUE),
              unifrac.pct = calc.distance2(phy,"unifrac",rarefy=TRUE,use.pct=TRUE),
              wunifrac = calc.distance2(phy,"wunifrac"),
              wunifrac.pct = calc.distance2(phy,"wunifrac",use.pct=TRUE),
              tax.wunifrac = calc.distance2(phy,"wunifrac",use.tax.tree=TRUE),
              tax.wunifrac = calc.distance2(phy,"wunifrac",use.tax.tree=TRUE,use.pct=TRUE),
              tax.dpcoa=calc.distance2(phy,"tax.dpcoa"),
              tax.dpcoa.rel = calc.distance2(phy,"tax.dpcoa",use.pct=TRUE),
              horn = calc.distance2(phy,"horn"),
              unfoldhorn = calc.unfold.distance(phy,"horn"),
              tax.horn = calc.mean.distance(phy,"horn"),
              taxhorn = calc.distance2(phy,"taxhorn"))

# dists$tax.horn <- calc.mean.distance(phy,"horn")
# dists$unfold.horn <- calc.unfold.distance(phy,"horn")



# violin ------------------------------------------------------------------

view_violin <- function(dist,title) {
  pairdata <- dist %>% get.pairwise() %>%
    inner_join(pairs,by=c("sample1","sample2"))
  ggplot(pairdata,aes(x=status,y=dist,fill=status)) + geom_violin() +
    geom_boxplot(alpha=0.2) + ggtitle(title) + expand_limits(y=c(0,1)) + theme(legend.position="none") +
    scale_y_continuous(str_glue("{title}\ndistance")) +
    scale_x_discrete("")
}

glist.violin <- dists %>% imap(~{
  view_violin(.x,.y)
})

inject(wrap_plots(!!!glist.violin,ncol=6))

# show reseqs -------------------------------------------------------------



make_species_taxa_names <- function(phy) {
  level <- rank_names(phy) %>% last() %>% sym()
  tax <- get.tax(phy) %>%
    group_by(!!level) %>%
    mutate(newotu_=paste(!!level,row_number())) %>%
    ungroup()
  taxa_names(phy) <- tax$newotu_
  return(phy)
}





s.reseq <- s %>% filter(has.reseq)
phy.reseq <- phy %>% filter(has.reseq) %>% make_species_taxa_names()
otu.reseq <- phy.reseq %>% get.otu.melt()
g.reseq.stack <- ggplot(otu.reseq,aes(x=sample,y=pctseqs,fill=otu,label=Species)) + 
  geom_taxonomy(label.split = TRUE) + 
  facet_grid(. ~ pt.day.samp,space="free_x",scales="free_x") +
  theme(axis.text.x=element_text(angle=-90,hjust=0))

g.reseq.stack

ggsave("plots/ying_reseq_samples.pdf",g.reseq.stack,width=16,height=10)
shell.exec("plots/ying_reseq_samples.pdf")

xx <- phy.reseq %>% 
  # phy.collapse() %>%
  yingtools2:::make_taxonomy_distinct.phyloseq() %>%
  get.otu.melt() %>%
  mutate(otu=fct_reordern(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species,otu),
         taxfacet=as.character(fct_lump_n(Phylum,5)))

yy <- xx %>% count(sample,pt.day.samp)

g.reseq.bar <- ggplot(xx) +
  geom_col(aes(x=otu,y=pctseqs,fill=otu),show.legend=FALSE) +
  # geom_text(data=yy,aes(x=0,y=Inf,label=pt.day.samp,color=pt.day.samp),hjust="inward",vjust="inward") +
  scale_y_continuous(trans=log_epsilon_trans(),limits=c(0,1)) +
  scale_fill_taxonomy(data=xx,fill=otu) +
  theme(axis.text.x=element_text(angle=-90,hjust=0),strip.text.y = element_text(angle=0)) +
  ggnewscale::new_scale_fill() +
  geom_rect(data=yy,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,fill=pt.day.samp),alpha=0.1,show.legend=FALSE) +
  scale_fill_discrete() +
  facet_grid(sample ~ taxfacet,scales="free_x",space="free_x")

g.reseq.bar
ggsave("plots/ying_reseq_bar.pdf",g.reseq.bar,width=80,height=15,limitsize = FALSE)
shell.exec("plots/ying_reseq_bar.pdf")

s.sameday <- s %>% filter(has.sameday)
phy.sameday <- phy %>% filter(has.sameday) %>% make_species_taxa_names()
otu.sameday <- phy.sameday %>% get.otu.melt()
g.sameday.stack <- ggplot(otu.sameday,aes(x=sample,y=pctseqs,fill=otu,label=Species)) + 
  geom_taxonomy(label.split = TRUE) + 
  facet_grid(. ~ pt.day,space="free_x",scales="free_x") +
  theme(axis.text.x=element_text(angle=-90,hjust=0))

# g.reseq.stack / g.sameday.stack

ggsave("plots/ying_sameday_samples.pdf",g.sameday.stack,width=16,height=10)
shell.exec("plots/ying_sameday_samples.pdf")


xx2 <- phy.sameday %>% 
  # phy.collapse() %>%
  yingtools2:::make_taxonomy_distinct.phyloseq() %>%
  get.otu.melt() %>%
  mutate(otu=fct_reordern(otu,Superkingdom,Phylum,Class,Order,Family,Genus,Species,otu),
         taxfacet=as.character(fct_lump_n(Phylum,5)))
yy2 <- xx2 %>% count(sample,pt.day)

g.sameday.bar <- ggplot(xx2) +
  geom_col(aes(x=otu,y=pctseqs,fill=otu),show.legend=FALSE) +
  # geom_text(data=yy,aes(x=0,y=Inf,label=pt.day.samp,color=pt.day.samp),hjust="inward",vjust="inward") +
  scale_y_continuous(trans=log_epsilon_trans(),limits=c(0,1)) +
  scale_fill_taxonomy(data=xx2,fill=otu) +
  theme(axis.text.x=element_text(angle=-90,hjust=0),strip.text.y = element_text(angle=0)) +
  ggnewscale::new_scale_fill() +
  geom_rect(data=yy2,aes(xmin=-Inf,xmax=Inf,ymin=-Inf,ymax=Inf,fill=pt.day),alpha=0.1,show.legend=FALSE) +
  scale_fill_discrete() +
  facet_grid(sample ~ taxfacet,scales="free_x",space="free_x")
ggsave("plots/ying_sameday_bar.pdf",g.sameday.bar,width=80,height=15,limitsize = FALSE)
shell.exec("plots/ying_sameday_bar.pdf")


# view hclust -------------------------------------------------------------


view.hclust <- function(dist,physub,title="") {
  
  # physub
  # dist <- dists$manhattan
  # dist <- calc.distance(physub,dist)
  hc <- hclust(dist)
  tr <- as.phylo(hc)

  droptips <- setdiff(tr$tip.label,sample_names(physub))
  tr <- drop.tip(tr,droptips)
  gt <- ggtree(tr) %<+% get.samp(physub)

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
  otu <- physub %>% get.otu.melt() %>%
    left_join(map,by="sample")

  g.tax <- ggplot(otu,aes(x=y,y=pctseqs,fill=Species,label=Species)) +
    geom_taxonomy() +
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
    scale_x_reverse() +
    coord_flip(ylim=xlim,expand=FALSE) +
    ggtitle(title)


  gg.stack2(g.hclust,g.groups,g.tax,heights=c(1.5,1,4))
  # gg.stack2(g.hclust,g.tax,heights=c(1.5,4),align.xlim=FALSE)
}

physub <- phy %>% subset_samples(pt!="PT3" & day<100)
physub <- phy %>% filter(has.sameday)

# view.hclust(dists$manhattan,physub)
# view.hclust(dists$tax.dpcoa,physub)

glist.hclust <- dists %>% imap(~{
    view.hclust(.x,physub,.y)
  })


pdf("plots/ying_hclust.pdf",width=20,height=12)
glist.hclust
dev.off()

shell.exec("plots/ying_hclust.pdf")





# manual calc -------------------------------------------------------------



lookup <- function(dist,s1,s2) {
  as.matrix(dist)[s1,s2]
}

comparedist <-function(samp1,samp2,pct=FALSE,eps=NULL,.phy=phy) {
  # samp1 <- "HV_d0_samp1_2";samp2 <- "HV_d0_samp1_7"; .phy=phy; eps=NULL;pct=TRUE
  physub <- prune_samples(c(samp1,samp2),.phy) %>% 
    prune_unused_taxa()
  tax <- get.tax(physub)
  otu <- physub %>% get.otu(as.matrix=FALSE) %>%
    rename(sample1=!!sym(samp1),sample2=!!sym(samp2))
  n1 <- sum(otu$sample1)
  n2 <- sum(otu$sample2)
  if (pct) {
    otu <- otu %>% mutate(sample1=sample1/sum(sample1),
                          sample2=sample2/sum(sample2))
  }
  otu <- otu %>% 
    mutate(otu=fct_reordern(otu,sample1,sample2),
           diff=abs(sample1-sample2),
           diff.squared=diff^2,
           sum=sample1+sample2,
           product=sample1*sample2)
  
  quants <- otu %>%
    select(-sample1,-sample2) %>%
    pivot_longer(cols=-otu)
  dcalc <- otu %>% 
    summarize(manhattan=sum(diff),
              bray=sum(diff)/(sum(sum)),
              euclidean=sqrt(sum(diff.squared)),
              lambda1=sum(sample1^2)/sum(sample1)^2,
              lambda2=sum(sample2^2)/sum(sample2)^2,
              horn=1-(2*sum(product)/((lambda1+lambda2)*sum(sample1)*sum(sample2))))
distance()
  title <- str_glue("manhattan={pretty_number(dcalc$manhattan)}; bray={pretty_number(dcalc$bray)}; euclidean={pretty_number(dcalc$euclidean)}; horn={pretty_number(dcalc$horn)}")
  if (is.null(eps)) {
    eps <- c(otu$sample1,otu$sample2) %>% quantile(probs=0.4) %>% as.vector() %>% log(base=10) %>% round() %>% {10^.}  
  }

  g.abundances <- ggplot(otu) +
    # geom_col(aes(x=otu,y=sample1),fill="orange") +
    # geom_col(aes(x=otu,y=-sample2),fill="steelblue") +
    geom_col(aes(x=otu,y=sample1,fill=otu),show.legend=FALSE) +
    geom_col(aes(x=otu,y=-sample2,fill=otu),show.legend=FALSE) +
    annotate("text",x=-Inf,y=Inf,label=str_glue("sample1={pretty_number(n1)}"),hjust="inward",vjust="inward") +
    annotate("text",x=-Inf,y=-Inf,label=str_glue("sample2={pretty_number(n2)}"),hjust="inward",vjust="inward") +
    geom_hline(yintercept=0) +
    theme(axis.text.x=element_blank()) +
    scale_fill_taxonomy(data=tax,fill=otu) +
    scale_y_continuous(trans=log_epsilon_trans(eps))
  g.quants <- ggplot(quants,aes(x=otu,y=value,fill=name)) + geom_col() +
    facet_grid(name ~ .,scales="free_y") +
    ggtitle(title)

  gg.stack2(g.quants,g.abundances)
}

# same samp
# comparedist("HV_d0_samp1_2","HV_d0_samp1_7",pct=FALSE)
comparedist("HV_d0_samp1_2","HV_d0_samp1_7",pct=TRUE)

# same samp, but terrible absolute bray
# comparedist("HV_d0_samp1_4","HV_d0_samp1_6",pct=FALSE)
comparedist("HV_d0_samp1_4","HV_d0_samp1_6",pct=FALSE)
comparedist("HV_d0_samp1_4","HV_d0_samp1_6",pct=TRUE)

# same samp, but terrible bray
# comparedist("HV_d0_samp1_1","HV_d0_samp1_6",pct=FALSE)
comparedist("HV_d0_samp1_1","HV_d0_samp1_6",pct=TRUE)


# same day but more diffs
comparedist("PT2_d0_samp24_1","PT2_d0_samp22_1",pct=TRUE)


comparedist("PT2_d9_samp40_1","PT2_d0_samp22_1",pct=TRUE,eps=0.001)



# view adjacent dists in a patient ----------------------------------------

plot.dist <- function(title,dst,yrow=1.1,ss,color="pink") {
  pw <- get.pairwise(dst) %>%
    mutate(pos1=ss$pos[match(sample1,ss$sample)],
           pos2=ss$pos[match(sample2,ss$sample)]) %>%
    filter(!is.na(pos1),!is.na(pos2)) %>%
    mutate(x.mid=midpoint(pos1,pos2),
           size=scales::rescale(dist,from = c(0,pmax(1,max(dist))),to=c(0,1)),
           dist.lbl=pretty_number(dist)) %>%
    filter(abs(pos1-pos2)==1)
  list(geom_vline(data=pw,aes(xintercept=x.mid),
                  alpha=0.8,
                  color="gray"),
       geom_point(data=pw,aes(x=x.mid,y=yrow,size=size),color=color,show.legend=FALSE,alpha=0.8),
       geom_text(data=pw,aes(x=x.mid,y=yrow,label=dist.lbl),angle=90),
       annotate("text",x=-1,y=yrow,label=title))  
}


plot.pt <- function(p="PT1") {
  ssub <- s %>% filter(pt==p) %>%
    arrange(day,sample) %>%
    mutate(pos=1:n())
  physub <- ssub$sample %>% 
    prune_samples(phy) %>%
    phy.collapse()
  otu <- physub %>% get.otu.melt() %>%
    left_join(select(ssub,sample,pos),by="sample")
  dups <- ssub %>% group_by(pt.day) %>% filter(n()>1)
  g.tax <- ggplot() +
    plot.dist("Unifrac",dists$unifrac.pct,yrow=1.5,color="darkolivegreen3",ss=ssub) +
    plot.dist("Euclidean",dists$euclidean.pct,yrow=1.4,color="steelblue",ss=ssub) +
    plot.dist("Bray",dists$bray.pct,yrow=1.3,color="pink",ss=ssub) +
    plot.dist("Horn",dists$horn,yrow=1.2,color="darkred",ss=ssub) +
    plot.dist("unfold-Horn",dists$unfold.horn,yrow=1.1,color="red",ss=ssub) +
    geom_taxonomy(data=otu,aes(x=pos,y=pctseqs,fill=Species,label=Species),
                  tax.palette = yt.palette2) +
    geom_text(data=ssub,aes(x=pos,y=0,label=day),vjust=1.15) +
    geom_bracket(data=dups,aes(x=pos-0.45,xend=pos+0.45,y=-0.03,group=pt.day,label=""),flip=TRUE,tip = "round",bracket.width = unit(0.35,'char')) +
    geom_bracket(data=dups,aes(x=pos-0.45,xend=pos+0.45,y=1.02,group=pt.day,label=""),tip = "round",bracket.width=unit(0.35,'char')) +
    scale_size_continuous(range=c(1,10)) +
    theme(axis.text.x=element_blank())
  g.tax  
}


plot.pt("PT1")
plot.pt("HV")
plot.pt("PT2")
plot.pt("PT3")





#############
comparedist("HV_d0_samp1_6","HV_d0_samp1_1",pct=TRUE)




# different people, both diverse
comparedist("PT3_d0_samp3_3","HV_d0_samp1_1",pct=TRUE,eps=0.001)




# totally different
# comparedist("HV_d0_samp1_2","PT1_d18_samp84_1",pct=FALSE,eps=100)
comparedist("HV_d0_samp1_2","PT1_d18_samp84_1",pct=TRUE,eps=0.001)

# different people but similar in distance
comparedist("PT2_d24_samp59_1","PT1_d22_samp90_1",pct=TRUE,eps=0.001)



lookup(dists$manhattan,"HV_d0_samp1_1","HV_d0_samp1_6")
lookup(dists$bray.rel,"HV_d0_samp1_1","HV_d0_samp1_6")
lookup(dists$horn.rel,"HV_d0_samp1_1","HV_d0_samp1_6")





dists$taxhorn %>% get.pairwise() %>%
  inner_join(pairs,by=c("sample1","sample2")) %>% 
  filter(status %like% "1") %>%
  arrange(desc(dist))

  


dists$horn.rel %>% get.pairwise() %>%
  # filter(sample1 %like% "HV_d0_samp1",
  #        sample2 %like% "HV_d0_samp1") %>%
  arrange(desc(dist))



otu <- prune_samples(c(samp1,samp2),phy) %>% 
  prune_unused_taxa() %>% get.otu(as.matrix=FALSE) %>%
  rename(sample1=!!sym(samp1),sample2=!!sym(samp2)) %>% 
  mutate(sample1.pct=sample1/sum(sample1),
         sample2.pct=sample2/sum(sample2))


dcalc <- otu %>% mutate(diff=abs(sample1-sample2),
                        diff.squared=diff^2,
                        sum=sample1+sample2,
                        diff.pct=abs(sample1.pct-sample2.pct),
                        sum.pct=sample1.pct+sample2.pct) %>% 
  summarize(manhattan=sum(diff),
            bray=sum(diff)/(sum(sum)),
            euclidean=sqrt(sum(diff.squared))) %>%
  ungroup()

dcalc
lookup(dists$manhattan,samp1,samp2)
lookup(dists$bray,samp1,samp2)
lookup(dists$euclidean,samp1,samp2)




otu %>%
  mutate(otu=fct_reordern(otu,sample1.pct)) %>% 
  ggplot() +
  geom_col(aes(x=otu,y=sample1.pct),fill="orange") +
  # geom_col(aes(x=otu,y=sample1)) +
  geom_col(aes(x=otu,y=-sample2.pct),fill="steelblue") +
  theme(axis.text.x=element_blank()) +
  scale_y_continuous(trans=log_epsilon_trans(0.001))



# # create various distances, including the customized ------------------------------------------------
# 
# # usual distances 
# dist.manhattan <- distance(phy,"manhattan")
# dist.bray <- distance(phy,"bray")
# dist.euclidean <- distance(phy,"euclidean")
# dist.horn <- distance(phy,"horn")
# dist.unifrac <- distance(phy,"unifrac")
# dist.wunifrac <- distance(phy,"wunifrac")
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
# 
# # view some of the pairs to get an idea ----------------------------------------------------------
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





