library(patchwork)
library(rlang)

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
    imap(get.pairwise) %>%
    list_rbind(names_to = "taxlevel") %>%
    left_join(wts,by="taxlevel")
  pairwise.calcdist <- pairwise.all %>% group_by(sample1,sample2) %>%
    summarize(dist=sum(dist*weight)/sum(weight),
              .groups="drop")
  taxdist <- get.dist(pairwise.calcdist)
  taxdist
}



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

calc.unfold.distance <- function(phy,method) {
  phy.unfold <- phy.unfold.taxranks(phy)
  distance(phy.unfold,method=method)
}

nphy <- transform_sample_counts(phy, function(x) x / sum(x) )

dist.list <- list(
  "horn"=distance(nphy,"horn"),
  "bray"=distance(nphy,"bray"),
  "euclidean"=distance(nphy,"euclidean"),
  "tax.horn"=calc.mean.distance(nphy,"horn",weights=c(1,1,1,1,1,1,1,1)),
  "tax.bray"=calc.mean.distance(nphy,"bray",weights=c(1,1,1,1,1,1,1,1)),
  "tax.euclidean"=calc.mean.distance(nphy,"euclidean",weights=c(1,1,1,1,1,1,1,1)),
  "unfold.horn"=calc.unfold.distance(nphy,"horn"),
  "unfold.bray"=calc.unfold.distance(nphy,"bray"),
  "unfold.euclidean"=calc.unfold.distance(nphy,"euclidean")
)


view_violin <- function(dist,title) {
  pairdata <- dist %>% get.pairwise() %>%
    inner_join(pairs,by=c("sample1","sample2"))
  ggplot(pairdata,aes(x=status,y=dist,fill=status)) + geom_violin(position = position_dodge(10)) +
    geom_boxplot(alpha=0.1) + ggtitle(title) + expand_limits(y=c(0,1)) +
    theme(legend.position= "none") 
}


violin.list <- dist.list %>% imap(~{view_violin(.x,.y)})
glist <- inject(wrap_plots(!!!violin.list,ncol=3))
glist









