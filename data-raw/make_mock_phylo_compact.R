
# get a subset of the data (run this) --------------------------------------------------------------

library(tidyverse)
library(yingtools2)
library(phyloseq)
rm(list=ls())
load("data/mock_phylo.RData")

if (FALSE) {
  s <- get.samp(mock.phy)
  #get tyler samps
  dist.tyler <- s %>% filter(source=="tyler") %>% pull(sample) %>% prune_samples(mock.phy) %>%
    prune_unused_taxa() %>% distance("bray")
  samps.tyler <- get.pairwise(dist.tyler) %>% arrange(desc(dist)) %>%
    slice(round(seq(1,n(),length.out=6))) %>%
    {c(.$sample1,.$sample2)} %>%
    unique()
  samps.bmt <- s %>% group_by(pt,pt.day.samp) %>%
    filter(source=="bmt",n()>2) %>%
    summarize() %>%
    slice(1:2) %>%
    select(pt) %>% inner_join(s,by="pt") %>% pull(sample)
  samp.subset <- c(samps.bmt,samps.tyler)
  copy.as.Rcode(samp.subset)
}


samp.subset <- c("FMT.0092A..pool624.FMT.Hiseq", "FMT.0092A..pool627.FMT", "FMT.0092U..pool627.FMT", "FMT.0092A..pool643.Hiseq", "FMT.0092U..pool643.Hiseq",
                 "FMT.0092AA..pool675", "FMT.0092DD..pool675", "FMT.0092EE..pool675", "FMT.0092FF..pool675", "FMT.0092GG..pool675", "FMT.0092B..pool706",
                 "FMT.0092CC..pool706", "FMT.0092C..pool706", "FMT.0092D..pool706", "FMT.0092E..pool706", "FMT.0092F..pool706", "FMT.0092G..pool706",
                 "FMT.0092HH..pool706", "FMT.0092II..pool706", "FMT.0092I..pool706", "FMT.0092J..pool706", "FMT.0092K..pool706", "FMT.0092L..pool706",
                 "FMT.0092M..pool706", "FMT.0092N..pool706", "FMT.0092O..pool706", "FMT.0092P..pool706", "FMT.0092Q..pool706", "FMT.0092R..pool706",
                 "FMT.0092S..pool706", "FMT.0092T..pool706", "FMT.0092W..pool706", "FMT.0092X..pool706", "FMT.0092Y..pool706", "FMT.0092Z..pool706",
                 "FMT.0092JJ..pool709", "FMT.0092KK..pool777", "FMT.0092MM..pool777", "FMT.0092NN..pool802", "FMT.0092OO..pool802", "FMT.0091A..pool624.FMT.Hiseq",
                 "FMT.0091A..pool627.FMT", "FMT.0091A..pool643.Hiseq", "FMT.0091Y..pool643.Hiseq", "FMT.0091B..pool658", "FMT.0091C..pool658",
                 "FMT.0091D..pool658", "FMT.0091E..pool658", "FMT.0091F..pool658", "FMT.0091G..pool658", "FMT.0091H..pool658", "FMT.0091I..pool658",
                 "FMT.0091J..pool658", "FMT.0091AA..pool662", "FMT.0091BB..pool662", "FMT.0091CC..pool662", "FMT.0091DD..pool662", "FMT.0091EE..pool662",
                 "FMT.0091FF..pool662", "FMT.0091GG..pool662", "FMT.0091K..pool662", "FMT.0091O..pool662", "FMT.0091Q..pool662", "FMT.0091V..pool662",
                 "FMT.0091W..pool662", "FMT.0091X..pool662", "FMT.0091HH..pool663.Hiseq", "FMT.0091II..pool663.Hiseq", "FMT.0091KK..pool665.Hiseq",
                 "FMT.0091LL..pool675", "FMT.0091MM..pool675", "FMT.0091NN..pool675", "FMT.0091OO..pool706", "FMT.0091PP..pool706", "FMT.0091QQ..pool706",
                 "FMT.0091S..pool706", "FMT.0091Z..pool706", "FMT.0091RR..pool709", "FMT.0091SS..pool736", "FMT.0091TT..pool756", "FMT.0091UU..pool781",
                 "FMT.0089A..pool624.FMT.Hiseq", "FMT.0089A..pool627.FMT", "FMT.0089A..pool643.Hiseq", "FMT.0089P..pool643.Hiseq", "FMT.0089B..pool661",
                 "FMT.0089E..pool661", "FMT.0089H..pool661", "FMT.0089I..pool661", "FMT.0089J..pool661", "FMT.0089K..pool661", "FMT.0089L..pool661",
                 "FMT.0089M..pool661", "FMT.0089N..pool661", "FMT.0089O..pool663.Hiseq", "FMT.0089Q..pool663.Hiseq", "FMT.0089R..pool663.Hiseq",
                 "FMT.0089S..pool663.Hiseq", "FMT.0089C..pool732", "FMT.0089D..pool732", "FMT.0089F..pool732", "FMT.0089G..pool732", "4A.-20.D3",
                 "10B.RT.D11", "10A.RT.D11", "1A", "5A.-80.D3", "7A.4C.D8", "3A.4C.D3", "3B.4C.D3", "5B.-80.D3", "8B.-20.D8", "8A.-20.D8")

#this is 15 mb instead.
smaller.phy <- prune_samples(samp.subset,mock.phy) %>% prune_unused_taxa()


s <- smaller.phy %>% get.samp() %>%
  arrange(pt,day,Sample_ID) %>% 
  group_by(pt) %>%
  mutate(id=cur_group_id(),
         day=day-min(day)) %>%
  ungroup() %>%
  group_by(Sample_ID) %>%
  mutate(sample.id=paste0("samp",cur_group_id()),
         replicate=row_number()) %>%
  ungroup() %>%
  mutate(
    is.healthy=substr(MRN,1,2)=="ty",
    pt=ifelse(is.healthy,"HV",paste0("PT",id)),
         pt.day=str_glue("{pt}_d{day}"),
         pt.day.samp=str_glue("{pt.day}_{sample.id}"),
         new.sample=str_glue("{pt.day.samp}_{replicate}"))

sample_names(smaller.phy) <- s$new.sample[match(sample_names(smaller.phy),s$sample)]
s.final <- s %>% select(sample=new.sample,pt,day,pt.day,pt.day.samp,sample.id,replicate,source,is.healthy) %>% 
  group_by(pt.day) %>% mutate(has.sameday=n()>1) %>% 
  group_by(pt.day.samp) %>% mutate(has.reseq=n()>1) %>%
  ungroup()
sample_data(smaller.phy) <- s.final %>% set.samp()

library(msa)
library(seqinr)
library(ape)
seq <- refseq(smaller.phy)
msa <- msa::msa(seq)
align <- msaConvert(msa,type="seqinr::alignment")
d <- dist.alignment(align,"identity")
tree <- njs(d)
phy_tree(smaller.phy) <- tree

saveRDS(smaller.phy,file="data/mock_phylo_compact.rds")
# smaller.phy <- readRDS("data/mock_phylo_compact.rds")
# s <- get.samp(smaller.phy)
# new.s <- s

