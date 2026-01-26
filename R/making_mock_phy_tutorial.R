


# phy <- phy.combined.all
# subset_samplenames <- sample_names(phy)[1:5]
# subset_taxanames <- taxa_names(phy)[1:10]
# fewer.samples.phy <- prune_samples(subset_samplenames,phy)
# fewer.samples.and.fewer.taxa <- prune_taxa(subset_taxanames,fewer.samples.phy)
# fewer.samples.and.fewer.taxa

# wizardry with dots
# phy.combined.all %>% get.samp()
# get.samp(phy.combined.all)
# phy.combined.all %>% get.samp(.)

phy <- phy.combined.all %>%
  prune_samples(sample_names(.)[1:2],.) %>%
  prune_taxa(taxa_names(.)[1:10],.)
sample_names(phy) <- paste0("sample",1:2)


## task: replace sample1 and sample2 with sample3, zero abundance 
# make new otu (add sample3, remove sample1 and sample2)
otu <- phy %>% get.otu()
df <- otu %>% as.data.frame()
new_df <- df %>% mutate(sample3=0) %>% select(sample3)
new_otu <- new_df %>% as.matrix()
otable <- otu_table(new_otu,taxa_are_rows=TRUE)
# keep old tax table 
oldttable <- tax_table(phy)
newphy <- phyloseq(otable,oldttable)

newotu <- phy %>% 
  get.otu(as.matrix=FALSE) %>%
  mutate(sample3=0) %>% 
  select(otu,sample3) %>%
  set.otu()
newtax <- phy %>% 
  get.tax() %>%
  select(-full_taxonomy) %>%
  set.tax() %>%
  arrange(taxid)
# rownames(newtax) <- letters[1:10]
newphy <- phyloseq(newotu,newtax)



## task keep phy, but add new bacteria, with 0 abundances
oldotu <- get.otu(phy,as.matrix=FALSE)
newbact <- oldotu %>% 
  mutate(otu=paste("bacteria",row_number()),
         sample1=0,sample2=0)
newotu <- rbind(oldotu,newbact)

oldtax <- get.tax(phy)
newtaxrows <- oldtax %>%
  select(-full_taxonomy) %>%
  mutate(otu=paste("bacteria",row_number()),
         taxid=paste0("x",taxid),
         Species=paste("Species",row_number()),
         Genus=ifelse(Genus=="Clostridium","Clost",Genus),
         Family=as.character(rnorm(n())))
newtax <- bind_rows(oldtax,newtaxrows)
newphy <- phyloseq(set.otu(newotu),set.tax(newtax))


## combining 2 phys with no overlap in bacteria or samples

phy1 <- phy
sample_names(phy1) <- paste0(sample_names(phy1),"A")
taxa_names(phy1) <- paste0(taxa_names(phy1),"A")
phy2 <- phy
sample_names(phy2) <- paste0(sample_names(phy2),"B")
taxa_names(phy2) <- paste0(taxa_names(phy2),"B")


otu1 <- phy1 %>% get.otu(as.matrix=FALSE)
otu2 <- phy2 %>% get.otu(as.matrix=FALSE)

newotu <- bind_rows(otu1,otu2) %>%
  mutate(sample1A=ifelse(is.na(sample1A),0,sample1A),
         sample2A=ifelse(is.na(sample2A),0,sample2A),
         sample1B=ifelse(is.na(sample1B),0,sample1B),
         sample2B=ifelse(is.na(sample2B),0,sample2B))
# faster ways
# newotu <- bind_rows(otu1,otu2) %>%
#   replace_na(list(sample1A=0,sample1B=0,sample2A=0,sample2B=0))
# newotu <- bind_rows(otu1,otu2) %>%
#   mutate(across(.cols=-otu,~coalesce(.,0)))
tax1 <- phy1 %>% get.tax()
tax2 <- phy2 %>% get.tax()
newtax <- bind_rows(tax1,tax2)
newphy <- phyloseq(set.otu(newotu),set.tax(newtax))
newphy %>% get.otu()

# actually do this super fast way
testphy <- merge_phyloseq(phy1,phy2)














