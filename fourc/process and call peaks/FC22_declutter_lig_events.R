#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#@ script for declulttering ligation event data from q4C experiment - correct barcode errors

# load libraries and functions
library(tidyverse)
library(gridExtra)
source(file = "FC20_functions.R")

# run starts here

# load data
dataset.path <- <PATH/TO/EXTRACTED/LIG/EVENTS>
is_chr <- "chr22"         # for example	
is_position <- 44323198   # for example

  
ligeventfiles <- list.files( path = dataset.path, pattern = "ligevent", full.names = TRUE )

ligeventsall <- list()
for(i in 1:length(ligeventfiles) )
{
  samplename <- gsub(".+/|_.+", "", ligeventfiles[i])
  ligeventsall[[samplename]] <- read_tsv(ligeventfiles[i])
}
ligeventsall <- bind_rows(ligeventsall, .id = "samplename")
rm(ligeventfiles, i, samplename, dataset.path)  


# declutter by event ID (can also delutter by ligation site unless all come from the same original clone)
ligeventsdeclut <- ligeventsall %>%
  group_by(EventID) %>%
  filter(ReadCount == max(ReadCount)) %>%
  ungroup() %>% arrange(samplename, EventID)

ligeventsdeclut <- ligeventsdeclut %>%
  filter(!EventID %in% (count(ligeventsdeclut, EventID) %>% filter(n > 1) %>% .[["EventID"]]))

ligeventsdeclut <- ligeventsdeclut %>%
  mutate(EventID = str_replace_all(EventID, "\\*\\*\\*", "\\*_\\*_\\*")) %>% # correct to allow three columns
  separate(EventID, c( "lig.chr", "lig.pos", "lig.strand", "shear.chr", "shear.pos", "shear.strand"), 
           sep = ":|_", remove = FALSE)

ligeventsdeclutclean <- correct.lig.site(ligeventsdeclut, is_chr, is_position)
save(ligeventsall, ligeventsdeclutclean, file = "declutligevents.Rdata")



