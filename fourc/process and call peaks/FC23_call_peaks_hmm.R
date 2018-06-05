#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

#@ script for Calling 4C peaks using HMM and depMixS4 

# load libraries and functions

library(tidyverse)
library(GenomicRanges)
source("FC20_functions.R")
source("FC21_hmm_functions.R")
source("FC23_additional_functions.R")

#### load required data ####
load("declutligevents.Rdata") # produced using FC22, or use example Rdata object
SAMPLEPATH <- <PATH/TO/SAMPLE/LIST> # see supplementary file 1 in the manuscript. 

all_samples <- read_tsv(SAMPLEPATH)

#### fit to hmm ####

WIN = 5000; STEP = 1000
fitted_clones = fit_all_clones(res, AllSamples)

#### identify and aggregate all intervals  ####

intervalmerge <- list()
keepclonestates <- list()
keepaggintervals <-list()

for(clone in names(fitted_clones)) {
  clonestates <- find_clonestates_for_clone(clone)
  aggintervals <- calculate_agg_intervals(clonestates = clonestates)
  keepaggintervals[[clone]] <- aggintervals
  intervalmerge[[clone]] <- merge_intervals(aggintervals = aggintervals)
  keepclonestates[[clone]] <- choose_data_to_keep(clonestates, aggintervals)
  rm(clonestates, aggintervals)
}  

keepaggintervals <- bind_rows(keepaggintervals, .id = "clone")
keepclonestates <- bind_rows(keepclonestates)
intervalmerge <- bind_rows(intervalmerge, .id = "clone")

#### filter peaks ####

ggplot(intervalmerge) + 
  geom_density(aes(x = width), fill = "blue", alpha = 0.1) + 
  theme_bw() + 
  geom_vline(xintercept = 5e4, linetype = "dotted")

hcpeaks <- intervalmerge %>% 
  filter(width < 5e4) %>% ungroup() %>%
  mutate(wtclone = case_when(grepl("ED", clone) ~ "ED", 
                             grepl("TBX4B ", clone) ~ "TBX4B",
                             TRUE ~ clone)) %>%
  mutate(clone = factor(clone)) %>%
  dplyr::rename(chr = seqnames)

# find self overlaps between peaks to identify 
peaks1.gr <- GRanges(hcpeaks)

ol <- findOverlaps(peaks1.gr)
dups <- ol[which(peaks1.gr[queryHits(ol)]$wtclone != peaks1.gr[subjectHits(ol)]$wtclone)]
dups <- unique(queryHits(dups[peaks1.gr[queryHits(dups)]$peakheight < peaks1.gr[subjectHits(dups)]$peakheight]))

hcpeaks <- hcpeaks[-dups, ] %>%
  mutate(chr = ifelse(chr != 23, paste0("chr", chr), "chrX"), 
         note = "")

is.gr <- AllSamples %>% 
  select(clonename, chr = Chr, position = Position) %>% 
  distinct() %>% 
  filter(position > 0)

is.gr <- GRanges(seqnames = is.gr$chr, IRanges(start = is.gr$position, width = 1))

peaks1.gr <- GRanges(hcpeaks)

ol <- suppressWarnings(findOverlaps(is.gr, peaks1.gr))

hcpeaks$note[subjectHits(ol)] <- "IS peak"











