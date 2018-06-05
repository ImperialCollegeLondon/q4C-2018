#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

find_clonestates_for_clone <- function(clonename) {
  
  message(clonename)
  
  dataforclone <- fitted_clones[[clonename]]
  
  # if there was not enough data to call states, move on. 
  if (is.null(dataforclone$clone1$state)) {return(NULL)}
  if (test_clone_states(clonename) == 0) {return(NULL)}
  
  # extract all ligation events per window (including <2)
  winforclone <- get_clone_counts(res = res, clone = dataforclone$clone_name, steps = 1000, win = WIN)
  
  clonestates <- list()  
  for (smnum in 1:2) {
    clonestates[[smnum]] <- 
      subset_fitted_data_for_sample(dataforclone = dataforclone, 
                                    winforclone = winforclone, smnum = smnum) %>%
      select(chr, pos, ligevents, state, counts1) 
  }
  clonestates <- bind_rows(clonestates, .id = "smnum") %>%
    mutate(clone = clonename)
  
  return(clonestates)
  
}

test_clone_states <- function(clonename){
  
  dataforclone <- fitted_clones[[clonename]]
  
  a <- bind_cols(dataforclone$clone1$data, dataforclone$clone1$state) %>%
    filter(state != 1) %>% select(chr, pos)
  
  b <- bind_cols(dataforclone$clone2$data, dataforclone$clone2$state) %>%
    filter(state != 1) %>% select(chr, pos)
  
  return(inner_join(a, b, by = c("chr", "pos")) %>% nrow())
  
}

subset_fitted_data_for_sample <- function(dataforclone, winforclone, smnum = 1){
  # retrieve data from hmm fitted object
  data1 <- dataforclone[[paste0("clone", smnum)]]$data
  state1 <- dataforclone[[paste0("clone", smnum)]]$state
  data1 <- cbind(data1, state1) 
  
  # retrieve complete window information
  wins_sam1 <- winforclone[[paste0("wins_sam", smnum)]][-c(24, 25)] %>%
    lapply(FUN = function(x) 
      x <- as.data.frame(cbind(pos = as.numeric(names(x)), ligevents = x))) %>%
    bind_rows(.id = "chr") %>% 
    mutate(chr = as.numeric(chr)) %>% arrange(chr, pos)
  
  # merge the two to add the <2 read windows back into the hmm data. These will be state = 4
  clonestates <- left_join(wins_sam1, data1, by = c("chr", "pos"))
  ah <- which(!is.na(clonestates$state))
  clonestates$state[-ah] <- 0
  
  # reorder the states by frequency. 
  clonestates$state <- factor(clonestates$state, 
                              levels = rev(names(table(clonestates$state))[order(table(clonestates$state))]), 
                              labels = (0:3)[1:length(unique(clonestates$state))])
  
  # in case not all states are represented, add them into factor labels
  clonestates$state <- forcats::fct_expand(clonestates$state, 0:3)
  
  clonestates
}

calculate_agg_intervals <- function(clonestates){
  
  if(is.null(clonestates)) {return(NULL)}
  
  aggintervals <- list()
  
  for (chrtodo in 1:23) {
    message(paste0(chrtodo, " "))
    aggintervals[[chrtodo]] <- list()
    clonechr <- clonestates[clonestates$chr == chrtodo, ]
    check <- count(clonechr, smnum, state) %>% 
      mutate(state = as.numeric(as.character(state))) %>% 
      group_by(smnum) %>% 
      summarise(max = max(state)) %>%
      summarise(min = min(max)) %>% .[["min"]]
    if (check < 2) {next} 
    agg_states <- calculate_agg_states(clonechr)
    
    for (smnumtodo in 1:2) {
      aggintervals[[chrtodo]][[smnumtodo]] <- 
        find_top_intervals_per_chr(clonechr[clonechr$smnum == smnumtodo, ], 
                                   valuestoadd = agg_states)
    }
    
    aggintervals[[chrtodo]] <- bind_rows(aggintervals[[chrtodo]], .id = "smnum")
  }
  
  aggintervals <- bind_rows(aggintervals) %>%
    mutate(col = "potential peaks")
  
  if (nrow(aggintervals) > 0) {   
    return(aggintervals)
  }else{return(NULL)}
}

calculate_agg_states <- function(clonechr, cutoff = 2){
  # calculate a window score based on both replicates. 
  clonechrjoin <- inner_join(clonechr[clonechr$smnum == 1,], 
                             clonechr[clonechr$smnum == 2,], 
                             by = c("clone", "chr", "pos"), suffix = c("1", "2")) %>%
    mutate(state1 = as.numeric(as.character(state1)), 
           state2 = as.numeric(as.character(state2))) %>%
    group_by(clone, chr, pos) %>% 
    summarise(minstate = min(state1, state2)) %>%
    filter(minstate >= cutoff)
  
  return(clonechrjoin)
}

find_top_intervals_per_chr <- function(clonechr, splineval = 0.25, valuestoadd = agg_states){
  
  # a - clonestates data for a particular sample, for a given chromosome.
  a <- clonechr
  
  # spline function used to fit all windows > 2. skip if no such windows exist. 
  spline.a <- spline(x = a$pos, y = a$ligevents, n = splineval * length(a$pos))
  valley <- which(diff(sign(diff(spline.a$y))) == 2) + 1

  # split ligation data to intervals by the valleys called by the spline. 
  a$interval <- as.numeric(as.character(cut(
    a$pos, 
    breaks = c(-Inf, spline.a$x[valley], Inf), 
    labels = c(0:length(valley)),
    include.lowest = TRUE)))
  
  localmaxlist <-  a %>% 
    group_by(smnum, interval) %>%
    filter(ligevents == max(ligevents)) %>% 
    dplyr::slice(1) %>% ungroup() %>% 
    inner_join(valuestoadd, by = c("clone", "chr", "pos"))
  
  # keepintervals - keep those intervals with more than one ligation event level (exclude skyscrapers)
  # keepintervals <- a %>%
  #   group_by(smnum, interval) %>%
  #   inner_join(valuestoadd, by = c("clone", "chr", "pos")) %>%
  #   count(ligevents) %>%
  #   count() %>% 
  #   filter(nn > 1) %>% .[["interval"]]
  
  intervals <- tibble(
    id = 0:length(valley), 
    start = c(1, spline.a$x[valley]),
    end = c(spline.a$x[valley], max(a$pos))) %>% 
    inner_join(localmaxlist, by = c("id" = "interval")) %>%
    # filter(id %in% keepintervals) %>% 
    select(smnum, id, chr, start, end, peakpos = pos, peakheight = ligevents)
  
  return(intervals)
}  

merge_intervals <- function(aggintervals) {
  
  if(is.null(aggintervals)) {return(NULL)}
  
  peaks1.gr <- aggintervals %>% filter(smnum == 1) %>% GRanges()
  peaks2.gr <- aggintervals %>% filter(smnum == 2) %>% GRanges()
  
  intervalmerge <- bind_rows(as.data.frame(pintersect(findOverlapPairs(peaks1.gr, peaks2.gr))), 
                             as.data.frame(pintersect(findOverlapPairs(peaks2.gr, peaks1.gr))))
  
  intervalmerge <- 
    intervalmerge %>%
    filter(end - start > 1) %>%
    group_by(start) %>%
    filter(peakpos > start & peakpos < end) %>%
    count(start, end) %>% 
    filter(n > 1) %>%
    inner_join(intervalmerge, by = c("start", "end")) %>%
    mutate(col = "common peaks") %>%
    ungroup()
}

choose_data_to_keep <- function(clonestates, aggintervals) {
  chrstodo <- unique(aggintervals$chr)
  res <- list()
  for (chrtodo in chrstodo) {
    res[[chrtodo]] <- clonestates %>%
      filter(chr == chrtodo, pos > min(aggintervals$start) - 3e6, pos < max(aggintervals$end) + 3e6)
  }
  res <- res %>% bind_rows()
  return(res)
}

