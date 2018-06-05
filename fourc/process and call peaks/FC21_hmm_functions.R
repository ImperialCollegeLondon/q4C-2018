#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018

# HMM calling functions - written by Tomas Fitzgerald


# function to extract all infected chrosmomes across clones
# @ returns lists of read counts in windows across infected chrosmomes
get_infected_counts <- function(res, AllSamples) {
  pin = 1
  infe_list = NULL
  model_list = list()
  clones = unique(res$clonename)
  print("getting infected counts...")
  for (x in 1:length(clones)) {
    print(paste("clone", x))
    a = AllSamples[AllSamples$clonename == clones[x],]
    if (length(grep("chr", a$Chr))>0) {
      r = res[res$clonename==clones[x] & res$lig.chr==unique(a$Chr),]
      b = get_clone_counts(r, clone=clones[x])
      r$clonename = "control"
      infe_list = rbind(infe_list, r)
      model_list[[pin]] = b
      pin = pin+1
    }
  }
  infe_list$samplename = "control"
  refcounts = get_clone_counts(infe_list, clone="control")
  return(list("refc"=refcounts, "allc"=model_list))
}

# due to bugs in depMixS4 here we extract all valid models for each infected chromosome
# we are attempting to fit a 3 state HMM and if it works we save the model
# @ returns a list of models
get_all_valid_models <- function(res, AllSamples) {
  rcounts = get_infected_counts(res, AllSamples)
  b = rcounts$allc
  models = list()
  pin = 1
  
  print("testing models...")
  for(x in 1:length(b)) {
    set.seed(1)
    test = data.frame(1,(unlist(b[[x]]$wins_sam1)))
    test = test[test[,2]>1,]
    colnames(test) = c("pos", "counts1")
    #test[,2] = log(test[,2]+1)+5
    mod1 <- depmix(counts1~1, nstates = 3, family = poisson(), data=test)
    t = try( fit(mod1, verbose = TRUE) )
    if(class(t)!="try-error") {
      print("yey!")
      set.seed(1)
      models[[pin]] = fit(mod1, verbose = TRUE)
      pin = pin+1
    }
  }
  return(models)
}

# main function to calculate the read count per window across chrosmomes
# @ returns a list of readcount for sample1, sample2 and combined
get_clone_counts <- function(res, steps = STEP, win = WIN, clone="LGZ") {
  # setup clones and samples
  clones = unique(res$clonename)
  i = which(clones==clone)
  test_data = res[res$clonename==clones[i],]
  # window parameters
  keep = c(1)
  # keep only read passing declutter
  test_data = test_data[test_data$keep %in% keep,]
  # get counts across windows for sample1, sample2, and combined
  samples =  unique(test_data$samplename)
  chrs = unique(test_data$lig.chr)
  chrs = chrs[order(as.numeric(gsub("chr", "", chrs)))]
  chr_ranges = lapply(1:length(chrs), function(x) cbind(min(test_data$lig.pos[test_data$lig.chr == chrs[x]]),
                                                        max(test_data$lig.pos[test_data$lig.chr == chrs[x]])))
  win_all = sapply(1:length(chrs), 
                   function(x) countInWindow(
                     test_data$lig.pos[test_data$lig.chr == chrs[x] &
                                         test_data$samplename %in% samples], 
                     steps, win, wintype = "overlapping", 
                     range = unlist(chr_ranges[x]), report = "mid"))
  wins_sam1 = sapply(1:length(chrs), 
                     function(x) countInWindow(
                       test_data$lig.pos[test_data$lig.chr == chrs[x] & 
                                           test_data$samplename %in% samples[1]], 
                       steps, win, wintype = "overlapping", 
                       range = unlist(chr_ranges[x]), report = "mid"))
  wins_sam2 = sapply(1:length(chrs), 
                     function(x) countInWindow(
                       test_data$lig.pos[test_data$lig.chr == chrs[x] & 
                                           test_data$samplename %in% samples[2]], 
                       steps, win, wintype = "overlapping", 
                       range = unlist(chr_ranges[x]), report = "mid"))

  return(list("win_all" = win_all, "wins_sam1" = wins_sam1, "wins_sam2" = wins_sam2))
}

# this is the defencive function to fit HMM to a clone
# we simply try all models and take the first that works
# @ returns a list of input data, state definition and porbabilies
fit_clone <- function(res, clone, models, sample=1) {
  set.seed(1)
  rcounts = get_clone_counts(res, clone=clone)
  res = NULL
  test = NULL
  if(sample==1) {
    test = data.frame(as.numeric(names(unlist (rcounts$wins_sam1[-c(24, 25)])))
                      ,(unlist(rcounts$wins_sam1[-c(24, 25)]))
                      , unlist(lapply(1:length(rcounts$wins_sam1[-c(24, 25)]), 
                                      function(x) rep(x, length(rcounts$wins_sam1[[x]])))))
  } else {
    test = data.frame(as.numeric(names(unlist (rcounts$wins_sam2[-c(24, 25)])))
                      ,(unlist(rcounts$wins_sam2[-c(24, 25)]))
                      , unlist(lapply(1:length(rcounts$wins_sam2[-c(24, 25)]), 
                                      function(x) rep(x, length(rcounts$wins_sam2[[x]])))))
  }
  pin = 1
  colnames(test) = c("pos", "counts1", "chr")
  test = test[test[,2]>1,]
  mod <- depmix(counts1~1, nstates = 3, family = poisson(), data=test)
  hmm <- setpars(mod,getpars(models[[pin]]))
  t = try(viterbi(hmm))
  if(class(t)!="try-error") {
    v = viterbi(hmm)
    res = v
  } else {
    while (t=="try-error") {
      pin = pin+1
      hmm <- setpars(mod,getpars(models[[pin]]))
      t = try(viterbi(hmm))
      print(pin)
    }
    if(class(t)!="try-error") {
      v = viterbi(hmm)
      res = v
    }
  }
  return(list("data"=test, "state"=res))
}

# loop all clones and fit 3 state HMM to each sample
# @return list of clones, samples, input data and hmm result
fit_all_clones <- function(res, AllSamples) {
  require(depmixS4)
  print("Getting all valid models...")
  models = get_all_valid_models(res, AllSamples)
  # loop all clones
  resl = list()
  clones = unique(res$clonename)
  clones = clones[clones!="Jurkat"]
  print("looping across all clones...")
  for(x in 1:length(clones)) {
    fitted_1 = fit_clone(res, clones[x], models, 1)
    fitted_2 = fit_clone(res, clones[x], models, 2)
    resl[[clones[x]]] = list()
    resl[[clones[x]]]$clone1 = fitted_1
    resl[[clones[x]]]$clone2 = fitted_2
    resl[[clones[x]]]$clone_name = clones[x]
    print(x)
  }
  return(resl)
}

# make a per chrosmome plot for a clone
# show sample 1 and sample 2 and hmm state definitions
plot_fitted_clone_per_chromosome <- function(fitted_clone) {
  rchrs = as.numeric(unique(c(fitted_clone$clone1$data$chr, fitted_clone$clone2$data$chr)))
  rchrs = rchrs[order(rchrs)]
  par(mfrow=c(4, ceiling(length(rchrs)/4)))
  for(y in 1:length(rchrs)) {
    ind1 = which(fitted_clone$clone1$data$chr==rchrs[y])
    ind2 = which(fitted_clone$clone2$data$chr==rchrs[y])
    m = c(max(fitted_clone$clone1$data$counts1[ind1]), max(fitted_clone$clone2$data$counts1[ind2]))
    m=m[m>0]
    if(length(ind1)>0) {
      plot(fitted_clone$clone1$data$pos[ind1], fitted_clone$clone1$data$counts1[ind1], col=fitted_clone$clone1$state[ind1,1], ylim=c(-max(m, na.rm=T), max(m, na.rm=T)), pch="+", main=paste("chr", rchrs[y], sep=""), ylab="+- read counts", xlab="chromosome position")
      if(length(ind1)>0) {
        points(fitted_clone$clone2$data$pos[ind2], -fitted_clone$clone2$data$counts1[ind2], col=as.factor(fitted_clone$clone2$state[ind2,1]))
      }
    } else {
      plot(fitted_clone$clone2$data$pos[ind2], -fitted_clone$clone2$data$counts1[ind2], col=as.factor(fitted_clone$clone2$state[ind2,1]), ylim=c(-max(m, na.rm=T), max(m, na.rm=T)), pch="+", main=paste("chr", rchrs[y], sep=""), ylab="+- read counts", xlab="chromosome position")
    }
  }
  plot(0, pch="")
  legend("center", paste("clone:", fitted_clones[[x]]$clone_name))
}

# plot only the intergrated chromosome but by chromosome position and simply index
# NB. genomic windows with 1 or less read were excluded from the input to the hmm
plot_intergrated_chromosome <- function(fitted_clone, AllSamples) {
  info = AllSamples[AllSamples$clonename==fitted_clone$clone_name,]
  chr = unique(info$Chr)
  
  if(is.na(chr)) {
    t = table(fitted_clone$clone1$state[,1], fitted_clone$clone1$data[,3])
    chr = which.max(t[nrow(t),])
  }
  
  chr = gsub("chr", "", chr)
  if(chr=="X") {
    chr = "23"
  }
  par(mfrow=c(2,2))
  chr = as.numeric(chr)
  ind1 = which(fitted_clone$clone1$data$chr==chr)
  ind2 = which(fitted_clone$clone2$data$chr==chr)
  plot(fitted_clone$clone1$data$pos[ind1], fitted_clone$clone1$data$counts1[ind1], col=fitted_clone$clone1$state[ind1,1], ylab="read counts", xlab="chromosome position", main=paste("clone: ", fitted_clones[[x]]$clone_name, " - sample: ", info$SampleID[1], " - chr", chr, sep=""))
  plot(fitted_clone$clone1$data$counts1[ind1], col=fitted_clone$clone1$state[ind1,1], ylab="read counts", xlab="Index", , main=paste("clone: ", fitted_clones[[x]]$clone_name, " - sample: ", info$SampleID[1], " - chr", chr, sep=""))
  plot(fitted_clone$clone2$data$pos[ind2], fitted_clone$clone2$data$counts1[ind2], col=fitted_clone$clone2$state[ind2,1], ylab="read counts", xlab="chromosome position", , main=paste("clone: ", fitted_clones[[x]]$clone_name, " - sample: ", info$SampleID[2], " - chr", chr, sep=""))
  plot(fitted_clone$clone2$data$counts1[ind2], col=fitted_clone$clone2$state[ind2,1], ylab="read counts", xlab="Index", main=paste("clone: ", fitted_clones[[x]]$clone_name, " - sample: ", info$SampleID[2], " - chr", chr, sep=""))
}
