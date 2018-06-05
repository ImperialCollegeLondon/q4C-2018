#@ This code contains sample scripts for analysis of data as described in Melamed, Yaguchi et al., 
#@ "The human leukemia virus HTLV-1 alters the structure and transcription of host chromatin *in cis*" eLife, 2018


correct.lig.site <- function(x, ISchr, ISposition ) {

#@ correct.lig.site - function for correcting the ligation site to add the length of a provirus. 

  # x <- x[x$lig.chr == "HTLV-1" | x$lig.chr == ISchr, ]
  # x <- x[x$shear.chr == "HTLV-1" | x$shear.chr == ISchr, ]
  x <- x[x$shear.chr != "*", ]
  
  ISposition <- as.numeric(ISposition)
  x$lig.pos <- as.numeric(x$lig.pos)
  x$shear.pos <- as.numeric(x$shear.pos)
  
  x$fixligpos <- x$lig.pos
  x$fixligpos[x$lig.chr == ISchr & x$lig.pos < ISposition] <-
    x$lig.pos[x$lig.chr == ISchr & x$lig.pos < ISposition]
  x$fixligpos[x$lig.chr == ISchr & x$lig.pos >= ISposition] <-
    x$lig.pos[x$lig.chr == ISchr & x$lig.pos >= ISposition] + 9000
  x$fixligpos[x$lig.chr == "HTLV-1"] <-
    x$lig.pos[x$lig.chr == "HTLV-1"] + ISposition
  
  x$fixshearpos <- x$shear.pos
  x$fixshearpos[x$shear.chr == ISchr & x$shear.pos < ISposition] <-
    x$shear.pos[x$shear.chr == ISchr & x$shear.pos < ISposition]
  x$fixshearpos[x$shear.chr == ISchr & x$shear.pos >= ISposition] <-
    x$shear.pos[x$shear.chr == ISchr & x$shear.pos >= ISposition] + 9000
  x$fixshearpos[x$shear.chr == "HTLV-1"] <-
    x$shear.pos[x$shear.chr == "HTLV-1"] + ISposition
  
  x$fixchr <- x$lig.chr
  x$fixchr[x$lig.chr == ISchr | x$lig.chr == "HTLV-1"] <- ISchr
  
  x$class <- NA
  x$class[x$lig.chr == ISchr | x$lig.chr == "HTLV-1"] <-
    assignLigSiteClass(x[x$lig.chr == ISchr | x$lig.chr == "HTLV-1", ],
                       ligChrCol = "fixchr",
                       shearChrCol = "fixchr",
                       ligSiteCol = "fixligpos",
                       shearSiteCol = "fixshearpos")
  x$class[x$lig.chr != ISchr & x$lig.chr != "HTLV-1"] <-
    assignLigSiteClass(x[x$lig.chr != ISchr & x$lig.chr != "HTLV-1", ])
 
  x <- x[x$class == "convergent", ]

  return(x)
}

########

assignLigSiteClass <- function(a, ligSiteCol = NULL, ligChrCol = NULL, 
                               ligStrandCol = NULL, shearSiteCol = NULL, 
                               shearChrCol = NULL, shearStrandCol = NULL){

#@ function for classification of ligation events based on the relative orientation 
#@ read1 and read2 - noise is drastically reduced when limiting to convergent orientation

  if(is.null(ligSiteCol)){ligSiteCol <- "lig.pos"}
  if(is.null(ligChrCol)){ligChrCol <- "lig.chr"}
  if(is.null(ligStrandCol)){ligStrandCol <- "lig.strand"}
  if(is.null(shearSiteCol)){shearSiteCol <- "shear.pos"}
  if(is.null(shearChrCol)){shearChrCol <- "shear.chr"}
  if(is.null(shearStrandCol)){shearStrandCol <- "shear.strand"}

  a$class <- rep( NA, dim( a )[1] )
  
  a$class[a[, shearChrCol] == "*" ] <- "ambiguous"
  
  a$class[a[, shearChrCol] != "*" & a[, ligChrCol] != a[, shearChrCol] ] <- "alt chr"
  
  a$class[a[, ligChrCol] == a[, shearChrCol] &
            a[, ligStrandCol] == a[, shearStrandCol] ] <- "tandem"
   
  a$class[a[, ligChrCol] == a[, shearChrCol] &
            a[, ligStrandCol] != a[, shearStrandCol] &
            a[, ligStrandCol] == 0 &
            a[, shearSiteCol] >= a[, ligSiteCol] ] <- "convergent"
  
  a$class[a[, ligChrCol] == a[, shearChrCol] &
            a[, ligStrandCol] != a[, shearStrandCol] &
            a[, ligStrandCol] == 16 &
            a[, shearSiteCol] <= a[, ligSiteCol] ] <- "convergent"
  
  a$class[a[, ligChrCol] == a[, shearChrCol] &
            a[, ligStrandCol] != a[, shearStrandCol] &
            a[, ligStrandCol] == 0 &
            a[, shearSiteCol] < a[, ligSiteCol] ] <- "divergent"
  
  a$class[a[, ligChrCol] == a[, shearChrCol] &
            a[, ligStrandCol] != a[, shearStrandCol] &
            a[, ligStrandCol] == 16 &
            a[, shearSiteCol] > a[, ligSiteCol] ] <- "divergent"

  return(a$class)
}


#########
## rounddown
## input - number (vector) and rounding factor
## output - number rounded
rounddown <- function(x = c(), by = 10)
{
  return(x-x%%by)
}

########
countInWindow <- function( dataset = NULL, step = NULL, win = NULL, wintype = NULL, 
                           range = NULL, report = c( "min", "max", "mid" ) ){
#@ bin ligation event data to bin by overlapping or discrete windows
#@ output used for plotting and for identification of peaks

  require(GenomicRanges)
  
  if( is.null( dataset ) ){ stop( "Define dataset to scan\n" ) }
  if( is.null( step ) ){ stop( "Define step length\n" ) }
  if( is.null( win ) ){ stop( "Define window width\n" ) }
  if( is.null( wintype ) ){ stop( "Choose window type. available options: \"discrete\", \"overlapping\"\n" ) }
  
  if( is.data.frame( dataset ) ){
    if( is.null( dataset$Position ) ){
      stop( "dataset must either be a dataframe with a 'Position' column or a vector of positions" )
    }
    positions <- dataset$Position
  }else if( is.vector( dataset ) ){
    positions <- dataset
  }else{
    stop( "dataset must either be a dataframe with a 'Position' column or a vector of positions" )
  }

  if( is.null( range ) )
  {
    winmin <- min( positions )
    winmax <- max( positions )
  }
  else
  {
    winmin <- range[1]
    winmax <- range[2]
  }
  
  report <- report[1]
  
  # discrete windows
  if( wintype == "discrete" ){
    
    winlow <- seq( winmin, winmax, win )
    
    winhigh <- winlow+win-1
    
  }else if( wintype == "overlapping" ){
    
    winlow <- seq( winmin, 
                   winmax, step )
    
    winhigh <- winlow+win
    
  }else{
    stop( "Choose window type. available options: \"discrete\", \"overlapping\"" )
  }
  
  window.rl <- IRanges( winlow, winhigh )
  
  sites.rl <- IRanges( positions, positions )
  
  overlaps <- countOverlaps( window.rl, sites.rl )
  
  names(overlaps) <- switch( report, min = winlow, max = winhigh, mid = (winlow+winhigh)/2 )
  
  return ( overlaps )
  
}

######




