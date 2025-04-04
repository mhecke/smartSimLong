#' @title DataCharacteristics
#'
#' @description data characteristics based on existing Smart-seq3 data
#'
#' @param inDir directory containing output from characterizeData
#' @return object of class "DataCharacteristics"
#'
#' @importFrom dplyr %>% group_by n_distinct summarise
#' @export
DataCharacteristics <- function(inDir){

  longReads <- FALSE
  
  #read data
  UMIcounts <- read.table(paste0(inDir, "/", "UMIcounts.tsv"), sep = "\t", header = TRUE)
  fragLen <- read.table(paste0(inDir, "/", "fragmentLength.tsv"), sep = "\t", header= TRUE)
  fragLen_UMI <- read.table(paste0(inDir, "/", "fragmentLength_UMI.tsv"), sep = "\t", header = TRUE)
  fragLen_nonUMI <- read.table(paste0(inDir, "/", "fragmentLength_nonUMI.tsv"), sep = "\t", header = TRUE)
  readQual <- read.table(paste0(inDir, "/", "readQual.tsv"), sep = "\t", header = TRUE, comment.char="", quote = NULL)
  bcQual <- read.table(paste0(inDir, "/", "readQual.tsv"), sep = "\t", header = TRUE, comment.char="", quote = NULL)
  
  if (file.exists(paste0(inDir, "/", "long_fragmentLength.tsv"))){
    longReads <- TRUE
    longFragLen <- read.table(paste0(inDir, "/", "long_fragmentLength.tsv"), sep = "\t", header= TRUE)
    longUMIcounts <- read.table(paste0(inDir, "/", "long_UMIcounts.tsv"), sep = "\t", header = TRUE)
    longReadQual <- read.table(paste0(inDir, "/", "long_readQual.tsv"), sep = "\t", header = TRUE, comment.char="", quote = NULL)
  }

  # number of unique UMIs per gene
  UMIsPerGene <- UMIcounts[UMIcounts$UMI != "noUMI",] %>%
    group_by(cell, gene) %>%
    summarize(numUMI = n_distinct(UMI))
  q01 <- quantile(UMIsPerGene$numUMI, 0.01)
  q99 <- quantile(UMIsPerGene$numUMI, 0.99)
  UMIsPerGene_hist <- hist(UMIsPerGene$numUMI[UMIsPerGene$numUMI >= q01 & UMIsPerGene$numUMI <= q99],
                           breaks = seq(1, max(UMIsPerGene$numUMI)),
                           plot = FALSE)

  #number of reads per UMI
  readsPerUMI <- UMIcounts[UMIcounts$UMI != "noUMI",]
  q01 <- quantile(readsPerUMI$numReads, 0.01)
  q99 <- quantile(readsPerUMI$numReads, 0.99)
  readsPerUMI_hist <- hist(readsPerUMI$numReads[readsPerUMI$numReads >= q01 & readsPerUMI$numReads <= q99],
                           breaks = seq(1, max(readsPerUMI$numReads)),
                           plot = FALSE)

  #fragment length for UMI-containing reads
  fragLen_UMI <- fragLen_UMI[order(fragLen_UMI$val), ]
  fragLen_UMI$cumSum <- cumsum(fragLen_UMI$count)
  tot <- fragLen_UMI$cumSum[nrow(fragLen_UMI)]
  fragLen_UMI <- fragLen_UMI[fragLen_UMI$cumSum >= 0.01*tot, ] #q01
  fragLen_UMI <- fragLen_UMI[fragLen_UMI$cumSum <= 0.99*tot, ] #q99
  fragLen_UMI_hist <- list("mids" = fragLen_UMI$val, "counts" = fragLen_UMI$count)

  #fragment length for non UMI-containing reads
  fragLen_nonUMI <- fragLen_nonUMI[order(fragLen_nonUMI$val), ]
  fragLen_nonUMI$cumSum <- cumsum(fragLen_nonUMI$count)
  tot <- fragLen_nonUMI$cumSum[nrow(fragLen_nonUMI)]
  fragLen_nonUMI <- fragLen_nonUMI[fragLen_nonUMI$cumSum >= 0.01*tot, ] #q01
  fragLen_nonUMI <- fragLen_nonUMI[fragLen_nonUMI$cumSum <= 0.99*tot, ] #q99
  fragLen_nonUMI_hist <- list("mids" = fragLen_nonUMI$val, "counts" = fragLen_nonUMI$count)


  #distribution of number of nonUMI reads vs. UMI reads
  numNonUMIvsUMI <- UMIcounts %>%
    group_by(cell, gene) %>%
    summarize(numNonUMIvsUMI = sum(numReads[UMI == "noUMI"]) / sum(numReads[UMI != "noUMI"]))
  numNonUMIvsUMI <- numNonUMIvsUMI$numNonUMIvsUMI[! is.infinite(numNonUMIvsUMI$numNonUMIvsUMI)]
  q01 <- quantile(numNonUMIvsUMI, 0.01)
  q99 <- quantile(numNonUMIvsUMI, 0.99)
  nonUMIvsUMI_hist <- hist(numNonUMIvsUMI[numNonUMIvsUMI >= q01 & numNonUMIvsUMI <= q99],
                           breaks = seq(0, max(numNonUMIvsUMI), by=0.1),
                           plot = FALSE)

  #quality
  readQual$pos <- readQual$pos + 1 #R uses 1-based indexing
  bcQual$pos <- bcQual$pos + 1 #R uses 1-based indexing

  if (longReads){
    UMIcounts_merged <- merge(UMIcounts, longUMIcounts,by=c("cell", "gene", "UMI"), suffixes = c("_short", "_long"), 
                              all = TRUE)
    
    UMIs_longVsShort <- UMIcounts_merged %>%
      group_by(cell, gene) %>%
      summarise(
        unique_long = n_distinct(numReads_long, na.rm = TRUE),  # Count unique values in col_long
        unique_short = n_distinct(numReads_short, na.rm = TRUE),  # Count unique values in col_short
        ratio = unique_long / unique_short,  # Compute ratio
        .groups = "drop"
      )
    UMIs_longVsShort <- UMIs_longVsShort$ratio
    UMIs_longVsShort <- UMIs_longVsShort[is.finite(UMIs_longVsShort)]
    q01 <- quantile(UMIs_longVsShort, 0.01)
    q99 <- quantile(UMIs_longVsShort, 0.99)
    UMIs_longVsShort_hist <- hist(UMIs_longVsShort[UMIs_longVsShort >= q01 & UMIs_longVsShort <= q99],
                             breaks = seq(0, max(UMIs_longVsShort), by=0.01),
                             plot = FALSE)
    
    UMIs_LongNew <- UMIcounts_merged %>%
      group_by(cell, gene) %>%
      summarise(newUMIs = sum(is.na(numReads_short))/sum(!is.na(numReads_long)), .groups = "drop")
    UMIs_LongNew <- UMIs_LongNew[! is.na(UMIs_LongNew$newUMIs),]
    q01 <- quantile(UMIs_LongNew$newUMIs, 0.01)
    q99 <- quantile(UMIs_LongNew$newUMIs, 0.99)
    UMIs_longNew_hist <- hist(UMIs_LongNew$newUMIs[UMIs_LongNew$newUMIs >= q01 & UMIs_LongNew$newUMIs <= q99],
                             breaks = seq(0, max(UMIs_LongNew$newUMIs), by = 0.01),
                             plot = FALSE)
    
    readsPerUMI_long <- longUMIcounts$numReads
    q01 <- quantile(readsPerUMI_long, 0.01)
    q99 <- quantile(readsPerUMI_long, 0.99)
    readsPerUMI_long_hist <- hist(readsPerUMI_long[readsPerUMI_long >= q01 & readsPerUMI_long <= q99],
                         breaks = seq(1, max(readsPerUMI_long)),
                         plot = FALSE)

    #get percent of genelen covered
    longFragLen_hist <- list("mids" = longFragLen$val, "counts" = longFragLen$count)
    
    longReadQual$pos <- longReadQual$pos + 1 #R uses 1-based indexing
  } else{
    UMIs_longVsShort_hist = NULL
    UMIs_longNew_hist = NULL
    readsPerUMI_long_hist = NULL
    longFragLen_hist = NULL
    longReadQual = NULL
  }
  
  #create object
  obj <- list("fragLen_UMI" = fragLen_UMI_hist,
              "fragLen_nonUMI" = fragLen_nonUMI_hist,
              "UMIsPerGene" = UMIsPerGene_hist,
              "readsPerUMI" = readsPerUMI_hist,
              "nonUMIvsUMI" = nonUMIvsUMI_hist,
              "readQual" = readQual,
              "bcQual" = bcQual, 
              "longVsShortUMI"= UMIs_longVsShort_hist, 
              "newLongUMI" = UMIs_longNew_hist,
              "readsPerUMI_long" = readsPerUMI_long_hist,
              "fragLen_long" = longFragLen_hist,
              "readQual_long"=longReadQual
              )
  class(obj) <- "DataCharacteristics"
  return(obj)

}

#' @title sampleUMIcounts
#'
#' @description sample n elements from the distribution of UMIcounts
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#' @param minCount minimum number of counts
#'
#' @export
sampleUMIcounts <- function(dataChar, n, minCount = 0){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$UMIsPerGene$mids[dataChar$UMIsPerGene$mids > minCount] - 0.5
  counts <- dataChar$UMIsPerGene$counts[dataChar$UMIsPerGene$mids > minCount]
  
  return(sample(midpoints, size = n, replace = TRUE, prob = counts))
}

#' @title sampleNumReadsPerUMI
#'
#' @description sample n elements from the distribution of number of reads per UMI
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#'
#' @export
sampleNumReadsPerUMI <- function(dataChar, n){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$readsPerUMI$mids - 0.5
  counts <- dataChar$readsPerUMI$counts

  return(sample(midpoints, size = n, replace = TRUE, prob = counts))
}

#' @title sampleUMIfragLen
#'
#' @description sample n elements from the distribution of fragment length for reads containing UMIs
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#' @param minLen minumum fragment length
#'
#' @export
sampleUMIfragLen <- function(dataChar, n, minLen){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$fragLen_UMI$mids
  counts <- dataChar$fragLen_UMI$counts

  return(sample(midpoints[midpoints >= minLen], size = n, replace = TRUE, prob = counts[midpoints >= minLen]))

}

#' @title sampleNonUMIfragLen
#'
#' @description sample n elements from the distribution of fragment length for reads not containing UMIs
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#' @param minLen minumum fragment length
#'
#' @export
sampleNonUMIfragLen <- function(dataChar, n, minLen){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$fragLen_nonUMI$mids
  counts <- dataChar$fragLen_nonUMI$counts

  return(sample(midpoints[midpoints >= minLen], size = n, replace = TRUE, prob = counts[midpoints >= minLen]))

}

#' @title sampleNonUMIvsUMI
#'
#' @description sample n elements from the distribution of number of nonUMI reads vs. UMI reads
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#'
#' @export
sampleNonUMIvsUMI <- function(dataChar, n){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$nonUMIvsUMI$mids
  counts <- dataChar$nonUMIvsUMI$counts

  return(sample(midpoints, size = n, replace = TRUE, prob = counts))

}

#' @title sampleReadQual
#'
#' @description sample read quality for read of length rLen
#'
#' @param dataChar DataCharacteristics object
#' @param rLen read length
#' @param n number of read qualities to sample
#'
#' @export
sampleReadQual <- function(dataChar, rLen, n){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(rLen) | ! rLen%%1 == 0 | ! rLen>0) stop("value pos = ", rLen, " is not valid")

  maxPos <- max(dataChar$readQual$pos)

  if (rLen > maxPos) warning("original data did not contain reads with length ", rLen)

  qualStr <- matrix("", nrow = rLen, ncol = n)
  for (pos in 1:rLen){

    if (pos <= maxPos){
      midpoints <- dataChar$readQual$qual[dataChar$readQual$pos == pos]
      counts <- dataChar$readQual$count[dataChar$readQual$pos == pos]
    }else{
      midpoints <- dataChar$readQual$qual
      counts <- dataChar$readQual$count
    }
    qualStr[pos,] = sample(midpoints, size = n, prob = counts, replace = TRUE)
  }

  return(apply(qualStr, 2, paste0, collapse = ""))

}


#' @title sampleBCQual
#'
#' @description sample barcode quality for barcode of length rLen
#'
#' @param dataChar DataCharacteristics object
#' @param rLen length of barcode
#' @param n number of barcode qualities to sample
#'
#' @export
sampleBCQual <- function(dataChar, rLen, n){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(rLen) | ! rLen%%1 == 0 | ! rLen>0) stop("value pos = ", rLen, " is not valid")

  maxPos <- max(dataChar$bcQual$pos)

  if (rLen > maxPos) warning("original data did not contain barcodes with length ", rLen)

  qualStr <- matrix("", nrow = rLen, ncol = n)
  for (pos in 1:rLen){

    if (pos <= maxPos){
      midpoints <- dataChar$bcQual[dataChar$bcQual$pos == pos, "qual"]
      counts <- dataChar$bcQual[dataChar$bcQual$pos == pos, "count"]
    }else{
      midpoints <- dataChar$bcQual$qual
      counts <- dataChar$bcQual$count
    }
    qualStr[pos, ] = sample(midpoints, size = n, prob = counts, replace = TRUE)
  }
  return(apply(qualStr, 2, paste0, collapse = ""))

}

#' @title sampleLongVsShort
#'
#' @description sample n elements from the distribution of number of long UMIs vs. short UMIs
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#'
#' @export
sampleLongVsShort <- function(dataChar, n){
  
  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")
  
  midpoints <- dataChar$longVsShortUMI$mids
  counts <- dataChar$longVsShortUMI$counts
  
  return(sample(midpoints, size = n, replace = TRUE, prob = counts))
  
}

#' @title sampleNumLongReads
#'
#' @description sample n elements from the distribution of number of long reads per UMI
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#'
#' @export
sampleNumLongReads <- function(dataChar, n){
  
  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")
  
  midpoints <- dataChar$readsPerUMI_long$mids
  counts <- dataChar$readsPerUMI_long$counts
  
  return(sample(midpoints, size = n, replace = TRUE, prob = counts))
  
}

#' @title sampleLongReadLen
#'
#' @description sample n elements from the distribution of long read length
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#'
#' @export
sampleLongReadLen <- function(dataChar, n){
  
  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")
  
  midpoints <- dataChar$fragLen_long$mids
  counts <- dataChar$fragLen_long$counts
  
  return(sample(midpoints, size = n, replace = TRUE, prob = counts))
  
}

#' @title sampleLongReadQual
#'
#' @description sample long read quality for reads with length rLen
#'
#' @param dataChar DataCharacteristics object
#' @param rLen read length vector
#'
#' @export
sampleLongReadQual <- function(dataChar, rLen){
  
  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  
  maxPos <- max(dataChar$readQual_long$pos)
  
  if(any(rLen > maxPos)) warning("original data did not contain reads with length > ", maxPos)
  
  qualStr <- matrix("", nrow = max(rLen), ncol = length(rLen))
  for (pos in 1:max(rLen)){
    if (pos <= maxPos){
      midpoints <- dataChar$readQual_long$qual[dataChar$readQual_long$pos == pos]
      counts <- dataChar$readQual_long$count[dataChar$readQual_long$pos == pos]
    }else{
      midpoints <- dataChar$readQual_long$qual
      counts <- dataChar$readQual_long$count
    }
    qualStr[pos,] = sample(midpoints, size = length(rLen), prob = counts, replace = TRUE)
  }
  
  allQualStr = apply(qualStr, 2, paste0, collapse = "")
  
  return(mapply(substr, allQualStr, 1 , rLen))
  
}
