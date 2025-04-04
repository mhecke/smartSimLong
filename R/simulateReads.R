
### helper functions

#add random error to fragments
#code from Polyester package (https://github.com/alyssafrazee/polyester-release/blob/master/R/add_error.R)
add_error <- function(tFrags, error_rate = 0.005){
  adj_error = error_rate*4/3
  #^so you don't have to choose *another* nucleotide for an error: just
  # choose *a* nucleotide.

  allSeq = unlist(tFrags)
  insertLocs = sample(c(TRUE,FALSE), size = length(allSeq),
                      replace=TRUE, prob = c(adj_error, 1-adj_error))

  newletters = DNAString(
    paste(sample(c("A", "C", "G", "T"), sum(insertLocs),
                 replace=TRUE), collapse="") )
  allSeq = replaceLetterAt(allSeq, insertLocs, newletters)

  eFrags = DNAStringSet(allSeq,
                        start=c(1, (cumsum(width(tFrags))+1)[-length(tFrags)]),
                        width=width(tFrags))

  names(eFrags) = names(tFrags)
  return(eFrags)
}

generateRandomSeq <- function(len){
  paste0(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
}
###

#' @title simulateReads()
#'
#' @description simulate Smart-seq3 reads for selected transcripts
#'
#' @param transcriptExpression expression for each transcript and cell
#' @param transcript_seq sequence for each transcript
#' @param dataChar data characteristics
#' @param out output directory to write reads to
#' @param barcodeLen length of barcode [default = 8]
#' @param UMIlen length of UMI [default = 10]
#' @param readLen length of read [default = 150]
#' @param tag tag to indicate UMI containing read [default = "ATTGCGCAATG"]
#' @param PAF PCR amplification factor [default = 2]
#' @param errorRate errorRate per bp [default = 0.005]
#' @param long also generate long reads [default = long]
#' @param longErrorRate error rate for long reads [default = 0.05]
#' @param seed for reporducibility [default = 1234]
#'
#' @importFrom Biostrings subseq width replaceLetterAt BStringSet
#' @importFrom ShortRead ShortReadQ writeFastq
#' @export
simulateReads <- function(transcriptExpression, transcript_seq, dataChar, out, barcodeLen = 8,
                          UMIlen = 10, readLen = 150, tag = "ATTGCGCAATG", PAF = 2, errorRate = 0.005,
                          long = FALSE, longErrorRate = 0.05, seed = 1234){

  print("simulating reads")
  dir.create(out, showWarnings = FALSE)
  set.seed(seed)

  #get number of PCR cycles from data (number of reads per UMI)
  nPCR <- sampleNumReadsPerUMI(dataChar, n = 1) * PAF

  I1_all <- NULL; I2_all <- NULL
  R1_all <- NULL; R2_all <- NULL
  cell2barcode <- NULL
  longBarcodes <- NULL

  #remove existing files if they exist
  if (file.exists(paste0(out, "/unassigned_R1.fq.gz"))) file.remove(paste0(out, "/unassigned_R1.fq.gz"))
  if (file.exists(paste0(out, "/unassigned_R2.fq.gz"))) file.remove(paste0(out, "/unassigned_R2.fq.gz"))
  if (file.exists(paste0(out, "/unassigned_I1.fq.gz"))) file.remove(paste0(out, "/unassigned_I1.fq.gz"))
  if (file.exists(paste0(out, "/unassigned_I2.fq.gz"))) file.remove(paste0(out, "/unassigned_I2.fq.gz"))
  
  if (long){
    if (file.exists(paste0(out, "/long.fq.gz"))) file.remove(paste0(out, "/long.fq.gz"))
  }

  #generate reads separately per cell
  for (cell in unique(transcriptExpression$cell)){
    print(cell)

    #generate random barcode
    barcode_1 <- generateRandomSeq(barcodeLen)
    barcode_2 <- generateRandomSeq(barcodeLen)
    barcode <- paste0(barcode_1, barcode_2)
    cell2barcode <- rbind(cell2barcode, c(cell, barcode))
    longBarcodes <- rbind(longBarcodes, barcode_2)

    cellExpression <- transcriptExpression[transcriptExpression$cell == cell,]

    #generate random UMIs
    randomUMIs <- replicate(sum(cellExpression$transcriptExpression),
                            generateRandomSeq(UMIlen))

    #generate tagged + UMI transcripts
    thisTranscripts <- rep(transcript_seq[cellExpression$transcript_id], cellExpression$transcriptExpression)
    taggedTranscripts <- DNAStringSet(paste0(tag, randomUMIs, "GGG", thisTranscripts))
    names(taggedTranscripts) <- names(thisTranscripts)

    #remove transcripts that have length < readLen
    taggedTranscripts = taggedTranscripts[width(taggedTranscripts) >= readLen]

    #PCR amplification
    taggedTranscripts <- rep(taggedTranscripts, nPCR)

    #tagmentation

    #UMI containing reads
    UMIfragLen <- sampleUMIfragLen(dataChar, n = length(taggedTranscripts), minLen = readLen)
    UMIfragLen <- pmax(pmin(UMIfragLen, width(taggedTranscripts)), readLen)
    UMIfrags <- subseq(taggedTranscripts, start = 1, end = UMIfragLen)
    remainingSeq <- subseq(taggedTranscripts, start = UMIfragLen + 1, end = width(taggedTranscripts))

    #remove remaining sequences with length < readlen
    remainingSeq <- remainingSeq[width(remainingSeq) >= readLen]

    #number of internal reads
    remainingUMIcounts <- table(names(remainingSeq))
    nInternal <- round(sampleNonUMIvsUMI(dataChar, n = length(remainingUMIcounts)) * remainingUMIcounts)
    #get sequences to take internal reads from
    remainingSeq <- split(remainingSeq, names(remainingSeq))
    remainingSeq <- mapply(function(seqs, n) seqs[sample(length(seqs), min(n, length(seqs)), replace = TRUE)],
                           remainingSeq, nInternal[names(remainingSeq)], SIMPLIFY = FALSE)
    #convert back to DNAStringSet
    names(remainingSeq) <- NULL
    remainingSeq <- do.call(c, remainingSeq)

    internalFragLen <- sampleNonUMIfragLen(dataChar, n = length(remainingSeq), minLen = readLen)
    internalFragLen <- pmax(pmin(internalFragLen, width(remainingSeq)), readLen)
    validStarts <- width(remainingSeq) - internalFragLen + 1
    randomStarts <- sapply(validStarts, function(x) { sample(1:x, size = 1) })
    nonUMIfrags <- subseq(remainingSeq, start = randomStarts, width = internalFragLen)

    #combine UMI and nonUMI frags
    frags <- c(UMIfrags, nonUMIfrags)

    #add random read errors
    frags <- add_error(frags, errorRate)

    #subsample
    frags <- sample(frags, round(length(frags)/PAF), replace = FALSE)

    #get reads
    R1 <- subseq(frags, start = 1, end = readLen)
    R2 <- reverseComplement(subseq(frags, start = width(frags)-readLen+1, end = width(frags)))

    names(R1) <-  paste0(barcode, "_", names(R1), ":read", seq(1, length(R1)))
    names(R2) <-  paste0(barcode, "_", names(R2), ":read", seq(1, length(R2)))

    I1 <- DNAStringSet(rep(barcode_1, length(R1)))
    I2 <- DNAStringSet(rep(barcode_2, length(R2)))
    names(I1) <- paste0(names(R1), " ", "1:N:0", barcode_1, "+", barcode_2)
    names(I2) <- paste0(names(R2), " ", "2:N:0", barcode_1, "+", barcode_2)

    # Write to FASTQ
    fastq_Q1 <- BStringSet(sampleReadQual(dataChar, readLen, length(R1)))
    fastq_Q2 <- BStringSet(sampleReadQual(dataChar, readLen, length(R2)))
    fastq_R1 <- ShortReadQ(sread = R1, quality = fastq_Q1, id = BStringSet(names(R1)))
    fastq_R2 <- ShortReadQ(sread = R2, quality = fastq_Q2, id = BStringSet(names(R2)))

    writeFastq(fastq_R1, paste0(out, "/unassigned_R1.fq.gz"), compress = TRUE, mode = "a")
    writeFastq(fastq_R2, paste0(out, "/unassigned_R2.fq.gz"), compress = TRUE, mode = "a")

    fastq_Q1 <- BStringSet(sampleBCQual(dataChar, barcodeLen, length(I1)))
    fastq_Q2 <- BStringSet(sampleBCQual(dataChar, barcodeLen, length(I1)))
    fastq_I1 <- ShortReadQ(sread = I1, quality = fastq_Q1, id = BStringSet(names(I1)))
    fastq_I2 <- ShortReadQ(sread = I2, quality = fastq_Q2, id = BStringSet(names(I2)))

    writeFastq(fastq_I1, paste0(out, "/unassigned_I1.fq.gz"), compress = TRUE, mode = "a")
    writeFastq(fastq_I2, paste0(out, "/unassigned_I2.fq.gz"), compress = TRUE, mode = "a")

    #generate long reads
    if(long){
      fwd_primer_libr = "AATGATACGGCGACCACCGAGATCTACAC"
      fwd_primer = "TCGTCGGCAGCGTCAGATGTGTATA"
      TSO = "AGAGACAG"
      spacer = "AT"
      
      longExpr <- round(cellExpression$transcriptExpression * sampleLongVsShort(dataChar, length(cellExpression$transcriptExpression)))
      
      longUMIs <- replicate(sum(longExpr),
                           generateRandomSeq(UMIlen))
      longTranscripts <- rep(transcript_seq[cellExpression$transcript_id], longExpr)
      taggedLongTranscripts <- DNAStringSet(paste0(fwd_primer_libr,barcode_2, fwd_primer, TSO, tag, longUMIs, spacer, "GGG", longTranscripts))
      names(taggedLongTranscripts) <- names(longTranscripts)
      
      #number of reads per long UMI
      nReadsLong = round(sampleNumLongReads(dataChar, length(taggedLongTranscripts)))
      taggedLongTranscripts <- rep(taggedLongTranscripts, nReadsLong)
      
      #long read lengths
      longReadLen <- round(sampleLongReadLen(dataChar, length(taggedLongTranscripts)) * width(taggedLongTranscripts))
      longReads <- subseq(taggedLongTranscripts, start = 1, end = longReadLen)
      names(longReads) <-  paste0(barcode, "_", names(longReads), ":read", seq(1, length(longReads)))
      
      # Write to FASTQ
      fastq_Q <- BStringSet(sampleLongReadQual(dataChar, longReadLen))
      fastq_R <- ShortReadQ(sread = longReads, quality = fastq_Q, id = BStringSet(names(longReads)))
      
      writeFastq(fastq_R, paste0(out, "/long.fq.gz"), compress = TRUE, mode = "a")
      
    }
    
    
  }

  cell2barcode <- as.data.frame(cell2barcode)
  colnames(cell2barcode) <- c("cell", "barcode")

  #write concatenated barcodes to file
  write.table(cell2barcode$barcode, paste0(out, "/barcodes.txt"),
              row.names = F, col.names = F, quote = F)
  write.table(cell2barcode,
              paste0(out, "/cell2barcode.tsv"),
              sep = "\t", row.names = F, quote = F)
  if (long){
    write.table(longBarcodes, paste0(out, "/longBarcodes.txt"), 
                row.names = F, col.names = F, quote = F)
    
    write.table(cbind(cell2barcode$barcode, longBarcodes), paste0(out, "/short2longBarcodes.txt"), 
                row.names = F, col.names = F, quote = F)
  }
}

