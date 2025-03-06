
### helper functions

#add de novo exon for single gene
addNewExonSingleGene <- function(existingTranscripts, ref_seq, canonical){

  #all transcripts should be from the same gene
  if (length(unique(existingTranscripts$gene_id)) != 1)
    stop("make sure all transcripts in existingTranscripts are from the same gene")
  if (nrow(existingTranscripts[existingTranscripts$type == "gene",]) != 1)
    stop("existingTranscripts does not contain exactly 1 gene row")

  strand = existingTranscripts$strand[1]
  seqnames = existingTranscripts$seqnames[1]
  gene_id = existingTranscripts$gene_id[1]
  gene_type = existingTranscripts$gene_type[1]
  gene_name = existingTranscripts$gene_name[1]

  #get exons
  exons <- existingTranscripts[existingTranscripts$type == "exon",]
  exons$exon_number <- as.numeric(exons$exon_number)
  minExonLen <- min(exons$width)
  maxExonLen <- max(exons$width)

  #get sequence for whole gene
  geneStart <- existingTranscripts[existingTranscripts$type == "gene", "start"]
  geneEnd <- existingTranscripts[existingTranscripts$type == "gene", "end"]
  gene_seq <- subseq(ref_seq, start = geneStart, end = geneEnd)

  if (strand == "+"){
    leftBorder = "AG" #acceptor
    rightBorder = "GT" #donor
  }else{
    leftBorder = "AC"
    rightBorder = "CT"
  }

  if (canonical){
    #get all possible start and end positions for novel exon (canonical splice junctions)
    allStartPos <- unlist(end(vmatchPattern(leftBorder, gene_seq)) + 1) + geneStart - 1
    allEndPos <- unlist(start(vmatchPattern(rightBorder, gene_seq)) - 1) + geneStart - 1
  }else{
    allStartPos <- geneStart:geneEnd
    allEndPos <- geneStart:geneEnd
  }

  #only select internal novel exons
  allStartPos <- allStartPos[allStartPos > min(exons$end)]
  allEndPos <- allEndPos[allEndPos < max(exons$start)]

  #loop over all start positions (in random order)
  allStartPos <- allStartPos[sample(length(allStartPos))] #random shuffle
  for (startPos in allStartPos){

    #select an existing exon that ends before startPos, such that it forms a canonical junction with the new exon
    prevExons <- exons[exons$end < startPos,]
    if (canonical)
      prevExons <- prevExons[unlist(mapply(function(p) as.character(subseq(ref_seq, p+1, p+2)), prevExons$end)) == rightBorder, ]
    if (nrow(prevExons) == 0)
      next
    prevExon <- prevExons[sample(nrow(prevExons), 1), ] #select random exon

    #get all possible end positions that fit the start position
    posEndPos <- allEndPos[allEndPos > (startPos + minExonLen) & allEndPos < (startPos + maxExonLen)]
    posEndPos <- posEndPos[sample(length(posEndPos))]

    for (endPos in posEndPos){
      #check if exon does not yet exist
      if (any(exons$start == startPos | exons$end == endPos))
        next #only add exons that do not have same exon boundary as existing exon

      #select an existing exon that starts after endPos, such that it forms a canonical junction with the new exon
      nextExons <- exons[exons$start > endPos,]
      if(canonical)
        nextExons <- nextExons[unlist(mapply(function(p) as.character(subseq(ref_seq, p-2, p-1)), nextExons$start)) == leftBorder, ]
      if (nrow(nextExons) == 0)
        next
      nextExon <- nextExons[sample(nrow(nextExons), 1), ] #select random exon

      #make transcript using this exon
      if(strand == "+"){
        firstExons <- exons[exons$transcript_id == prevExon$transcript_id & exons$exon_number <= prevExon$exon_number,]
        lastExons <- exons[exons$transcript_id == nextExon$transcript_id & exons$exon_number >= nextExon$exon_number,]
      }else{
        firstExons <- exons[exons$transcript_id == nextExon$transcript_id & exons$exon_number <= nextExon$exon_number,]
        lastExons <- exons[exons$transcript_id == prevExon$transcript_id & exons$exon_number >= prevExon$exon_number,]
      }

      newExon <- data.frame(matrix(NA, nrow = 1, ncol = ncol(existingTranscripts)))
      colnames(newExon) <- colnames(existingTranscripts)

      newExon$start <- startPos
      newExon$end <- endPos
      newExon$width <- endPos - startPos + 1
      newExon$type <- "exon"
      newExon$exon_id <-  paste0(gene_id, "novelExon")

      newTranscript <- rbind(rbind(firstExons, newExon), lastExons)
      newTranscript$exon_number <- 1:nrow(newTranscript)

      #add transcript header
      transcriptRow <- data.frame(matrix(NA, nrow = 1, ncol = ncol(newTranscript)))
      colnames(transcriptRow) <- colnames(newTranscript)
      transcriptRow$start = min(newTranscript$start)
      transcriptRow$end = max(newTranscript$end)
      transcriptRow$width = transcriptRow$end - transcriptRow$start + 1
      transcriptRow$type <- "transcript"
      newTranscript <- rbind(transcriptRow,
                             newTranscript)

      newTranscript$seqnames <- seqnames
      newTranscript$strand <- strand
      newTranscript$gene_id <- gene_id
      newTranscript$gene_type <- gene_type
      newTranscript$gene_name <- gene_name
      newTranscript$transcript_id <- paste0(gene_id, "novelExon")
      newTranscript$transcript_type <- "novelExon"

      newTranscript[, ! colnames(newTranscript) %in% c("seqnames", "start", "end", "width", "strand", "type",
                                                       "gene_id", "gene_type", "gene_name", "transcript_id",
                                                       "exon_number", "exon_id")] <- NA

      return(rbind(existingTranscripts, newTranscript))

    }
  }


  warning("no new exon could be added for gene ", gene_id)
  return(existingTranscripts)

}

###

#' @title addNewExons()
#'
#' @description generate a new transcript containing a novel exon for each gene in gene_ids
#'
#' @param annot_gtf table containing gtf info for annotated transcripts
#' @param gene_ids genes for which transcript with novel exon should be added
#' @param canonical new exons are canonical (GT donor site, AG acceptor site) [default = TRUE]
#' @param seed seed for reproducability [default = 1234]
#'
#' @return the updated table containing gtf info for the annotated transcripts and the novel transcripts
#' @export
#' @importFrom Biostrings DNAString subseq vmatchPattern
#' @importFrom IRanges end start
addNewExons <- function(annot_gtf, gene_ids, fasta, canonical = TRUE, seed = 1234){

  if (! all(gene_ids %in% annot_gtf$gene_id))
    stop("not all gene_ids are present in annot_gtf")

  print("adding de novo exons")
  set.seed(seed)

  #get ref fasta
  ref_seq = readDNAStringSet(fasta)

  #reset rownames
  rownames(annot_gtf) <- 1:nrow(annot_gtf)

  for (gene_id in gene_ids){

    #get all gtf info for this gene
    gene_gtf <- annot_gtf[annot_gtf$gene_id == gene_id, ]

    #get all transcripts for this gene, including additional transcript with 1 novel exon
    geneTranscripts <- addNewExonSingleGene(gene_gtf, ref_seq[unique(gene_gtf$seqnames)], canonical)

    #replace old gene gtf info with novel gene gtf info
    annot_gtf <- rbind(annot_gtf[as.numeric(rownames(annot_gtf)) < min(as.numeric(rownames(gene_gtf))), ],
                       geneTranscripts,
                       annot_gtf[as.numeric(rownames(annot_gtf)) > max(as.numeric(rownames(gene_gtf))), ])

    #reset rownames
    rownames(annot_gtf) <- 1:nrow(annot_gtf)

  }

  return(annot_gtf)

}



