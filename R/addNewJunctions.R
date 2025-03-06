
### helper functions
getJunction <- function(row1, row2){

  #make sure exons are not overlapping
  if (max(row1$start, row2$start) <= min(row1$end, row2$end))
    return(NULL)

  start <- min(row1$end, row2$end)
  end <- max(row2$start, row1$start)

  if (row1$start < row2$start)
    return(list("start" = start, "end" = end, "startExon" = row1, "endExon" = row2))
  else
    return(list("start" = start, "end" = end, "startExon" = row2, "endExon" = row1))

}

#add de novo junction for single gene
addNewJunctionSingleGene <- function(existingTranscripts){

  #all transcripts should be from the same gene
  if (length(unique(existingTranscripts$gene_id)) != 1){
    stop("make sure all transcripts in existingTranscripts are from the same gene")
  }
  strand = existingTranscripts$strand[1]
  seqnames = existingTranscripts$seqnames[1]
  gene_id = existingTranscripts$gene_id[1]
  gene_type = existingTranscripts$gene_type[1]
  gene_name = existingTranscripts$gene_name[1]

  #get exons
  exons <- existingTranscripts[existingTranscripts$type == "exon",]
  exons$exon_number <- as.numeric(exons$exon_number)

  #get existing junctions
  junctions <- exons %>%
    arrange(transcript_id, start) %>%
    group_by(transcript_id) %>%
    reframe(
      junctionStart = end[-n()],
      junctionEnd = start[-1]
    )

  exonPairs <- t(combn(nrow(exons), 2))
  exonPairs <- exonPairs[sample(nrow(exonPairs)), , drop = FALSE]

  for (i in 1:nrow(exonPairs)){
    ex1 <- exons[exonPairs[i, 1],]
    ex2 <- exons[exonPairs[i, 2],]

    junction <- getJunction(ex1, ex2)
    if(is.null(junction))
      next

    #check if corresponding junction exists in ref
    if (any(junctions$junctionStart == junction$start & junctions$junctionEnd == junction$end))
      next #only add junctions that do not yet exist

    #make transcript using this junction
    if(strand == "+"){
      firstExons <- exons[exons$transcript_id == junction$startExon$transcript_id & exons$exon_number <= junction$startExon$exon_number,]
      lastExons <- exons[exons$transcript_id == junction$endExon$transcript_id & exons$exon_number >= junction$endExon$exon_number,]
    }else{
      firstExons <- exons[exons$transcript_id == junction$endExon$transcript_id & exons$exon_number <= junction$endExon$exon_number,]
      lastExons <- exons[exons$transcript_id == junction$startExon$transcript_id & exons$exon_number >= junction$startExon$exon_number,]
    }

    newTranscript <- rbind(firstExons, lastExons)
    newTranscript$exon_number <- 1:nrow(newTranscript)

    #add transcript header
    transcriptRow <- data.frame(matrix(NA, nrow = 1, ncol = ncol(newTranscript)))
    colnames(transcriptRow) <- colnames(newTranscript)
    transcriptRow$seqnames <- seqnames
    transcriptRow$start = min(newTranscript$start)
    transcriptRow$end = max(newTranscript$end)
    transcriptRow$width = transcriptRow$end - transcriptRow$start + 1
    transcriptRow$gene_name = gene_name
    transcriptRow$strand <- strand
    transcriptRow$type <- "transcript"
    transcriptRow$gene_id <- gene_id
    transcriptRow$gene_type <- gene_type
    newTranscript <- rbind(transcriptRow,
                           newTranscript)

    newTranscript$transcript_id <- paste0(gene_id, "newJunction")
    newTranscript$transcript_type <- "novelJunction"

    newTranscript[, ! colnames(newTranscript) %in% c("seqnames", "start", "end", "width", "strand", "type",
                                                     "gene_id", "gene_type", "gene_name", "transcript_id",
                                                     "exon_number", "exon_id")] <- NA

    return(rbind(existingTranscripts, newTranscript))

  }

  warning("no new junction could be added for gene ", gene_id)
  return(existingTranscripts)

}

###


#' @title addNewJunctions()
#'
#' @description generate a new transcript containing a novel junction for each gene in gene_ids
#'
#' @param annot_gtf table containing gtf info for annotated transcripts
#' @param gene_ids genes for which transcript with novel junction should be added
#' @param seed seed for reproducability [default = 1234]
#'
#' @return the updated table containing gtf info for the annotated transcripts and the novel transcripts
#' @export
#' @importFrom dplyr %>% arrange n reframe
addNewJunctions <- function(annot_gtf, gene_ids, seed = 1234){

  if (! all(gene_ids %in% annot_gtf$gene_id))
    stop("not all gene_ids are present in annot_gtf")

  print("adding de novo junctions")
  set.seed(seed)

  #reset rownames
  rownames(annot_gtf) <- 1:nrow(annot_gtf)

  for (gene_id in gene_ids){

    #get all gtf info for this gene
    gene_gtf <- annot_gtf[annot_gtf$gene_id == gene_id, ]

    #get all transcripts for this gene, including additional transcript with 1 novel junction
    geneTranscripts <- addNewJunctionSingleGene(gene_gtf)

    #replace old gene gtf info with novel gene gtf info
    annot_gtf <- rbind(annot_gtf[as.numeric(rownames(annot_gtf)) < min(as.numeric(rownames(gene_gtf))), ],
                       geneTranscripts,
                       annot_gtf[as.numeric(rownames(annot_gtf)) > max(as.numeric(rownames(gene_gtf))), ])

    #reset rownames
    rownames(annot_gtf) <- 1:nrow(annot_gtf)

  }

  return(annot_gtf)

}



