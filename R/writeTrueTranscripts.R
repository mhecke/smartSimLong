#' @title writeTrueTranscripts()
#'
#' @description write ground truth to files (gtf annotation and fasta with sequence for each transcript)
#'
#' @param gtf gtf table containing ground truth transcripts
#' @param fasta fastaFile containing the reference genome
#' @param out outdir to write resulting gtf and fasta file
#' @return transcript sequences as DNASTringSet object
#'
#' @importFrom rtracklayer export
#' @importFrom Biostrings readDNAStringSet DNAStringSet reverseComplement writeXStringSet
#'
#' @export
writeTrueTranscripts <- function(gtf, fasta, out){
  
  print("writing true transcripts to file")
  dir.create(out, showWarnings = FALSE)

  #write gtf file
  export(gtf, paste0(out, "/true.gtf"))

  #get sequence info from transcripts in gtf
  ref_seq = readDNAStringSet(fasta)
  
  #in case fasta file contains space separated chromosome names and chromosome numbers
  ref_seq_names = as.data.frame(do.call(rbind, strsplit(names(ref_seq), split = " ")))
  selectCol <- apply(ref_seq_names, 2, function(col) all(unique(gtf$seqnames) %in% col))
  if (! any(selectCol)){
    stop("fasta file does not contain all chromosomes in gtf.")
  }
  names(ref_seq) <- ref_seq_names[,selectCol]
  
  #get transcript sequences
  ## code adapted from polyester (https://github.com/alyssafrazee/polyester-release/blob/master/R/seq_gtf.R)
  chrs = unique(gtf$seqname)
  seqlist = lapply(chrs, function(chr){
    dftmp = gtf[(gtf$seqname==chr) & (gtf$type == "exon"),]
    fullseq = ref_seq[which(names(ref_seq) == chr)]
    these_seqs = subseq(rep(fullseq, times=nrow(dftmp)), start=dftmp$start,
                        end=dftmp$end)
    names(these_seqs) = dftmp$transcript_id
    these_seqs
  })
  #get exons for all chromosomes
  full_list = do.call(c, seqlist)
  
  #reverse complement for neg strand
  full_list[names(full_list) %in% gtf$transcript_id[gtf$strand == "-"]] = reverseComplement(full_list[names(full_list) %in% gtf$transcript_id[gtf$strand == "-"]])

  #join sequences per transcript
  split_list = split(full_list, names(full_list))
  transcript_fasta <- DNAStringSet(lapply(split_list, unlist))
  
  writeXStringSet(transcript_fasta, paste0(out, "/transcripts.fa"))
  return(transcript_fasta)
}
