#' @title getRandomClusterPSI()
#'
#' @description generate random cluster PSI values
#'
#' @param gtf table containing gtf info for all transcripts
#' @param nClusters number of clusters
#' @param fixedGenes genes that have the same PSI values across clusters
#' @param diffGenes genes that have different PSI values across clusters
#' @param seed for reproducibility [default = 1234]
#'
#' @return matrix containing PSI value for each cluster and transcript
#'
#' @export
#' @importFrom dplyr %>% group_by mutate bind_rows
getRandomClusterPSI <- function(gtf, nClusters, fixedGenes, diffGenes, seed){

  set.seed(seed)

  #check if genes in gtf
  if (! all(fixedGenes %in% gtf$gene_id))
    stop("not all fixedGenes are present in gtf: missing ", fixedGenes[! fixedGenes %in% gtf$gene_id])
  if (! all(diffGenes %in% gtf$gene_id))
    stop("not all diffGenes are present in gtf: missing ", diffGenes[! diff %in% gtf$gene_id])
  if (length(intersect(fixedGenes, diffGenes)) > 0)
    stop("fixedGenes and diffGenes should not overlap, intersection contains: ", intersect(fixedGenes, diffGenes))


  #fixed transcripts have same PSI for all clusters
  fixed_transcripts <- gtf[gtf$type == "transcript" & gtf$gene_id %in% fixedGenes, c("gene_id", "transcript_id")]
  fixed_transcripts$PSI <- runif(nrow(fixed_transcripts)) #generate random values
  fixed_transcripts <- fixed_transcripts %>%
      group_by(gene_id) %>%
      mutate(PSI = PSI/sum(PSI))
  fixed_transcripts <- bind_rows(replicate(nClusters, fixed_transcripts, simplify = FALSE), .id = "cluster") #repeat for all clusters

  diff_transcripts <- gtf[gtf$type == "transcript" & gtf$gene_id %in% diffGenes, c("gene_id", "transcript_id")]
  diff_transcripts <- bind_rows(replicate(nClusters, diff_transcripts, simplify = FALSE), .id = "cluster") #repeat for each cluster
  diff_transcripts$PSI <- runif(nrow(diff_transcripts)) #generate random values
  diff_transcripts <- diff_transcripts %>%
      group_by(cluster, gene_id) %>%
      mutate(PSI = PSI/sum(PSI))

  return(rbind(fixed_transcripts, diff_transcripts))

}

