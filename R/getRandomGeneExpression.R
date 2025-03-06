#' @title getRandomGeneExpression()
#'
#' @description generate random gene expression based on gene expression from data.
#'   gene expression for cells within the same cluster is only different by cellSize
#'
#' @param nClusters number of clusters
#' @param nCellsPerCluster number of cells per cluster
#' @param genes geneIDs to generate read counts for
#' @param dataChar data characteristics object
#' @param minCount minimum gene expression [default = 0]
#' @param seed for reproducability [default = 1234]
#'
#' @return matrix containing gene expression for each cell
#'
#' @export
#get random gene expression for each cell
getRandomGeneExpression <- function(nClusters, nCellsPerCluster, genes, dataChar, minCount = 0, seed = 1234){

  print("getting random gene expression")
  set.seed(seed)

  if (nClusters < 1 | nClusters%%1 != 0)
    stop("nClusters = ", nClusters, " is not a valid number")
  if (nCellsPerCluster < 1 | nCellsPerCluster%%1 != 0)
    stop("nCellsPerCluster = ", nCellsPerCluster, " is not a valid number")

  #get gene expression for cluster
  clusterExpr <- sampleUMIcounts(dataChar, n = length(genes) * nClusters, minCount = minCount)
  clusterExpr <- as.data.frame(clusterExpr, nrow = length(genes) * nClusters, ncol = 1)
  colnames(clusterExpr) <- c("clusterExpression")
  clusterExpr$gene_id <- rep(genes, nClusters)
  clusterExpr$cluster <- rep(1:nClusters, each = length(genes))

  #sample cellsize for each cell within the same cluster (cellSize is between 1/4 and 4)
  cellSize <- as.data.frame(2^runif(nCellsPerClust*nClusters, min = -2, max = 2),
                            nrow = nCellsPerClust * nClusters, ncol = 1)
  colnames(cellSize) <- "cellSize"
  cellSize$cell <- paste0("cluster", rep(1:nClusters, each = nCellsPerClust), "_cell", rep(1:nCellsPerClust, nClusters))
  cellSize$cluster <- rep(1:nClusters, each = nCellsPerClust)
  
  geneExpr <- merge(clusterExpr, cellSize, by = "cluster")
  geneExpr$geneExpression <- geneExpr$clusterExpr * geneExpr$cellSize
  geneExpr <- geneExpr[, c("cluster", "cell", "gene_id", "geneExpression")]

  return(geneExpr)

}
