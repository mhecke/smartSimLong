
#set path to smartSim
ssPath <- "~/Documents/smartSimLong" #TODO: change to path/to/smartSim

#install from source
install.packages(ssPath, repos = NULL, type="source")

library(smartSimLong)
library(reticulate)
library(R.utils)
setwd(ssPath)

# use smartSim conda environment
use_condaenv("smartSim")

#input data
shortBam = "/home/bioit/mhecke/Documents/jurkat_241209/zUMIs/jurkat.filtered.Aligned.GeneTagged.sorted.bam"
longBamDir = "/home/bioit/mhecke/Documents/jurkat_241209/UMI_tools/UMIs_longBarcodes"
gtfFile <- "/home/bioit/mhecke/Documents/refGenomeNoReadThrough/noReadThrough.gtf"
fastaFile <- "/home/bioit/mhecke/Documents/refGenomeNoReadThrough/GRCh38.primary_assembly.genome.fa"

#output directories
dataCharDir <- "example_data_long/dataChar"
trueTranscriptDir <- "example_data_long/trueTranscripts"
readDir <- "example_data_long/reads"

#parameters
nClust = 2            #number of cell populations
nCellsPerClust = 10   #number of cells per cell population
nGenes = 10           #number of genes to simulate
nTx = 2               #number of transcripts per gene
barcodeLen = 10       #length of cell barcode
UMIlen = 8            #length of UMI
readLen = 150         #length of 1 end of PE read
tag = "ATTGCGCAATG"   #Smart-seq3 tag used to identify UMI containing reads

#for reproducibility
set.seed(1234)

if (! file.exists(fastaFile) & file.exists(paste0(fastaFile, ".gz"))) {
    gunzip(paste0(fastaFile, ".gz"))
}
if(! file.exists(gtfFile) & file.exists(paste0(gtfFile, ".gz"))){
  gunzip(paste0(gtfFile, ".gz"))
}

#characterize existing dataset
dataChar <- characterizeData(shortBam,
                             gtfFile,
                             fastaFile,
                             dataCharDir,
                             longBam = longBamDir,
                             nCPU = 4,
                             numReads = 1000000,
                             seed = 1234,
                             stage = "geneFilt",
                             phred = 33)

# this may take a while.
# To save time you can also read in the datacharacteristics provided, by uncommenting the following line:
#dataChar <- DataCharacteristics(dataCharDir)

#we get 2 transcripts for 10 random genes (with at least 2 annotated transcripts and 1000 bp length) from the reference annotation
gtf_simul <- getTranscripts(gtf = gtfFile,
                            fasta = fastaFile,
                            nTx = nTx,
                            out = trueTranscriptDir,
                            nRandomGenes = nGenes,
                            minTxLen = 1000,
                            seed = 1234)

#add novel junctions for the first 2 genes
genes_nj <- unique(gtf_simul$gene_id)[1:2]
gtf_simul <- addNewJunctions(gtf_simul,
                             genes_nj,
                             seed = 1234)

# add novel exons for the next 2 genes
genes_ne <- unique(gtf_simul$gene_id)[3:4]
gtf_simul <- addNewExons(annot_gtf = gtf_simul,
                         fasta = fastaFile,
                         gene_ids = genes_ne,
                         canonical = TRUE,
                         seed = 1234)

#write ground truth to file
transcript_seq = writeTrueTranscripts(gtf = gtf_simul,
                                      fasta = fastaFile,
                                      out = trueTranscriptDir)

#get gene expression similar to data
geneExpression <- getRandomGeneExpression(nClusters = nClust,
                                          nCellsPerCluster = nCellsPerClust,
                                          genes = unique(gtf_simul$gene_id),
                                          dataChar = dataChar,
                                          minCount = nTx,
                                          seed = 1234)

#assign random PSI values
geneIDs <- unique(gtf_simul$gene_id)
fixedGenes <- geneIDs[1:(length(geneIDs)/2)]
diffGenes <- geneIDs[((length(geneIDs)/2)+1):length(geneIDs)]
clusterPSI <- getRandomClusterPSI(gtf = gtf_simul,
                                  nClusters = nClust,
                                  fixedGenes = fixedGenes,
                                  diffGenes = diffGenes,
                                  seed = 1234)

#create transcript count matrix
transcriptExpression <- merge(clusterPSI, geneExpression, by = c("gene_id", "cluster"))
transcriptExpression$transcriptExpression <- round(transcriptExpression$geneExpression * transcriptExpression$PSI)
write.table(transcriptExpression, paste0(trueTranscriptDir, "/transcriptExpression.tsv"), 
            sep = "\t")

#simulate reads
simulateReads(transcriptExpression = transcriptExpression,
              transcript_seq = transcript_seq,
              dataChar = dataChar,
              out = readDir,
              barcodeLen = barcodeLen,
              UMIlen = UMIlen,
              readLen = readLen,
              tag = tag,
              PAF = 2,
              errorRate = 0.005,
              long = TRUE,
              longErrorRate = 0.05,
              seed = 1234)


