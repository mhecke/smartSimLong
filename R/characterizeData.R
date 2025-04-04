#' @title characterizeData()
#'
#' @description characterize Smart-seq3 dataset
#'
#' @param shortBam a single bam file, or a directory containing several short read bam files used to characterize data
#' @param gtf gtfFile containing the reference annotation
#' @param fasta fastaFile containing the reference genome
#' @param out output directory
#' @param longBam a single bam file, or a directory containing several long read bam files, used to characterize long read data (when available)
#' @param nCPU number of parallel CPU processes [default = 4]
#' @param numReads max number of reads in short bam file to be used for data characterization, by default all reads are used
#' @param nBam number of short bam files used for data characterization, only used if multiple bam files are used, by default all bam files are used
#' @param seed random seed for reproducibility [default = 1234]
#' @param stage stage from where to start running: geneFilt, readBams, plotStats [default = geneFilt]
#' @param phred phred encoding [default = 33]
#'
#' @return object of class "DataCharacteristics"
#'
#' @export
characterizeData <- function(shortBam, gtf, fasta, out, longBam = NULL, nCPU = 4, numReads = NULL, nBam = NULL, seed = 1234, stage = "geneFilt", phred = 33) {

  #checks
  if (! stage %in% c("geneFilt", "readBams", "plotStats")){
    print(paste0(stage, " is not a valid value for stage."))
    print("The options are geneFilt, readBams or plotStats")
    stop(cat(paste0(stage, " is not a valid value for stage.", "\n", "The options are geneFilt, readBams or plotStats")))
  }
  if (stage == "readBams" | stage == "plotStats"){

    if (! file.exists(paste0(out, "/exons.bin")) | ! file.exists(paste0(out, "/seqGC.bin")) | ! file.exists(paste0(out, "/geneGC.bin")) |
        ! file.exists(paste0(out, "/geneLen.bin")) | ! file.exists(paste0(out , "/uniqueTranscriptGenes.bin")))
      stop("geneFilt stage did not run successfully, please run characterizeData with stage = geneFilt")

    if (stage == "plotStats"){
      if (! file.exists(paste0(out, "/fragmentLength.tsv")))
        stop(cat(paste0("missing fragmentLength.tsv", "\n",
            "readBams stage did not run successfully, please run characterizeData with stage = readBams")))
      if (! file.exists(paste0(out, "/fragmentLength_UMI.tsv")))
        stop(cat(paste0("missing fragmentLength_UMI.tsv", "\n",
                        "readBams stage did not run successfully, please run characterizeData with stage = readBams")))
      if (! file.exists(paste0(out, "/fragmentLength_nonUMI.tsv")))
        stop(cat(paste0("missing fragmentLength_nonUMI.tsv", "\n",
                        "readBams stage did not run successfully, please run characterizeData with stage = readBams")))
      if (! file.exists(paste0(out, "/fragGeneLength.tsv")))
        stop(cat(paste0("missing fragGeneLength.tsv", "\n",
                        "readBams stage did not run successfully, please run characterizeData with stage = readBams")))
      if (! file.exists(paste0(out, "/GCperc.tsv")))
        stop(cat(paste0("missing GCperc.tsv", "\n",
                        "readBams stage did not run successfully, please run characterizeData with stage = readBams")))

    }
  } else{
    dir.create(out, showWarnings = FALSE)
  }

  if (! dir.exists(shortBam) & ! file.exists(shortBam)){
    print("provide valid shortBam file/directory")
    stop()
  }

  python_cmd <- system.file("scripts", "characterizeData.py", package = "smartSim")

  if (dir.exists(shortBam)){
    python_cmd = paste(python_cmd,
                       "--bamDir", shortBam)
  } else{
    python_cmd = paste(python_cmd,
                       "--bamFile", shortBam)
  }
  if(! is.null(longBam)){
    if (dir.exists(longBam)){
      python_cmd = paste(python_cmd, 
                         "--longBamDir", longBam)
    } else{
      python_cmd = paste(python_cmd, 
                         "--longBamFile", longBam)
    }
  }
  
  if (! is.null(nBam)){
    python_cmd = paste(python_cmd,
                       "--nBam", nBam)
  }
  if (! is.null(numReads)){
    python_cmd = paste(python_cmd,
                       "--numReads", format(numReads, scientific = FALSE, trim = TRUE))
  }

  python_cmd = paste(python_cmd,
                  "--gtf", gtf,
                  "--fasta", fasta,
                  "--out", out,
                  "--nCPU", nCPU,
                  "--seed", seed,
                  "--stage", stage,
                  "--phred", phred,
                  ">", paste0(out, "/out.log"))

  print(python_cmd)

  conda_env = "smartSim"
  system2("conda", args = c("run", "-n", conda_env, "python", python_cmd))

  #read data characteristics & return object
  return(DataCharacteristics(dataCharDir))

}


