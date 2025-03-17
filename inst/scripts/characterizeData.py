#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pickle
import sys
import HTSeq
import copy
from os import listdir
from os import makedirs
import random
import re
import multiprocessing
import argparse
import warnings
from collections import Counter
from functools import partial
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import numpy as np
import statsmodels.api as sm
import csv
#import git

lowess = sm.nonparametric.lowess

#remove warnings for mate pairs that are not found in BAM files
warnings.filterwarnings(action='ignore', category=UserWarning) 
warnings.filterwarnings(action='ignore', category=FutureWarning)

def getStandardGenes(gtfFile, outDir):

    print("getting non-overlapping genes")

    keepFile = outDir + "/genesToKeep.tsv"
    
    #read GTF
    gtf = HTSeq.GFF_Reader(gtfFile)
    
    #known chromosomes 
    chroms = []
    for i in range(1,23):
        chroms.append("chr" + str(i))
    chroms.append("chrX")
    chroms.append("chrY")
    
    #read in annotation
    geneIntervals = {}
    allGenes = set()
    gas = HTSeq.GenomicArrayOfSets(chroms, stranded=True)
    for feature in gtf:
        if feature.type == "gene":
            try:
                gas[feature.iv] += feature.name
                geneIntervals[feature.name] = feature.iv
                allGenes.add(feature.name)
            except KeyError:
                #ignore genes that are not located on a known chromosome
                pass
                
    
    #find genes that are completely contained in another gene 
    genesContained = set()
    for iv, step_set in gas.steps():
        if len(step_set) > 1:
            for g1 in step_set:
                for g2 in step_set:
                    if (g1 != g2 and geneIntervals[g1].contains(geneIntervals[g2])):
                        genesContained.add(g2)
    
    #remove genes that are completely contained in another gene
    for gene in genesContained:
        gene_iv = geneIntervals[gene]
        for iv, step_set in gas[gene_iv].steps():
            new_step_set = copy.deepcopy(step_set)
            new_step_set.remove(gene)
            gas.__setitem__(iv, new_step_set)
    
    #find genes that are still overlapping and only keep longest one
    genesToRemove = set() 
    for iv, step_set in gas.steps():
        if len(step_set) > 1:
            #keep longest one
            maxLen = 0
            maxGene = ""
            for gene in step_set:
                if geneIntervals[gene].length > maxLen:
                    maxGene = gene
            genesToRemove.update(step_set.difference({maxGene}))
    
    #remove genes that are still overlapping with other genes
    for gene in genesToRemove:
        gene_iv = geneIntervals[gene]
        for iv, step_set in gas[gene_iv].steps():
            new_step_set = copy.deepcopy(step_set)
            new_step_set.remove(gene)
            gas.__setitem__(iv, new_step_set)
            
            
    print("Total number of genes: ", len(allGenes))
    print("Number of genes that are completely contained within another gene: ", len(genesContained))
    print("Number of genes that were removed because they still show overlap with another gene: ", len(genesToRemove))      
    print("In total, we removed ", len(genesContained) + len(genesToRemove), "genes")
    
    genesToKeep = allGenes.difference(genesContained, genesToRemove)
    
    #write to file
    f = open(keepFile, "w")
    for gene in genesToKeep:
        f.write(gene + "\n")
    f.close()
    
    return genesToKeep

def getNonOverlappingGenes(gtfFile, geneSet):
    #read GTF
    gtf = HTSeq.GFF_Reader(gtfFile)  
    
    #read gtf 
    gasUnstranded = HTSeq.GenomicArrayOfSets(chroms = "auto", stranded = False)
    for feature in gtf:
        if feature.type == "exon":
            #only look at genes to keep
            if not feature.name in geneSet:
                continue
            gasUnstranded[feature.iv] += feature.name
    
    #get genes that overlap with another gene on the complementary strand
    overlappingGenes = set()
    for iv, step_set in gasUnstranded.steps():
        if len(step_set) > 1:
            overlappingGenes.update(step_set)
            
    print(len(overlappingGenes), " out of ", len(geneSet), " genes show overlap with a gene on the complementary strand")
    
    nonOverlappingGenes = geneSet.difference(overlappingGenes)
    return nonOverlappingGenes

def getUniqueTranscriptGenes(gtfFile, geneSet):
    
    uniqueTranscriptGenes = set()
    
    #read GTF
    gtf = HTSeq.GFF_Reader(gtfFile)   
    prev_gene = ""
    numTranscripts = 0
    for feature in gtf:
        #only look at genes to keep
        if not feature.name in geneSet:
            continue
        
        #count number of transcripts per gene
        if feature.type == "gene":
            #check if previous gene has 1 unique transcript
            if numTranscripts == 1:
                uniqueTranscriptGenes.add(prev_gene)
            prev_gene = feature.name
            numTranscripts = 0
        if feature.type == "transcript":
            numTranscripts += 1
    
    print("Number of non-overlapping genes with only 1 transcript: ", len(uniqueTranscriptGenes))
    return uniqueTranscriptGenes

def getGCcontent(gtfFile, fastaFile, geneSet):

    geneGC = {gene:0 for gene in geneSet}   #GC % per gene
    geneLen = {gene:0 for gene in geneSet}  #total exon length per gene

    #read reference fasta file
    ref_fasta = HTSeq.FastaReader(fastaFile)
    seqDir = {}
    for chrom in ref_fasta:
        seqDir[chrom.name] = chrom.seq
    
    #read GTF
    gtf = HTSeq.GFF_Reader(gtfFile)  
    exons = HTSeq.GenomicArrayOfSets(chroms = "auto", stranded = False)
    seqGC = HTSeq.GenomicArray(chroms = "auto", stranded = False, typecode = 'b') #true if base is G/C
    for feature in gtf:

        if feature.type == "exon":
            #add exon to reference genomic array
            exons[feature.iv] += feature.name

            #get GC
            #only look at selection of genes
            if not feature.name in geneSet:
                continue

            #get GC for each
            for pos in range(feature.iv.start, feature.iv.end, 1):
               s = seqDir[feature.iv.chrom][pos:pos+1].decode("utf-8")
               unStrandedPos = HTSeq.GenomicPosition(feature.iv.chrom, pos, strand = '.')

               seqGC[unStrandedPos] = (s == "G" or s == "C")
               geneGC[feature.name] += (s == "G" or s == "C")

            #total exon length per gene
            geneLen[feature.name] += feature.iv.length
    
    #remove seqDir from memory: typically large
    del seqDir

    #GC % per gene
    for gene, count in geneGC.items():
        geneGC[gene] = (count / geneLen[gene]) * 100

    return exons, seqGC, geneGC, geneLen

def readBams(exons, seqGC, geneGC, geneLen, bamFiles, numReads, outDir, nCPU, geneSet):

    # total fragment length distr.
    allFragLen = Counter()
    allFragGeneLen = Counter()
    UMIFragLen = Counter()
    nonUMIFragLen = Counter()

    #total GC distribution
    allGC = Counter()

    #observed GC values
    obsGC = {gene:[] for gene in geneSet}
    
    #UMI counts
    allUMIcounts = Counter()
    
    #quality scores per position
    allReadQual = []
    allBcQual = []
    
    #statsLock = multiprocessing.Lock()
        
    print("reading bam files")


    #read bamfiles in parallel
    nCPU = min(len(bamFiles), nCPU)
    with multiprocessing.Pool(nCPU) as pool:
        res = pool.imap(partial(readSingleBam, exons = exons, seqGC = seqGC, geneSet = geneSet, geneLen = geneLen, numReads = numReads), bamFiles)
        pool.close()
        pool.join()

        #merge counts from different files
        for fraglen, fragGeneLen, fragLen_UMI, fragLen_nonUMI, GC, GC_geneAvg, UMIcounts, readQual, bcQual in res:
            allFragLen.update(fraglen)
            allFragGeneLen.update(fragGeneLen)
            UMIFragLen.update(fragLen_UMI)
            nonUMIFragLen.update(fragLen_nonUMI)
            allGC.update(GC)
            for gene, avgGC in GC_geneAvg.items():
                obsGC[gene].append(avgGC)
            allUMIcounts.update(UMIcounts)
            
            if (len(allReadQual) < len(readQual)):
                allReadQual += [Counter() for _ in range(len(readQual)-len(allReadQual))]
            for i, qualCtr in enumerate(readQual):
                allReadQual[i].update(qualCtr)
            
            if (len(allBcQual) < len(bcQual)):
                allBcQual += [Counter() for _ in range(len(bcQual)-len(allBcQual))]
            for i, qualCtr in enumerate(bcQual):
                allBcQual[i].update(qualCtr)

    #write fragmentlength to file
    f = open(outDir + "/fragmentLength.tsv", "w")
    f.write("val" + "\t" + "count" + "\n")
    for val, count in allFragLen.items():
        f.write(str(val) + "\t" + str(count) + "\n")
    f.close()
    
    f = open(outDir + "/fragGeneLength.tsv", "w")
    f.write("fragLen" + "\t" + "geneLen" + "\t" + "count" + "\n")
    for (fLen, gLen), count in allFragGeneLen.items():
        f.write(str(fLen) + "\t" + str(gLen) + "\t" + str(count) + "\n" )
    f.close()
    
    f = open(outDir + "/fragmentLength_UMI.tsv", "w")
    f.write("val" + "\t" + "count" + "\n")
    for val, count in UMIFragLen.items():
        f.write(str(val) + "\t" + str(count) + "\n")
    f.close()
    
    f = open(outDir + "/fragmentLength_nonUMI.tsv", "w")
    f.write("val" + "\t" + "count" + "\n")
    for val, count in nonUMIFragLen.items():
        f.write(str(val) + "\t" + str(count) + "\n")
    f.close()

    #write GC contents
    f = open(outDir + "/GCperc.tsv", "w")
    f.write("val" + "\t" + "count" + "\n")
    for val, count in allGC.items():
        f.write(str(val) + "\t" + str(count) + "\n")
    f.close()

    #write average observed GC % per gene
    f = open(outDir + "/obsGeneGC.tsv", "w")
    f.write("geneGC" + "\t" + "obsGC" + "\n")
    for gene, GClist in obsGC.items():
        for obsGC in GClist:
            f.write(str(geneGC[gene]) + "\t" + str(obsGC) + "\n")
    f.close()
    
    f = open(outDir + "/UMIcounts.tsv", "w")
    f.write("cell"+ "\t"+ "gene"+ "\t"+ "UMI"+ "\t"+ "numReads"+ "\n")
    for (BC, gene, UMI), count in allUMIcounts.items():
        f.write(BC+ "\t" + gene + "\t" + UMI + "\t" + str(count) + "\n")
    f.close()

    f = open(outDir + "/readQual.tsv", "w")
    f.write("pos" + "\t" + "qual" + "\t" + "count" + "\n")
    for pos, qualCtr in enumerate(allReadQual):
        for q, count in qualCtr.items():
            f.write(str(pos)  +"\t" + q + "\t" + str(count) + "\n")

    f = open(outDir + "/bcQual.tsv", "w")
    f.write("pos" + "\t" + "qual" + "\t" + "count" + "\n")
    for pos, qualCtr in enumerate(allBcQual):
        for q, count in qualCtr.items():
            f.write(str(pos)  +"\t" + q + "\t" + str(count) + "\n")


    
def readSingleBam(file, exons, seqGC, geneSet, geneLen, numReads):

    print("reading bam file: " + file)

    bam_reader = HTSeq.BAM_Reader(file)
    fragLen_counter = Counter()     #counts per fragment length in sample
    fragGeneLen_counter = Counter() #count fragment length/gene length pairs
    UMIFragLen_counter = Counter()      #counts per fragment length for UMI containing reads
    nonUMIFragLen_counter = Counter()   #counts per fragment length for non UMI containing reads
    GC_counter = Counter()          #counts per GC percentage in sample
    gene_counter = Counter()        #nr of reads mapping to gene
    GC_geneAvg = {}
    UMIcounts = Counter()
    fragPos = {}
    
    readQual = []
    bcQual = []

    ctr = 0
    for r1,r2 in HTSeq.pair_SAM_alignments_with_buffer(bam_reader):

        #check if both reads exist
        if (r1 is None or r2 is None):
            continue
            
        #skip non primary alignments
        if (r1.not_primary_alignment or r2.not_primary_alignment):
            continue
        
        BC = r1.optional_field('BC')
        if BC == "":
            break
        
        #limit number of reads used for data characterization
        if (numReads is not None) and (ctr > numReads):
            break
        ctr += 1
            
        UMI = r1.optional_field('UB')
        if UMI == "":
            UMI = "noUMI"
        
        try:
            #assign to gene (get non-overlapping set of genes)
            iset = assignToGene(r1,r2,exons)
                        
            #skip reads that don't map to a gene
            if iset is None or len(iset) == 0: 
                continue
            
            #only 1 gene assigned
            gene = iset.pop()

            #remove genes not in geneSet
            if gene not in geneSet:
                continue

            gene_counter[gene] += 1
            
            #get edges of fragment
            minPos = min(r1.iv.start, r2.iv.start)
            maxPos = max(r1.iv.end, r2.iv.end)
            
            if gene not in fragPos:
                fragPos[gene] = {"UMI" : [], "noUMI" : []}
            if UMI == "noUMI":
                fragPos[gene]["noUMI"].append(str(minPos) + "-" + str(maxPos))
            else: 
                fragPos[gene]["UMI"].append(str(minPos) + "-" + str(maxPos))
                
            
            #compute fragmentlength
            fragment_iv = HTSeq.GenomicInterval(r1.iv.chrom, minPos, maxPos, ".")
            fragment_length = 0
            GC_count = 0
            for iv, step_set in exons[fragment_iv].steps():
                if gene in step_set:
                    fragment_length += iv.length
                    GC_count += sum(seqGC[iv])
            
            #keep fragment length
            fragLen_counter[fragment_length] += 1
            fragGeneLen_counter[(fragment_length, geneLen[gene])] += 1
            
            if UMI == "noUMI":
                nonUMIFragLen_counter[fragment_length] += 1
            else:
                UMIFragLen_counter[fragment_length] += 1
            
            #keep GC percentage
            GC_perc = round(GC_count/fragment_length * 100,2)
            GC_counter[GC_perc] += 1

            if gene not in GC_geneAvg:
                GC_geneAvg[gene] = 0
            GC_geneAvg[gene] += GC_perc
            
            UMIcounts[tuple((BC, gene, UMI))] += 1
            
            #keep quality per position for read
            thisReadLen = max(len(r1.read), len(r2.read))
            if thisReadLen > len(readQual):
                readQual += [Counter() for _ in range(thisReadLen - len(readQual))]
            for i, q in enumerate(r1.read.qualstr.decode('utf-8')):
                readQual[i][q] += 1
            for i, q in enumerate(r2.read.qualstr.decode('utf-8')):
                readQual[i][q] += 1
                
            #keep quality per position
            if (r1.optional_field("QB") != "" and r2.optional_field("QB") != ""):
                thisBCLen = len(r1.optional_field("QB")) + len(r2.optional_field("QB")) 
                if thisBCLen > len(bcQual):
                    bcQual += [Counter() for _ in range(thisBCLen - len(bcQual))]
                for i, q in enumerate(r1.optional_field("QB")):
                    bcQual[i][q] += 1
                for i, q in enumerate(r2.optional_field("QB")):
                    bcQual[i][q] += 1
                
        except KeyError:
            #chrom of read not in reference
            pass

    for gene in GC_geneAvg:
        GC_geneAvg[gene] /= gene_counter[gene]

    return fragLen_counter, fragGeneLen_counter, UMIFragLen_counter, nonUMIFragLen_counter, GC_counter, GC_geneAvg, UMIcounts, readQual, bcQual

def assignSingleReadToGene(read, exons):
    iset = None
    
    for cigop in read.cigar:
        #skip parts of alignment that are not alignment match, sequence match or sequence mismatch
        if not (cigop.type == "M" or cigop.type == "=" or cigop.type == "X") :
            continue
        for iv,step_set in exons[cigop.ref_iv].steps():
            if (len(step_set) > 0):
                if iset is None: #add to empty set
                    iset = step_set.copy()
                else:
                    iset.intersection_update(step_set) #only keep intersection of elements
    return iset

def assignToGene(read1, read2, exons):
    
    r1_genes = assignSingleReadToGene(read1, exons)
    r2_genes = assignSingleReadToGene(read2, exons)
    
    #both ends should be assigned to gene
    if (r1_genes is None) or (r2_genes is None):
        return None
    
    return (r1_genes.intersection(r2_genes))

def plotFragLen(inFile, outFile):
    #fragment length
    fragLen = pd.read_csv(inFile , sep = '\t')

    #sort values
    fragLen.sort_values(by='val', inplace=True)
    fragLen['cumCount'] = fragLen['count'].cumsum()

    #remove 5% largest fragments
    perc95 = 0.95 * fragLen['cumCount'].iloc[len(fragLen)-1]
    fragLen = fragLen.loc[fragLen['cumCount'] < perc95,:]
    
    #smoothing
    fragLen['smooth'] = (fragLen['count'].rolling(3, min_periods=1, center=True)
                                         .mean())

    #plot figure
    fig, ax = plt.subplots()
    ax.bar(fragLen['val'], fragLen['smooth'], width = 1)
    #ax.set_xlim(0,2500) 
    plt.title("fragment length")
    plt.xlabel("fragment length(bp)")
    plt.ylabel("number of occurences")

    #save figure
    plt.savefig(outFile)

def plotFragGeneLen(inFile, outFile):
    fragLen = pd.read_csv(inFile , sep = '\t')

    #sort values
    fragLen.sort_values(by='fragLen', inplace=True)
    fragLen['cumCount'] = fragLen['count'].cumsum()

    #remove 5% largest fragments
    perc95 = 0.95 * fragLen['cumCount'].iloc[len(fragLen)-1]
    fragLen = fragLen.loc[fragLen['cumCount'] < perc95,:]
    
    fit = lowess(fragLen['fragLen'], fragLen['geneLen'])

    #plot figure
    fig, ax = plt.subplots()
    ax.scatter(fragLen['geneLen'], fragLen['fragLen'])
    ax.plot(fit[:,0], fit[:,1], color='tomato', label='LOWESS')
    ax.set_ylim(0,2500)
    plt.title("fragment length")
    plt.xlabel("gene length(bp)")
    plt.ylabel("fragment length(bp)")

    #save figure
    plt.savefig(outFile)

def boxplotStats(df):
    df = df.sort_values(by = "qual")
    df["cumSum"] = df["count"].cumsum()
    totCount = df["cumSum"].iloc[-1]
    
    whislo = df["qual"].iloc[np.searchsorted(df["cumSum"], 0.05*totCount)]
    q1 = df["qual"].iloc[np.searchsorted(df["cumSum"], 0.25*totCount)]
    median =df["qual"].iloc[np.searchsorted(df["cumSum"], 0.50*totCount)]
    q3 = df["qual"].iloc[np.searchsorted(df["cumSum"], 0.75*totCount)]
    whishi = df["qual"].iloc[np.searchsorted(df["cumSum"], 0.95*totCount)]
    
    return [whislo, q1, median, q3, whishi]

def plotQual(inFile, outFile, phred):
    #read data
    readQual = pd.read_csv(inFile, sep = '\t', quoting=csv.QUOTE_NONE)
    readQual["qual"] = readQual["qual"].apply(ord) - phred
    readQual["pos"] += 1
    maxQual = max(readQual["qual"])
    
    #get boxplot info
    boxplot_df= readQual.groupby("pos")[["qual","count"]].apply(lambda x : boxplotStats(x))
    boxplot_df = pd.DataFrame(boxplot_df.tolist(), index=boxplot_df.index, columns=["whislo", "q1", "med", "q3", "whishi"])
    boxplot_df["label"] = boxplot_df.index
    boxplot_df["label"] = boxplot_df["label"].astype(str)
    boxplot_df = boxplot_df.transpose().to_dict().values()

    #get mean
    mean_df = readQual.groupby("pos")[["qual", "count"]].apply(lambda x : np.average(x["qual"], weights = x["count"]))

    # Create a boxplot
    fig, ax = plt.subplots(figsize=(15, 5))
    ax.bxp(boxplot_df, widths=0.5, showfliers=False)
    ax.set_ylim(0, maxQual)
    plt.plot(mean_df.index, mean_df, '-', color = "red")
    plt.xlabel("position")
    plt.ylabel("quality")
    plt.savefig(outFile)

def plotUMIcounts(inFile, outDir):
    #read data
    UMIcounts = pd.read_csv(inFile, sep = '\t')
    
    ##get number of unique UMIs per gene
    uniqueUMIs = UMIcounts[UMIcounts["UMI"] != "noUMI"]
    uniqueUMIs = uniqueUMIs.groupby(["cell", "gene"])["UMI"].nunique().reset_index(name="nUMIs")
    q99 = np.quantile(uniqueUMIs["nUMIs"], 0.99)
    
    #plot figure
    fig, ax = plt.subplots()
    ax.hist(uniqueUMIs["nUMIs"], bins = list(range(max(uniqueUMIs["nUMIs"]))))
    ax.set_xlim(0,q99)
    plt.title("unique UMIs per gene")
    plt.xlabel("number of unique UMIs per gene")
    plt.ylabel("count")

    #save figure
    plt.savefig(outDir + "/uniqueUMIs.png")
    
    
    ##get number of internal reads per gene
    internal = UMIcounts[UMIcounts["UMI"] == "noUMI"]
    q99 = np.quantile(internal["numReads"], 0.99)
    
    #plot figure
    fig, ax = plt.subplots()
    ax.hist(internal["numReads"], bins = list(range(max(internal["numReads"]))))
    ax.set_xlim(0,q99)
    #ax.set_xlim(0,140)
    plt.title("internal reads per gene")
    plt.xlabel("number of internal reads per gene")
    plt.ylabel("count")

    #save figure
    plt.savefig(outDir + "/internalReads.png")
    
    
    def computeInternalVsUMI(g): 
        if (g.loc[g["UMI"] != "noUMI", "numReads"].sum() == 0):
            return(np.nan)
        else:
            return (g.loc[g["UMI"] == "noUMI", "numReads"].sum() / g.loc[g["UMI"] != "noUMI", "numReads"].sum())

    ##get number of internal reads vs UMI containing reads per gene
    internal_vs_UMI = UMIcounts.groupby(["cell", "gene"]).apply(computeInternalVsUMI, include_groups=False).reset_index(name="val")
    internal_vs_UMI = internal_vs_UMI[~ internal_vs_UMI["val"].isna()]
    q99 = np.quantile(internal_vs_UMI["val"], 0.99)
    
    #plot figure
    fig, ax = plt.subplots()
    ax.hist(internal_vs_UMI["val"], bins = range(round(max(internal_vs_UMI["val"]))+1))
    ax.set_xlim(0,q99)
    #ax.hist(internal_vs_UMI["val"], bins = range(20))
    #ax.set_xlim(0,20)
    plt.title("internal reads per UMI containing read")
    plt.xlabel("number of internal reads vs UMI containing read")
    plt.ylabel("count")

    #save figure
    plt.savefig(outDir + "/internalReadsPerUMI.png")




def plotBias(outDir, phred):
    
    plotUMIcounts(outDir + "/UMIcounts.tsv", outDir)
    
    plotFragLen(outDir + "/fragmentLength.tsv", outDir + "/fragmentLength.png")
    plotFragLen(outDir + "/fragmentLength_UMI.tsv", outDir + "/fragmentLength_UMI.png")
    plotFragLen(outDir + "/fragmentLength_nonUMI.tsv", outDir + "/fragmentLength_nonUMI.png")
    
    plotFragGeneLen(outDir + "/fragGeneLength.tsv", outDir + "/fragmentGeneLength.png")

    #GC bias
    GC = pd.read_csv(outDir + "/GCperc.tsv", sep = '\t')

    #sort values
    GC.sort_values(by='val', inplace=True)
    GC['cumCount'] = GC['count'].cumsum()

    #plot figure
    fig, ax = plt.subplots()
    ax.bar(GC['val'], GC['count'])
    ax.set_xlim(0,100)
    plt.title("GC content")
    plt.xlabel("GC content (%)")
    plt.ylabel("number of occurences")

    #save figure
    plt.savefig(outDir + "/GC.png")

    #real geneGC vs observed GC
    #read observed GC
    obsGC = pd.read_csv(outDir + "/obsGeneGC.tsv" , sep = '\t')
    fig, ax = plt.subplots()
    ax.scatter(x = obsGC["geneGC"],y = obsGC["obsGC"])
    plt.title("GC bias")
    plt.xlabel("true GC content per gene (%)")
    plt.ylabel("observed GC content per gene (%)")
    plt.xlim(0,100)
    plt.ylim(0,100)
    plt.savefig(outDir + "/GCbias.png")

    #fit linear model to GCbias
    x = np.array(obsGC["geneGC"]).reshape((-1,1))
    y = np.array(obsGC["obsGC"])
    model = LinearRegression().fit(x, y)

    #plot figure with fitted model
    fig, ax = plt.subplots()
    ax.scatter(x = obsGC["geneGC"],y = obsGC["obsGC"])
    plt.title("GC bias")
    plt.xlabel("true GC content per gene (%)")
    plt.ylabel("observed GC content per gene (%)")
    plt.xlim(0,100)
    plt.ylim(0,100)
    x_vals = np.array(ax.get_xlim())
    y_vals = model.intercept_ + model.coef_ * x_vals
    plt.plot(x_vals, y_vals, '-', color = "red")

    plt.savefig(outDir + "/GCbias_fit.png")
    
    #plot read quality
    plotQual(outDir + "/readQual.tsv" , outDir + "/readQual.png", phred)
    plotQual(outDir + "/bcQual.tsv", outDir + "/bcQual.png", phred)



def _parse_arguments():
    pa = argparse.ArgumentParser(
        add_help=False
    )
    args, argv = pa.parse_known_args()
    
    pa = argparse.ArgumentParser(
        parents=[pa]
    )
    
    pa.add_argument(
        "--gtf",
        dest="gtfFile",
        type=str,
        help="path/to/gtfFile containing the reference annotation"
    )
    
    pa.add_argument(
        "--fasta",
        dest="fastaFile",
        type=str,
        help="path/to/fastaFile containing the reference genome"
        )
    
    pa.add_argument(
        "--out",
        dest="outDir",
        type = str,
        help="path/to/out where the output files will be written"
    )
    
    pa.add_argument(
        "--nCPU",
        dest="nCPU",
        type=int,
        default=4,
        help="number of parallel CPU processes (default = 4)",
    )
    
    pa.add_argument(
        "--bamFile", 
        dest="bamFile",
        type=str,
        help="path/to/bamFile containing bam file")
    
    pa.add_argument(
        "--numReads",
        dest="numReads",
        type=int,
        help="max number of reads in bam file to be used for data characterization, by default all reads are used")
    
    pa.add_argument(
        "--bamDir",
        dest="bamDir",
        type=str,
        help="path/to/directory containing multiple STAR output bam files")
    
    pa.add_argument(
        "--nBam",
        dest="nBam",
        type=int, 
        help="number of bam files used for bias computation only used if multiple bam files are used (--bamDir), by default all bam files are used")

    pa.add_argument(
        "--seed",
        dest="seed",
        type=int,
        help="seed for random bam selection")

    pa.add_argument(
        "--stage",
        dest="stage",
        type=str,
        default="geneFilt",
        help="start at specific stage. options: geneFilt, readBams, plotStats [default = geneFilt]"
    )
    
    pa.add_argument(
            "--phred", 
            dest="phred", 
            type=int, 
            default=33, 
            help="phred encoding, default using Phred+33 [default = 33]"
    )
    
    pa.add_argument(
            "--noGeneFilt",
            dest="noGeneFilt",
            action="store_true",
            default=False,
            help="do not filter genes with more than 1 transcript"
    )
    
    args = pa.parse_args()
    
    assert(args.stage in ['geneFilt', 'readBams', 'plotStats'])
    
    assert(args.bamDir is not None or args.bamFile is not None)
    
    if(args.bamDir is not None and args.bamFile is not None):
        raise ValueError("Please only provide --bamDir or --bamFile, not both")
        
    if (args.bamDir is None and args.nBam is not None):
        warnings.warn("nBam argument ignored, to limit the number of reads used, set --numreads argument")
    
    return args

def main(): 
    
    args = _parse_arguments()
    
    #repo = git.Repo(search_parent_directories=True)
    #commit_hash = repo.head.object.hexsha
    
    print("### getting bias from Smart-Seq3 experiment ###")
    #print("git commit: " + commit_hash)
    print(" ".join(sys.argv))
    print("###############################################")

    if (args.stage == 'geneFilt'):
        makedirs(args.outDir, exist_ok = True)

        #seed for random selection of bam files
        if args.seed is not None:
            random.seed(args.seed)

        #get standard genes
        #requirements:
        #only genes on normal chromosomes (ch1-22 + X/Y)
        #genes that are not completely contained in another gene
        #genes that are not overlapping with another gene on the same strand (largest one retained)
        print("getting standard genes")
        genesToKeep = getStandardGenes(args.gtfFile, args.outDir)

        #get genes that are not overlapping with another gene on the opposite strand (only look at exons)
        print("getting non-overlapping genes")
        nonOverlappingGenes = getNonOverlappingGenes(args.gtfFile, genesToKeep)

        #get genes with only 1 transcript
        if not args.noGeneFilt:
            print("getting genes with only 1 transcript")
            uniqueTranscriptGenes = getUniqueTranscriptGenes(args.gtfFile, nonOverlappingGenes)
        else:
            uniqueTranscriptGenes = nonOverlappingGenes

        #get GC content and exons for each selected gene
        print("getting GC content of genes")
        exons, seqGC, geneGC, geneLen = getGCcontent(args.gtfFile, args.fastaFile, uniqueTranscriptGenes)

        #write intermediary files
        with open(args.outDir + "/exons.bin", 'wb') as f:
            pickle.dump(exons, f)
        with open(args.outDir + "/seqGC.bin", 'wb') as f:
            pickle.dump(seqGC, f)
        with open(args.outDir + "/geneGC.bin", 'wb') as f:
            pickle.dump(geneGC, f)
        with open(args.outDir + "/geneLen.bin", 'wb') as f:
            pickle.dump(geneLen, f)
        with open(args.outDir + "/uniqueTranscriptGenes.bin", 'wb') as f:
            pickle.dump(uniqueTranscriptGenes, f)
        args.stage = 'readBams'

    else:
        #read intermediary files
        with open(args.outDir + "/exons.bin" , 'rb') as f:
            exons = pickle.load(f)
        with open(args.outDir + "/seqGC.bin" , 'rb') as f:
            seqGC = pickle.load(f)
        with open(args.outDir + "/geneGC.bin" , 'rb') as f:
            geneGC = pickle.load(f)
        with open(args.outDir + "/geneLen.bin" , 'rb') as f:
            geneLen = pickle.load(f)
        with open(args.outDir + "/uniqueTranscriptGenes.bin", 'rb') as f:
            uniqueTranscriptGenes = pickle.load(f)

    #read alignment files
    if (args.stage == 'readBams'):
        
        if (args.bamDir is not None):
            #read all bamfiles
            bamFiles = []
            for file in listdir(args.bamDir):
                if re.search(".*\\.bam$", file):
                    bamFiles.append(args.bamDir  +"/" + file)
            #select nBam random files
            if (args.nBam is not None):
                bamFiles = random.sample(bamFiles, args.nBam)
            
        else:
            bamFiles = [args.bamFile]
            
        readBams(exons, seqGC, geneGC, geneLen, bamFiles, args.numReads, args.outDir, args.nCPU, uniqueTranscriptGenes)
        args.stage = 'plotStats'
        
    #plot figures
    if (args.stage == 'plotStats'):
        plotBias(args.outDir, args.phred)
    
if __name__ == "__main__":
    main()
















