


#get data from Hagemann-Jensen (E-MTAB-8735 / ERS4261296)
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.diySpike.R1.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.diySpike.I1.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.diySpike.R2.fastq.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-8735/Smartseq3.diySpike.I2.fastq.gz

wget ftp://ftp.sra.ebi.ac.uk/vol1/run/ERR383/ERR3835349/Smartseq3.diySpike.unmapped.ENA.bam
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-8735/Smartseq3.diySpike.readcounts.txt
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-8735/Smartseq3.diySpike.UMIcounts.txt
wget https://www.ebi.ac.uk/biostudies/files/E-MTAB-8735/Smartseq3.diySpike.sample_annotation.txt

#get barcode file from sample annotation
awk 'NR > 1 {print $5}' Smartseq3.diySpike.sample_annotation.txt > barcodes_concat.txt

#cell barcode and UMI assignment + alignement using zUMIs
mkdir zUMIs
~/Programs/zUMIs/zUMIs.sh -c -y zUMIs_config.yaml

#get alignment chr1
for file in zUMIs/zUMIs_output/demultiplexed/*.bam
do
    #get cell name
    cell=${file##*/SS3.}
    cell=${cell%.demx.bam}
    echo $cell
    
    samtools index $file
    samtools view -h -b $file chr1 > ~/Documents/smartSim/example_data/aligned/$cell.bam
    
done
samtools view -h -b zUMIs/SS3.filtered.Aligned.GeneTagged.UBcorrected.sorted.bam chr1 > ~/Documents/smartSim/example_data/aligned/aligned_chr1.bam





