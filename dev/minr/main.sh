kraken_db="/media/genomics/Data/DBs/downloaded_kraken2/k2_pluspfp_20210127"
bwa_index="../../DR_Amany/new/BWA_Index/GRCh38_latest_genomic.fna"
hisat2_index="/media/genomics/Data/DBs/Hisat2_human_genome/genome"


quality_control(){
    mkdir -p QC 
    echo -e "fastQC analysis for $file"  
    fastqc -t 70 -f fastq -noextract $file -o QC/
}
multiqc (){
    cd QC/ && multiqc -z -o . . # perform multi QC 
    cd .. 

}

trimmer_pe () { # takes 2 file names, R1 and R2 and counter
    f1=$1
    f2=$2

    logname=${f1%'.fastq.gz'}'.log'
    summaryname=${f1%'.fastq.gz'}'.summry'
    
    new_f1=${f1%'.fastq.gz'}'.trim.fastq.gz'
    newf1=${new_f1}
    newf_2=${f2%'.fastq.gz'}'.trim.fastq.gz'
    newf2=${newf_2}

    newf_1U=${f1%'.fastq.gz'}'.se.trim.fastq.gz'
    newf1U=${newf_1U}
    newf_2U=${f2%'.fastq.gz'}'.se.trim.fastq.gz'
    newf2U=${newf_2U}

    adap="$CONDA_PREFIX/share/trimmomatic-0.39-1/adapters"

    trimmomatic PE -threads 20 -phred33 -trimlog ${logname} \
            -summary ${summaryname}  $f1 $f2 $newf1 $newf1U $newf2 $newf2U \
            SLIDINGWINDOW:4:10 MINLEN:30 \
            ILLUMINACLIP:$adap/TruSeq3-PE-2.fa:2:30:10:1
    ## PE -> paired ended
    ## SLIDINGWINDOW: Performs a sliding window trimming approach.
    ## ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10:2:keepBothReads
}

trimmer () {
    f1=$1

    logname=${f1%'.fq.gz'}'.log'
    summaryname=${f1%'.fq.gz'}'.summry'
    
    new_f1=${f1%'fq.gz'}'trim.fq.gz'
    newf1=${new_f1}
 
    adap="$CONDA_PREFIX/share/trimmomatic-0.39-1/adapters"

    trimmomatic SE -threads 20 -phred33 \
                -trimlog ${logname} \
                -summary ${summaryname}  $f1  $newf1 \
                ILLUMINACLIP:$adap/TruSeq3-SE.fa:2:30:10 \
                SLIDINGWINDOW:4:10 MINLEN:30 

}

runhisat2 () {
    f1=$1
    f2=$2
    hisat2 -p 70 -x $hisat2_index -1 $f1 -2 $f2 | \
            samtools view -bS | \
            samtools sort -o ${f1%'.fastq.gz'}'.sorted.bam' - 
}

bwa_rm_human () {
    index=$1
    R1=$2
    R2=$3
    bwa mem -t 72 $index $R1 $R2 | samtools view -bS | samtools sort -o ${f1%'.fastq.gz'}'.sorted.bam' - 
}

extract_unmapped () {
    samtools view -b -f 12 -F 256 $1 > $2
    samtools sort -n -m 5G -@ 2 $2 -o $2.sorted.bam
    samtools fastq -@ 8 $2.sorted.bam \
            -1 SAMPLE_host_removed_R1.fastq.gz \
            -2 SAMPLE_host_removed_R2.fastq.gz \
            -0 /dev/null -s /dev/null -n
}

kraken_ () {
    dbname=$1
    numofT=$2
    reportname=$3
    R1=$4
    R2=$5
    kraken2 --db $dbname --threads $numofT \
            --report $reportname \
            --paired $R1 $R2 > ${3%'.txt'}'.out'
            
}
omar_extract (){
    samtools index 02_Mapping{}.bam 02_Mapping{}.bam.bai
    samtools flagstat --threads 72 02_Mapping/{}.bam > 02_Mapping/stats/{}.txt
    samtools fastq -1 03_fastq_files{}_1.fq -2 03_fastq_files{}_2.fq -0 /dev/null -s /dev/null -n -f 4 02_Mapping{}.bam
}

kraken_omar () {
    kraken2 --use-mpa-style --memory-mapping \
            --threads 60 --db /ramdisk/k2_eupathdb48_20201113 \
            --report 04_Kraken{}_MetaPhlan_Big.tax \
            --paired  03_fastq_files{}_1.fq 03_fastq_files{}_2.fq > 04_Kraken{}Metaphlan_Big.txt
}

kraken2korona () {
    output=$1
    samples=$2
    ktImportTaxonomy -q 2 -t 3 \
                    -o $output \
                    $samples 
}

