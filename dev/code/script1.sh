YEL='\x1b[1;33m'
GRE='\x1b[1;32m'
RED='\x1b[1;31m'
BLU='\x1b[1;34m'
NC='\e[0m'

# build-ref
# bowtie2-build --threads 35 refs/Enterobacter_mori.fna.gz refs/bowtie2_index/Enterobacter_mori


names=(EM)
Rs=(1 2)
echo -e "${YEL}Analysis started${NC}"
for i in ${names[@]} ; do 
    mkdir -p ${i}
    echo -e "${YEL}fastp ${i} started${NC}"
    fastp --in1 "sample/${i}_R1.fastq" --in2 "sample/${i}_R2.fastq" \
          --out1 "sample/${i}-filttered_R1.fastq" --out2 "sample/${i}-filttered_R2.fastq" \
          --thread 16 -h "${i}-report.html"
    echo -e "${YEL}unicycler ${i} started${NC}"
    unicycler -1 "sample/${i}-filttered_R1.fastq" -2 "sample/${i}-filttered_R2.fastq" \
              -o "${i}-Assambly" \
              --threads 35 --no_pilon --no_correct

    # QC assembly
    quast -o ${i}/quast ${i}-Assambly/assembly.fasta
    
    # align-reads
    bowtie2 --threads 60 -x refs/bowtie2_index/Enterobacter_mori \
            -1 sample/${i}-filttered_R1.fastq \
            -2 sample/${i}-filttered_R2.fastq \
            -S ${i}/${i}-filttered.sam

    echo -e "${YEL}Fix-mate${NC}"
    # fix mates and compres 
    samtools sort -n \
                  -O sam ${i}/${i}-filttered.sam | samtools fixmate \
                  -m -O bam - ${i}/${i}-filttered-fixmate.sam
    
    echo -e "${YEL}convert to sorted bam${NC}"
    # convert to bam file and sort
    samtools sort -O bam ${i}/${i}-filttered-fixmate.sam \
                  -o ${i}/${i}-fixmate.bam
    
    echo -e "${YEL}remove duplicates ${NC}"
    # remove duplicates 
    samtools markdup -r -S ${i}/${i}-fixmate.bam ${i}/${i}-fixmate.sorted.dedup.bam

    echo -e "${YEL}Mapping statistics ${NC}"
    # Mapping statistics 
    samtools flagstat ${i}/${i}-fixmate.sorted.dedup.bam > ${i}/${i}-mapping-stats.txt
    qualimap bamqc -bam ${i}/${i}-fixmate.sorted.dedup.bam -outdir ${i}/${i}-alignment-stats

    echo -e "${YEL}extract mapped reads ${NC}"
    # extract mapped reads 
    samtools view -h -b -f 3 ${i}/${i}-fixmate.sorted.dedup.bam > ${i}/${i}-fixmate.sorted.dedup.concordant.bam

    echo -e "${YEL}unmapped reads ${NC}"
    # unmapped reads 
    samtools view -b -f 4 ${i}/${i}-fixmate.sorted.dedup.bam > ${i}/${i}-fixmate.unmapped.bam
    samtools view -c ${i}/${i}-fixmate.unmapped.bam > ${i}/${i}-N-unmapped.txt
    samtools fastq -1 ${i}/${i}-unmapped_R1.fastq.gz \
                   -2 ${i}/${i}-unmapped_R2.fastq.gz \
                   -0 ${i}/${i}-unmapped_U.fastq.gz ${i}/${i}-fixmate.unmapped.bam

done


# kraken2 
kraken_db="/media/genomics/AlphaFold1/db/Kraken2/k2_pluspfp_20210127"
numofT=50
R1="sample/EM-filttered_R1.fastq"
R2="sample/EM-filttered_R2.fastq"
reportname='EM.report'
R1_u="EM/EM-unmapped_R1.fastq.gz"
R2_u="EM/EM-unmapped_R2.fastq.gz"
reportname_u='EM_u.report'
echo -e "${YEL}Kraken reads ${NC}"
kraken2 --db $kraken_db --threads $numofT \
            --report $reportname \
            --paired $R1 $R2 > 'EM.out'
      
echo -e "${YEL}kraken unmapped reads ${NC}"
kraken2 --db $kraken_db --threads $numofT \
            --report $reportname_u \
            --paired $R1_u $R2_u > 'EM_u.out'

echo -e "${GRE}Analysis DONE...${NC}"

