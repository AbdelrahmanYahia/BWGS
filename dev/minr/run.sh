# for file in *.gz; do if [[ "$file" == *"R1"* ]];then \
# f1=$file; elif [[ "$file" == *"R2"* ]];then f2=$file; \
# fi; if [ $((counter%2)) -eq 0 ];then trimmer_pe $f1 $f2 ; \
# else :; fi; counter=$(( counter + 1 )) ; done

# bwa mem -t 60  "../../DR_Amany/new/BWA_Index/GRCh38_latest_genomic.fna" \
# 1_S1_L001_R1_001.trim.fastq.gz 1_S1_L001_R2_001.trim.fastq.gz | \
# samtools view -bS | samtools sort -o 1_S1_L001.trim.sorted.bam -


# for file in SAMPLE_host_removed* ; do if [[ "$file" == *"R1"* ]];then \
# f1=$file; elif [[ "$file" == *"R2"* ]];then f2=$file; \
# fi; if [ $((counter%2)) -eq 0 ];then echo "performing analysis on $f1 and $f2 ... " ;\
# kraken_ $kraken_db 70 ${file%'.fastq.gz'}'.txt' $f1 $f2; \
# else :; fi; counter=$(( counter + 1 )) ; done


source ../main.sh
for file in *fastq.gz
do 
    if [[ "$file" == "SAMPLE"* ]] 
    then
        : 
    else 
        if [[ "$file" == *"R1"* ]]
        then 
            f1=$file
        elif [[ "$file" == *"R2"* ]]
        then 
            f2=$file
        fi
        
        if [ $((counter%2)) -eq 0 ]
        then 
            echo "performing analysis on $f1 and $f2 ... " 
            kraken_ $kraken_db 70 ${file%'.fastq.gz'}'.txt' $f1 $f2
        else
            :
        fi
    fi
    counter=$(( counter + 1 ))
done

kraken2 --db $kraken_db --threads 70 \
            --report "all_unmapped.txt" \
            all_unmapped.fastq > all_unmapped.out 
