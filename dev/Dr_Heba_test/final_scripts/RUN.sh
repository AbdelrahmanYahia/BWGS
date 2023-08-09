Sdirs="/home/abdo/Analysis/WGS/test_plasmids/S/DR_Heba_samples"
outdir_platon="/home/abdo/Analysis/WGS/test_plasmids/DR_Heba_samples/platon_out"
outdir_blast="/home/abdo/Analysis/WGS/test_plasmids/DR_Heba_samples/Blast_plus"
db_platon="/home/abdo/db/platon/db"
db_blast="/home/abdo/db/plsdb/plsdb.fna"

mkdir -p $outdir_platon
mkdir -p $outdir_blast
eval "$(conda shell.bash hook)"
conda deactivate

for i in ${Sdirs}/*; do
    sample=$i
    base="$(basename -- $sample)"
	id="${base%".fasta"}"

    echo -e "${id} Platon Running"
    conda activate test
    platon \
        --db ${db_platon} \
        --prefix ${id} \
        --output ${outdir_platon} \
        --threads 32 ${sample}

    echo -e "${id} blast+ Running"
    conda activate GUAP
    blastn \
        -query ${sample} \
        -db ${db_blast} -num_threads 32 \
        -outfmt "6 qseqid sseqid pident qlen slen qseq length mismatch gapopen qstart qend sstart send evalue bitscore" \
        -out "${outdir_blast}/${id}.txt"
done



