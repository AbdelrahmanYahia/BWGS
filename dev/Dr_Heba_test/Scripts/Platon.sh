S1="/home/abdo/Analysis/WGS/test_plasmids/S/Raw_Assembly/assembly.fasta"
outdir="/home/abdo/Analysis/WGS/test_plasmids/platon_out"
S2="/home/abdo/Analysis/WGS/test_plasmids/S/unmappedAssembly/assembly.fasta"
db="/home/abdo/db/platon/db"

echo -e "Sample One Running"
platon \
    --db ${db} \
    --prefix "Raw_Assembly" \
    --output ${outdir} \
    --threads 32 ${S1}

echo -e "Sample Two Running"
platon \
    --db ${db} \
    --prefix "unmappedAssembly" \
    --output ${outdir} \
    --threads 32 ${S2}
