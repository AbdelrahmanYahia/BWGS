args=$@ 
echo -e "Args: ${args}"

re='^[0-9]+$'
if ! [[ $args =~ $re ]] ; then
   echo "Input is string" 

else 
    echo "input is Number"
fi

esearch -db genome \
    -query "Agrobacterium tumefaciens [ORGN]"| efetch \
    -format docsum | xtract \
    -pattern DocumentSummary -element Id

esearch -db genome \
    -query 177 |elink \
    -target assembly|esummary|xtract \
    -pattern DocumentSummary \
    -element FtpPath_RefSeq | uniq -u -

esearch \
    -db genome -query 177 | elink \
    -target assembly | esummary | xtract \
    -pattern DocumentSummary \
    -element Organism | uniq -u -


Agrobacterium tumefaciens str. B6 (a-proteobacteria)


taxID=$(esearch -db genome \
    -query "Agrobacterium tumefaciens str. B6 (a-proteobacteria)"| efetch \
    -format docsum | xtract \
    -pattern DocumentSummary -element TaxId)

esearch -db genome \
    -query $taxID |elink \
    -target assembly|esummary|xtract \
    -pattern DocumentSummary \
    -element FtpPath_RefSeq | uniq -u -

esearch -db genome \
    -query 177 |elink \
    -target assembly|esummary|xtract \
    -pattern DocumentSummary \
    -element 'Id,AssemblyAccession,Organism,Coverage'


esearch -db taxonomy \
    -query "Agrobacterium tumefaciens [Organism]"|elink \
    -target nuccore|efilter \
    -query "refseq"|esummary|xtract \
    -pattern DocumentSummary \
    -element 'Id,Organism,AssemblyAccession,Coverage'


esearch -db genome \
    -query 177 |elink \
    -target assembly|esummary



esearch -db genome \
    -query 177 | efilter \
    -query "refseq" |elink \
    -target assembly|esummary|xtract \
    -pattern DocumentSummary \
    -element 'Id,Organism,AssemblyAccession,Coverage,Stats'



esearch -db genome -query 177 | elink -target assembly | efilter -query "refseq"


echo -e "Organism\tAssemblyAccession\tTaxid\tassembly-status\tcontig_count\tcontig_l50\tcontig_n50\ttotal_length" > complete-assembly-info.tsv
esearch -db assembly -query '"Agrobacterium tumefaciens" [ORGN] AND "latest refseq"[filter]' | esummary | xtract -pattern DocumentSummary \
        -def "NA" -element "Organism,AssemblyAccession,Taxid,assembly-status" -block Stat \
        -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
        -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
         -element Stat \
        >> complete-assembly-info.tsv


esearch -db assembly -query GCF_001541315.1 | elink -target nucleotide -name \
        assembly_nuccore_refseq | efetch -format fasta > GCF_001541315.1.fa


esearch -db genome -query "'taxid358' [ORGN]" | esummary | xtract -pattern DocumentSummary \
                -def NA -element Organism,AssemblyAccession,Taxid,assembly-status -block Stat \
                -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
                -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
                -element Stat

esearch -db genome -query "txid358" |elink -target assembly|esummary| xtract -pattern DocumentSummary \
                -def NA -element Organism,AssemblyAccession,Taxid,assembly-status -block Stat \
                -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
                -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
                -element Stat

tax=358

esearch -db genome -query "txid${tax}" |elink -target assembly|esummary| xtract -pattern DocumentSummary \
                -def NA -element Organism,AssemblyAccession,Taxid,assembly-status -block Stat \
                -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
                -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
                -element Stat