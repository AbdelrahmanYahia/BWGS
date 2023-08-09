acc=$1

esearch -db assembly -query $acc </dev/null \
    | esummary \
    | xtract -pattern DocumentSummary -element FtpPath_GenBank \
    | while read -r url ; do
        fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
        echo "\033[;33;1m\ntrying to download ${fname}\033[;39;m\n" ;
        wget "$url/$fname" ;
    done ;
done
