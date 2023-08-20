from curses import meta
import os 
import argparse
import subprocess
import pandas as pd 

parser = argparse.ArgumentParser(description='NCBI download reference genome for taxid of tax name',
                                 prog='NCBI_download_ref',
                                 usage='')
parser.version = '1.2'
parser.add_argument('--name', help= "Organism Name",type=str)  
parser.add_argument('--taxid', help= "Taxonomy ID",type=str)  
parser.add_argument('--kraken_report', help= "Kraken report path",type=str)  
parser.add_argument("-o","--output",help= "Output name")

args = parser.parse_args()

def retrive_ref(tax,taxid=False):
    print("\033[;32;1mSearching...\033[;39;m")
    if taxid == False:
        out = subprocess.Popen(f"esearch -db assembly -query '{tax}[ORGN]' | esummary | xtract -pattern DocumentSummary \
                -def NA -element FtpPath_GenBank,Organism,AssemblyAccession,Taxid,assembly-status -block Stat \
                -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
                -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
                -element Stat", shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8')
    else:
        tax=int(tax)
        out = subprocess.Popen(f"esearch -db genome -query 'txid{tax}' |elink -target assembly|esummary| xtract -pattern DocumentSummary \
                -def NA -element FtpPath_GenBank,Organism,AssemblyAccession,Taxid,assembly-status -block Stat \
                -if Stat@category -equals contig_count -or Stat@category -equals contig_l50 \
                -or Stat@category -equals contig_n50 -or Stat@category -equals total_length \
                -element Stat", shell=True, stdout=subprocess.PIPE).communicate()[0].decode('utf-8')

    raw_lst = []
    lst = out.strip().split("\n")
    for i in lst:
        raw_lst.append(i.split("\t"))
    df = pd.DataFrame(raw_lst)
    df.columns =["FtpPath_GenBank",'Organism', 'AssemblyAccession', 
                    'Taxid', 'assembly-status', 'contig_count',
                    'contig_l50', 'contig_n50', 'total_length']

    df["Organism"].unique()
    df = df.astype({"contig_count": int})
    return df


def filter_taxa_df(df):
    print("\033[;33;1mPlease Select an Organism: \033[;39;m")
    counter = 1
    uniques = df["Organism"].unique()
    for i in uniques:
        print(f"{counter}: {i}")
        counter += 1
    n = int(input("\n"))
    m = n-1
    out = df[df["Organism"] == df["Organism"].unique()[m]].sort_values("contig_count")
    
    print("\033[;33;1mPlease Select an Assembly level: \033[;39;m")
    counter = 1
    uniques = out["assembly-status"].unique()
    for i in uniques:
        print(f"{counter}: {i}")
        counter += 1
    n = int(input("\n"))
    m = n-1
    out2 = out[out["assembly-status"] == out["assembly-status"].unique()[m]].sort_values("contig_count")
    return out2

def download_ref(df):
    df2 = df[["FtpPath_GenBank", "AssemblyAccession", "total_length", "assembly-status", "contig_count", "contig_l50", "contig_n50"]]
    df2.index = df2.reset_index(drop = True).index + 1

    if len(df2.index) > 1:
        print("\033[;33;1mPlease Select a Number to Download: \033[;39;m")
        print(df2.loc[:, df2.columns != 'FtpPath_GenBank'].to_string())
        n = int(input("\n"))
        m = n-1
        acc = df2.iloc[m]["AssemblyAccession"]
        ftp = df2.iloc[m]["FtpPath_GenBank"]
    else:
        print(df2.to_string())
        acc = df2.iloc[0]["AssemblyAccession"]
        ftp = df2.iloc[0]["FtpPath_GenBank"]
    ftp_full = ftp + "/" + ftp.split("/")[-1] + "_genomic.fna.gz"
    os.system(f"wget -cO - {ftp_full} > {args.output}")
    print("Done")

def process_kraken(path):
    df = pd.read_table(path)
    df.columns =['perc', 'clade', 'assigned', 'rank', 'txid', 'name']
    out = df[df["rank"] == "S" ].sort_values("perc",ascending=False)
    out.index = out.reset_index(drop = True).index + 1
    name = str(out.iloc[0]["name"]).strip()
    return name

if args.name is None :
    if args.taxid is None:
        if args.kraken_report is None:
            print("\033[;31;1mNo --name or --taxid or --kraken_report supplied: \033[;39;m")
            exit(1)
        else:
            df = retrive_ref(process_kraken(args.kraken_report))
    else:
        df = retrive_ref(args.taxid,taxid=True)
else:
    if args.taxid is None:
        df = retrive_ref(args.name)

df_filterred = filter_taxa_df(df)
download_ref(df_filterred)