from curses import meta
import os 
import argparse
from re import T
import subprocess
import pandas as pd 

parser = argparse.ArgumentParser(description='NCBI download reference genome for taxid of tax name',
                                 prog='NCBI_download_ref',
                                 usage='')
parser.version = '1.2'
parser.add_argument('--kraken_report', help= "Kraken report path",type=str,required=True)  
parser.add_argument("-o","--output",help= "Output directory")

args = parser.parse_args()
os.chdir(args.output)

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
    uniques = df["Organism"].unique()
    m = 0
    out = df[df["Organism"] == df["Organism"].unique()[m]].sort_values("contig_count")
    
    uniques = out["assembly-status"].unique()

    m = 0
    out2 = out[out["assembly-status"] == out["assembly-status"].unique()[m]].sort_values("contig_count")
    return out2

def download_ref(df):
    df2 = df[["FtpPath_GenBank", "Organism", "AssemblyAccession", "total_length", "assembly-status", "contig_count", "contig_l50", "contig_n50"]]
    df2.index = df2.reset_index(drop = True).index + 1
    if len(df2.index) > 1:
        ftp = df2.iloc[0]["FtpPath_GenBank"]
    else:
        ftp = df2.iloc[0]["FtpPath_GenBank"]
    organme = df2.iloc[0]["Organism"]
    ftp_full = ftp + "/" + ftp.split("/")[-1] + "_genomic.fna.gz"
    print("\n" + "\033[;33;1mThe following will be downloaded:\n\033[;37;1m" + str(df2.loc[:1, df2.columns != 'FtpPath_GenBank'].to_string())+ "\033[;39;m\n\n")
    os.system(f"wget -cO - {ftp_full} > '{organme}.fna.gz'")
    print("Done")

def process_kraken(path):
    df = pd.read_table(path)
    df.columns =['perc', 'clade', 'assigned', 'rank', 'txid', 'name']
    out = df[df["rank"] == "S" ].sort_values("perc",ascending=False)
    out.index = out.reset_index(drop = True).index + 1
    name = str(out.iloc[0]["name"]).strip()
    return name


df = retrive_ref(process_kraken(args.kraken_report))
df_filterred = filter_taxa_df(df)
download_ref(df_filterred)