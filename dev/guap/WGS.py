import os
import sys
import argparse
import yaml 
import pandas as pd
from collections import defaultdict
import logging

GUAP_DIR = (os.path.abspath(__file__)).replace("/bin/WGS/WGS.py","")
sys.path.append(f"{GUAP_DIR}/bin/common")

from common_functions import *


parser = argparse.ArgumentParser(description='WGS analysis pipeline',
                                 prog='GUAP WGS',
                                 usage='GUAP WGS [-i|--input PATH] [-o|--output PATH] [-m|--metadata FILE]',
                                 add_help=False)
parser.version = '0.9'
parser.add_argument('-i', '--input', help= "Input directory path", metavar='path', type=str)  
parser.add_argument('-o', '--output', help= "Output directory path", metavar='path', type=str)

parser.add_argument('-t', '--threads', help= "Number of total no. of threads", type=int)
parser.add_argument('--threads-index', help= "Number of threads to use during indexing ref", type=int )
parser.add_argument('--threads-assemble', help= "Number of threads to use per sample assembly", type=int, default = 4)
parser.add_argument('--threads-align', help= "Number of threads to use per sample aliging", type=int, default = 4)
parser.add_argument('--threads-kraken', help= "Number of threads to use per sample aliging", type=int)

parser.add_argument('--index', help= "fasta file to index", type=str, metavar='path')
parser.add_argument('--index-name', help= "name of index to use", type=str)
parser.add_argument('--index-path', help= "path to index files", type=str, metavar='path')
parser.add_argument('--org-name', help= "",type=str)
parser.add_argument('--org-txid', help= "",type=str)
parser.add_argument('--org-id', help= "",type=str)
parser.add_argument('--org-download-direct', help= "Download first found record automaticaly", action='store_true')
parser.add_argument('--org-auto-download', help= "",action='store_true')

parser.add_argument('--kraken2-index', help= "", type=str, metavar='path')

parser.add_argument('--skip-QC', action='store_true', help="Skip  step")
parser.add_argument('--skip-check-contaminant', action='store_true', help="Skip  step")
parser.add_argument('--skip-assembly', action='store_true', help="Skip  step")
parser.add_argument('--skip-assembly-QC', action='store_true', help="Skip  step")
parser.add_argument('--pre-check-contaminant', action='store_true', help="Skip  step")

parser.add_argument('--continue', action='store_true', help="continue analysis when re-run")
parser.add_argument('--overwrite', action='store_true', help="overwrite output dir if exsits")

parser.add_argument("-h", "--help", action='store_true')

args = parser.parse_args()

# check if no arguments 
if not len(sys.argv) > 1:
    print(f"\033[;31;1mERROR: \033[;39;mNo Arguments supplied")
    print(help_message_WGS)
    exit(1)

# check for help
if args.help is True:
    print(help_message_WGS)
    exit(0)

# check n of threads and use all threads if not supplied 
if args.threads is None:
    logging.warning(f"\033[;33;1mYou did not specify any number of threads, I will use all ({all_threads}). \033[;39;m")
    args.threads = all_threads
if args.threads_index is None:
    args.threads_index = args.threads
if args.threads_kraken is None:
    args.threads_kraken = args.threads

if args.threads_assemble > args.threads:
    args.thread_assemble == args.threads

if args.threads_index > args.threads:
    args.thread_index == args.threads

if args.threads_align > args.threads:
    args.thread_align == args.threads

elif args.threads < 4:
    logging.error(f"\033[;31;1mGUAP 16s can't use less than 4 threads, would you please change the number of threads?  \033[;39;m")
    new_threads = int(input("Number of threads: (4 and above, enter any thing else to exit) "))
    try:
        int(new_threads)
        if int(new_threads) < 4:
            logging.error(f"\033[;31;1mValue is smaller than 4. exiting...")
            exit(1)
        else:
            args.threads = new_threads
    except ValueError:
        logging.error(f"\033[;31;1mexiting...")
        exit(1)

outpath = os.path.abspath(args.output)
# create output path if doesn't exsit 
if not os.path.exists(outpath) and not args.overwrite:
    print(f"\033[;33;1mcreating working directory @ {outpath} \033[;39;m")
    os.mkdir(outpath)
elif args.overwrite:
    if not os.path.exists(outpath):
        print(f"\033[;33;1mworking directory {outpath} was not found, will create\033[;39;m")
        os.mkdir(outpath)
    else:
        import shutil
        shutil.rmtree(outpath)
        print(f"\033[;33;1mOverwriting working directory @ {outpath} \033[;39;m")
        os.mkdir(outpath)

# create sample table 
samples = defaultdict(dict)
sample_names = set()
pattern = ""
path = os.path.abspath(args.input)
samples = defaultdict(dict)
allfiles = os.listdir(path)

for file in allfiles:
    if os.path.isfile(path+"/"+file):
        if "fastq" in file or "fq" in file:
            if full_match(file) is False:
                continue
            else:
                sample_name, filename, file_extension, matched, matched_groups = full_match(file)
                if "1" in matched_groups[4]:
                    f2 = file.replace("_R1","_R2")
                    f2 = f2.replace("_1","_2")
                    samples[matched_groups[1]]["Sample ID"] = matched_groups[2]
                    samples[matched_groups[1]]["Sample rest of name"] = matched_groups[3]
                    samples[matched_groups[1]]["Extension"] = file_extension
                    samples[matched_groups[1]]["R1"] = file
                    samples[matched_groups[1]]["R2"] = ""
                    if "_001." in file:
                        samples[matched_groups[1]]["tail"] = "_001"
                    else:
                        samples[matched_groups[1]]["tail"] = ""
                    if f2 in allfiles:
                        samples[matched_groups[1]]["R2"] = f2
                        samples[matched_groups[1]]["PE"] = True
                    else:
                        samples[matched_groups[1]]["R2"] = ""
                        samples[matched_groups[1]]["PE"] = False
                    samples[matched_groups[1]]["path"] = (path)
                    if "R" in matched_groups[4]:
                        samples[matched_groups[1]]["R"] = "R"
                    else:
                        samples[matched_groups[1]]["R"] = ""
                    
                elif "2" in matched_groups[4]:
                    continue
                else:
                    samples[matched_groups[1]]["Sample ID"] = matched_groups[2]
                    samples[matched_groups[1]]["Sample rest of name"] = matched_groups[3]
                    samples[matched_groups[1]]["Extension"] = file_extension
                    samples[matched_groups[1]]["R1"] = file
                    samples[matched_groups[1]]["R2"] = ""
                    samples[matched_groups[1]]["PE"] = False
                    if "R" in matched_groups[4]:
                        samples[matched_groups[1]]["R"] = "R"
                    else:
                        samples[matched_groups[1]]["R"] = ""

                    if "_001." in file:
                        samples[matched_groups[1]]["tail"] = "_001"
                    else:
                        samples[matched_groups[1]]["tail"] = ""
                    samples[matched_groups[1]]["path"] = (path)
        else:
            continue
m_samples = samples
samples= pd.DataFrame(samples).T
samples = samples[["Sample ID", "Sample rest of name", "Extension", "R1", "R2", "tail", "PE", "path", "R"]]
samples = samples.sort_values(by=['Sample ID'])

ext = str(check_extension(samples)[0])
tail = str(check_tail(samples)[0])
PE = bool(check_PE(samples))
R = str(check_R(samples)[0])
compressed = False
EXT = ext

# to perform gunzipping 
if ".gz" in ext:
    compressed = True
    EXT = ext.replace(".gz","")
    
# check for paired end files 
if PE is False:
    logging.error(f"\033[;31;1mSamples not in pairs... exiting...")
    exit(1)

# check if analysis run before and created sample table 
if os.path.exists(outpath+"/"+"samples.tsv"):
    logging.warning(f"\033[;33;1mFound an exsiting sample.tsv file in output directory, will not override.\033[;39;m")
else:
    samples.to_csv(outpath+"/"+"samples.tsv",sep='\t')    

# create config file 
with open('config.yaml', 'w') as yaml_file:
    yaml.safe_dump(vars(args), yaml_file, default_flow_style=False, sort_keys=False)

with open('config.yaml', 'a') as yaml_file:
    yaml_file.writelines(f"path: {path}\n")
    yaml_file.writelines(f"working_dir: {outpath}\n")
    yaml_file.writelines(f"ext: {ext}\n")
    yaml_file.writelines(f"tail: {tail}\n")
    yaml_file.writelines(f"R: {R}\n")
    yaml_file.writelines(f"R1_pattern: _{R}1{tail}{EXT}\n")
    yaml_file.writelines(f"R2_pattern: _{R}2{tail}{EXT}\n")
    yaml_file.writelines(f"compressed: {compressed}\n")
    yaml_file.writelines(f"total_mem: {all_mem}\n")
    yaml_file.writelines(f"GUAP_DIR: {GUAP_DIR}")

