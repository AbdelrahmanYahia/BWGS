import argparse
import pandas as pd
import json

parser = argparse.ArgumentParser(usage='python3 extract_plasmid_info.py --input platon_out.jason --output out.csv')
parser.add_argument('-i', '--input', help= "jason file", metavar='path', type=str)  
parser.add_argument('-o', '--output', help= "outfile", metavar='path', type=str)
args = parser.parse_args()

# reading the JSON data using json.load()
file = args.input
with open(file) as JSON_file:
    dict_JSON_file = json.load(JSON_file)

found_plamids = {}
id = ""

for i in dict_JSON_file:
    if len(dict_JSON_file[i]["plasmid_hits"]) == 0:
        print(f"{i} didn't contain any plasmids")
    else:
        id = str(i)
        plasmid_id = (dict_JSON_file[id]["plasmid_hits"][0]["plasmid"]["id"])
        df = dict_JSON_file[id]["plasmid_hits"][0]
        df["plasmid"] = plasmid_id
        found_plamids[id] = df

plasmid_info = pd.DataFrame(found_plamids)
plasmid_info = plasmid_info.T
plasmid_info.index.name = 'contig_id'
plasmid_info.reset_index(inplace=True)
fileout = args.output
plasmid_info.to_csv(fileout, index=False)

