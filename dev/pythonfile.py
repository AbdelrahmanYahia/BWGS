import os
import glob
krakens = glob.glob('/home/genomics/Documents/DR_Amany/new/04_Kraken/Kraken2krona/*.txt')
for kraken in krakens:
    with open(kraken) as rtxt:
        lines = rtxt.readlines()
    with open("{}_without_human.krona".format(kraken[:-4]), "w") as wtxt:
        for line in lines:
            if "(taxid 0)" not in line:
                if "Homo" not in line:
                    tmp = line.split("\t")[1:3]
                    read_id = tmp[0]
                    tax_id = tmp[1]
                    tax_id = tax_id[tax_id.index("taxid")+6:tax_id.index(")")]
                    tmp_w = "{}\t{}\n".format(read_id, tax_id)
                    wtxt.write(tmp_w)
