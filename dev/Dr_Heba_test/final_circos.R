for_circos <- read.delim("D:/Dr_heba/Prokka/H_EL%2E/for_circos.tsv", header=FALSE)
df <- for_circos[for_circos$V3 != "tRNA", ]


library(stringr)

for ( i in rownames(df)){
test_str <- df[i,"V9"]
df[i,"gene"] <- str_remove(str_extract(test_str, "gene=[\\w.-]+"), "gene=")
df[i,"EC_N"] <- str_remove(str_extract(test_str, "eC_number=[\\w.-]+"), "eC_number=")
}
df2 <- dplyr::filter(df,  !is.na(gene))
df2 <- dplyr::filter(df2,  !is.na(EC_N))

GENES <- unique(df2$gene)

BiocManager::install("org.EcK12.eg.db")
library(clusterProfiler)
EL <- search_kegg_organism('', by='scientific_name')
search_kegg_organism('ec', by='kegg_code')

library(AnnotationDbi)
library(org.EcK12.eg.db)

found_IDs <- AnnotationDbi::select(org.EcK12.eg.db, keys=GENES, columns='ENTREZID', keytype='SYMBOL')
found_IDs <- dplyr::filter(found_IDs,  !is.na(ENTREZID))

ENTERZIDS <- found_IDs$ENTREZID
library(clusterProfiler)

eg2np <- bitr_kegg(ENTERZIDS, fromType='ncbi-geneid', toType='kegg', organism='eco')
keggs <- eg2np$kegg

mkk <- enrichMKEGG(gene = keggs,
                   organism = 'eco',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)
head(mkk)      
mkk




gene.df <- bitr(ENTERZIDS, fromType = "ENTREZID",
                toType = "GENENAME",
                OrgDb = org.EcK12.eg.db)

ego2 <- enrichGO(gene         = gene.df$GENENAME,
                 OrgDb         = org.EcK12.eg.db,
                 keyType       = 'GENENAME',
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)
head(ego2, 3)   

goplot(ego2)


