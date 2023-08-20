library("optparse", quietly = T)

############################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input counts file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="input counts file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  write(paste0("\033[0;", 31, "m","No [INPUT] was supplied","\033[0m"), stderr())
  quit(save = "no", status = 1)
}

plotfig <- function(obj,file) {
  out <- tryCatch(
    {
      message("trying to export figure :D")
      png(file,width = 1500,height = 1000,units = "px",res = 200)
      print(obj)
      dev.off()
    },
    error=function(cond) {
      message("\033[0;Task failed Successfully :D\033[0m")
      message("Here's the original error message if you want to know how you failed:")
      message(cond)
      dev.off()
      
    },
    warning=function(cond) {
      message(paste("Task caused a warning!!!!!!"))
      message(cond)
      return(NULL)
    }
  )    
}

outdir <- opt$output

file <- read.table(opt$input,sep="\t",header = T,fill = T)
genes <- unique(file$gene)
genes <- genes[genes != ""]

library(clusterProfiler)
library(AnnotationDbi)
library(org.EcK12.eg.db)

found_IDs <- AnnotationDbi::select(org.EcK12.eg.db, keys=genes, columns='ENTREZID', keytype='SYMBOL')
found_IDs <- dplyr::filter(found_IDs,  !is.na(ENTREZID))
ENTERZIDS <- found_IDs$ENTREZID

eg2np <- bitr_kegg(ENTERZIDS, fromType='ncbi-geneid', toType='kegg', organism='eco')
keggs <- eg2np$kegg

mkk <- enrichMKEGG(gene = keggs,
                   organism = 'eco',
                   pvalueCutoff = 1,
                   qvalueCutoff = 1)

gene.df <- bitr(ENTERZIDS, fromType = "ENTREZID",
                toType = "GENENAME",
                OrgDb = org.EcK12.eg.db)

ego2 <- enrichGO(gene         = gene.df$GENENAME,
                 OrgDb         = org.EcK12.eg.db,
                 keyType       = 'GENENAME',
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

write.csv(mkk,paste0(outdir,"/mkk.csv"),quote = F,row.names = F)
write.csv(ego2,paste0(outdir,"/ego.csv"),quote = F,row.names = F)

tryCatch(
  {
    message("trying to create an object")
    go <- goplot(ego2)
  },
  error=function(cond) {
    message("Task failed Successfully")
    message("Here's how you failed:")
    message(cond)
  },
  warning=function(cond) {
    message(paste("Task caused a warning!!!!!!"))
    message(cond)
    return(NULL)
  }
)

gobar <- barplot(ego2)
mkkbar <- barplot(mkk)
heatmkk <- heatplot(mkk)
dotmkk <- dotplot(mkk)
heatgo <- heatplot(ego2)
dotgo <- dotplot(ego2)


plotfig(go,paste0(outdir,"/GO.png"))
plotfig(gobar,paste0(outdir,"/GO-bar.png"))
plotfig(mkkbar,paste0(outdir,"/mkk-bar.png"))
# plotfig(heatmkk,paste0(outdir,"/heat-mkk.png"))
plotfig(dotmkk,paste0(outdir,"/dot-mkk.png"))
# plotfig(heatgo,paste0(outdir,"/heat-go.png"))
plotfig(dotgo,paste0(outdir,"/dot-go.png"))


