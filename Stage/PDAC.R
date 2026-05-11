#--------------------------------#
# PROJECT: MARIE M2 PDAC PROJECT #
#        START: 2025/07/16       #
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
workdir <- "C:/Users/PC PEO MAX/Videos/stage/PDAC"
setwd(workdir)

fig.path <- file.path(workdir,"Figures")
res.path <- file.path(workdir,"Results")
data.path <- file.path(workdir,"InputData")
comAnn.path <- file.path(workdir,"Annotation") ## path containning annotation files created by Jack

if (!file.exists(res.path)) { dir.create(res.path) }
if (!file.exists(fig.path)) { dir.create(fig.path) }
if (!file.exists(data.path)) { dir.create(data.path) }
if (!file.exists(comAnn.path)) { dir.create(comAnn.path) }

# load R package
library(sva)
library(limma)
library(ClassDiscovery)
library(clusterProfiler)
library(MCPcounter)
library(tableone)
library(jstable)
library(ComplexHeatmap)
library(ggplot2)
library(ggpubr)
library(gplots)
library(RColorBrewer)
library(survival)
library(survminer)
library(circlize)
library(viridis)
library(readxl)
library(CMScaller)
library(estimate)
library(dplyr)
library(tidyr)
library(fmsb)
library(ggalluvial)
library(MOVICS)
library(stringr)
library(forestplot)
library(GSVA)
library(ConsensusClusterPlus)
library(ComplexHeatmap)
library(survival)
library(irr)

options(scipen = 999)

# customized functions
display.progress = function (index, totalN, breakN=20) {
  
  if ( index %% ceiling(totalN/breakN)  ==0  ) {
    cat(paste(round(index*100/totalN), "% ", sep=""))
  }
  
}    

standarize.fun <- function(indata=NULL, halfwidth=NULL, centerFlag=T, scaleFlag=T) {  
  outdata=t(scale(t(indata), center=centerFlag, scale=scaleFlag))
  if (!is.null(halfwidth)) {
    outdata[outdata>halfwidth]=halfwidth
    outdata[outdata<(-halfwidth)]= -halfwidth
  }
  return(outdata)
}

gmt2list <- function(annofile){
  if (!file.exists(annofile)) {
    stop("There is no such gmt file.")
  }
  
  if (tools::file_ext(annofile) == "xz") {
    annofile <- xzfile(annofile)
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
    close(annofile)
  } else if (tools::file_ext(annofile) == "gmt") {
    x <- scan(annofile, what="", sep="\n", quiet=TRUE)
  } else {
    stop ("Only gmt and gmt.xz are accepted for gmt2list")
  }
  
  y <- strsplit(x, "\t")
  names(y) <- sapply(y, `[[`, 1)
  
  annoList <- lapply(y, `[`, c(-1,-2))
}

countToFpkm <- function(counts, effLen){
  N <- sum(counts)
  exp( log(counts) + log(1e9) - log(effLen) - log(N) )
}

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

# set colors
blue   <- "#5bc0eb"
yellow <- "#fde74c"
green  <- "#9bc53d"
red    <- "#f25f5c"
purple <- "#531f7a"
grey   <- "#8693ab"
orange <- "#fa7921"
white  <- "#f2d7ee"
darkred   <- "#F2042C"
lightred  <- "#FF7FBF"
lightblue <- "#B2EBFF"
darkblue  <- "#1d00ff"
cherry    <- "#700353"
lightgrey <- "#dcddde"
nake <- "#F8C364"
gold <- "#ECE700"
cyan <- "#00B3D0"
sun  <- "#E53435"
peach  <- "#E43889"
violet <- "#89439B"
soil   <- "#EC7D21"
lightgreen <- "#54B642"
darkblue   <- "#21498D"
darkgreen  <- "#009047"
brown      <- "#874118"
seagreen   <- "#008B8A"
jco <- c("#2874C5","#EABF00","#868686","#C6524A","#80A7DE")
jama <- c("#3B4E55","#D69044","#44A0D5","#A94747","#81AF96")
npg <- c("#D94C37","#3AA086","#435388")
npg <- c("#E64B35","#4DBBD5","#00A087","#3C5488","#F39B7F")
heatmap.BlWtRd <- c("#6699CC","white","#FF3C38")
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
heatmap.GrWtRd <- c("#2b2d42","#8d99ae","#edf2f4","#ef233c","#d90429")
heatmap.L.BlYlRd <- c("#4281a4","#9cafb7","#ead2ac","#e6b89c","#fe938c")
heatmap.BlBkRd <- c("#54FEFF","#32ABAA","#125456","#000000","#510000","#A20000","#F30000")
heatmap.BlWtRd2 <- c("#183869","#4195C1","white","#CB5746","#62011D")
mycol <- brewer.pal(n = 12, "Paired")
# add cibersort data
col_end <- c("#D52D2D") # endo
col_fib <- c("#3276A8","#95C2D7") # fib
col_lym <- c("#FF8102","#F7B400","#DFC77F") # lym
#col_epi <- c("#4A9C46","#64A399","#1CBE4F","#29D0D0","#F09193") # epi
col_epi <- c("#4A9C46","#64A399","#1CBE4F") # epi
col_mye <- c("#674696","#CE6DBD","#332288")# mye
col_seu <- c(col_epi, col_fib, col_end, col_lym, col_mye)

# geneset
MSigDB.HMARK <- read.gmt(file.path(comAnn.path,"h.all.v2025.1.Hs.symbols.gmt"))
MSigDB.ALL <- read.gmt(file.path(comAnn.path,"msigdb.v2025.1.Hs.symbols.gmt"))
hmark.list <- gmt2list(file.path(comAnn.path,"h.all.v2025.1.Hs.symbols.gmt"))

# process expression data of 192 samples from 2 batches
Ginfo44 <- read.delim(file = file.path(comAnn.path,"gencode.v44.annotation.gene.txt"),sep = "\t",row.names = NULL,check.names = F,stringsAsFactors = F,header = T)
rownames(Ginfo44) <- substr(Ginfo44$gene_id,1,15)
df <- read_excel(file.path(data.path,"marie.RNA.2025Sep11.xlsx"), sheet = 1) %>% as.data.frame() 
df <- df[-c(1:3),]
rownames(df) <- df$gene_id
df <- df[,setdiff(colnames(df),c("gene_id","gene_name","gene_type"))]

batchinfo <- read_excel(file.path(data.path,"batchInfo.xlsx"), sheet = 1) %>% as.data.frame() 
rownames(batchinfo) <- batchinfo$sample
df <- df[,batchinfo$sample]


modcombat <- model.matrix(~1, data = batchinfo)
count192 <- as.data.frame(ComBat_seq(counts=as.matrix(df), batch=batchinfo$batch))
count192 <- count192[rowSums(count192) > 0,]

Ginfo44 <- Ginfo44[rownames(count192),]
fpkms <- apply(count192, 2, countToFpkm, effLen = Ginfo44$width)
tpms <- apply(fpkms,2,fpkmToTpm)
tpms <- as.data.frame(round(tpms,4))
tpms <- log2(tpms + 1)
tpms <- apply(tpms, 2, function(x) tapply(x, INDEX=factor(Ginfo44$gene_name), FUN=median, na.rm=TRUE))
###################################################

###############################################

Mids <- Ginfo44[which(Ginfo44$gene_type == "protein_coding"),"gene_name"]

# mcpcounter
mcp <- MCPcounter.estimate(expression = tpms,
                           featuresType = "HUGO_symbols",
                           probesets = read.table(file.path(comAnn.path,"probesets.txt"),sep="\t",stringsAsFactors=FALSE,colClasses="character"),
                           genes = read.table(file.path(comAnn.path,"genes.txt"),sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE))

# load clinical information
sinfo <- read_excel(file.path(data.path,"Base de donneés.xlsx"),sheet = 1) %>% as.data.frame()
sinfo <- sinfo[1:337,]
sinfo$Patient <- paste0("P",sinfo$`Numéro patient dans la cohorte`)
sinfo <- sinfo[which(sinfo$`Code tube Genomeast` != "/."),]

annCol <- sinfo[,c("Localisation tumorale","SEX (0 = M, 1 = F)","AJCC STAGE","ADK CANALAIRE","NEODJUVANT CHEMOTHERAPY","Patient","GRADE OF DIFFERENTIATION","TYPE OF INTERVENTION","FIBROSIS (0-1-2-3)","IMMUNE INFILTRATION (0-1)","% CT","LOCALIZATION")]
rownames(annCol) <- sinfo$`Code tube Genomeast`
all(is.element(colnames(tpms),rownames(annCol)))
#############################################CHEMO####################################################
chemo <- tolower(sinfo$TYPE...120)
names(chemo) <- sinfo$`Code tube Genomeast`

chemo[is.na(chemo) | chemo == ""] <- "none"

chemo_group <- rep("FOLFIRINOX", length(chemo))  # par défaut
names(chemo_group) <- names(chemo)

# FOLFOX uniquement (sans folfoxiri)
chemo_group[grepl("folfox", chemo) & !grepl("folfoxiri", chemo)] <- "FOLFOX"

# None
chemo_group[chemo == "none"] <- "None"

# annotation
annCol$Traitement <- chemo_group[rownames(annCol)]
annCol$Traitement[is.na(annCol$Traitement)] <- "None"
########################################################
# samples with low tumor content to remove
problem.samples <- c(
  "HWLG295","HWLG266","HWLG267","HWLG212","HWLG347",
  "HWLG187","HWLG236","HWLG275","HWLG161","HWLG337",
  "HWLG246","HWLG162","HWLG342","HWLG355","HWLG17",
  "HWLG248","HWLG125","HWLG119","HWLG296",
  "HWLG215","HWLG285","HWLG15","HWLG311","HWLG92",
  "HWLG185","HWLG255","HWLG95"
)

# remove samples from expression matrix
tpms <- tpms[, !colnames(tpms) %in% problem.samples]

# remove samples from annotation
annCol <- annCol[!rownames(annCol) %in% problem.samples, ]

# remove samples from clinical data
sinfo <- sinfo[!sinfo$`Code tube Genomeast` %in% problem.samples, ]
annCol <- annCol[colnames(tpms), ]
####################################################################

# change name of sex
annCol$`SEX (0 = M, 1 = F)` <- ifelse(annCol$`SEX (0 = M, 1 = F)` == "1","Female","Male")
colnames(annCol)[2] <- "Sex"

# change category of stage
table(annCol$`AJCC STAGE`)
annCol[which(annCol$`AJCC STAGE` %in% c("0","NA")),"AJCC STAGE"] <- "N/A"
annCol[which(annCol$`AJCC STAGE` %in% c("I","I (car par de T1c dans la classification AJCC en tant que tel)","IA","IB")),"AJCC STAGE"] <- "I"
annCol[which(annCol$`AJCC STAGE` %in% c("IIA","IIB","IIB/MSI-high")),"AJCC STAGE"] <- "II"
colnames(annCol)[3] <- "Stage"

# set N/A for missing
annCol[is.na(annCol)] <- "N/A"

# change names
colnames(annCol)[4] <- "ADK"
colnames(annCol)[7] <- "Differentiation"
colnames(annCol)[8] <- "Type of intervention"
colnames(annCol)[9] <- "FIBROSIS"
colnames(annCol)[10] <- "IMMUNE INFILTRATION"
colnames(annCol)[11] <- "CT"
colnames(annCol)[12] <- "Localisation"

#annCol$`INFILTRATION DUODELE` <- as.character(annCol$`INFILTRATION DUODELE`)
annCol$`IMMUNE INFILTRATION` <- as.character(annCol$`IMMUNE INFILTRATION`)
annCol$FIBROSIS <- as.character(annCol$FIBROSIS)
annCol$CT <- as.numeric(annCol$CT)
annCol[is.na(annCol$CT),"CT"] <- 0

# predict previous subtypes for ADK only
sam.adk <- intersect(rownames(annCol[which(annCol$ADK == "1"),]), colnames(tpms))

# health samples
sam.health <- sinfo[which(sinfo$`Localisation tumorale` == "Tissu sain"),]
annCol[sam.health$`Code tube Genomeast`,"Type of intervention"] <- "N/A"
annCol[sam.health$`Code tube Genomeast`,"Localisation"] <- "N/A"
annCol[sam.health$`Code tube Genomeast`,"Differentiation"] <- "N/A"
#annCol[sam.health$`Code tube Genomeast`,"INFILTRATION DUODELE"] <- "N/A"

sam.health <- intersect(sam.health$`Code tube Genomeast`, colnames(tpms))
annCol <- annCol[c(sam.adk, sam.health),]

annCol$Differentiation <- case_when(
  # Bien différencié group
  grepl("Bien", annCol$Differentiation, ignore.case = TRUE) &
    !grepl("Peu", annCol$Differentiation, ignore.case = TRUE) &
    !grepl("Indiff", annCol$Differentiation, ignore.case = TRUE) ~ "Bien différencié",
  
  # Modérément différencié group (various spellings)
  grepl("Mod", annCol$Differentiation, ignore.case = TRUE) &
    !grepl("Peu", annCol$Differentiation, ignore.case = TRUE) ~ "Modérément différencié",
  
  # Peu différencié group
  grepl("Peu", annCol$Differentiation, ignore.case = TRUE) ~ "Peu différencié",
  
  # Indifférencié/dédifférencié group
  grepl("indiff|dediff", annCol$Differentiation, ignore.case = TRUE) ~ "Indifférencié",
  
  # N/A or missing
  is.na(annCol$Differentiation) | annCol$Differentiation %in% c("NA", "N/A", "") ~ "N/A",
  
  TRUE ~ "N/A"  # fallback for any unexpected term
)

annCol$Localisation <- case_when(
  grepl("Tête|Uncus|Isthme", annCol$Localisation, ignore.case = TRUE) ~ "Head",
  grepl("Corps", annCol$Localisation, ignore.case = TRUE) & 
    !grepl("queue", annCol$Localisation, ignore.case = TRUE) ~ "Body",
  grepl("Queue|Corps et queue|Cors et queue", annCol$Localisation, ignore.case = TRUE) ~ "Tail",
  grepl("Tête, corps et queue", annCol$Localisation, ignore.case = TRUE) ~ "Diffuse",
  is.na(annCol$Localisation) | annCol$Localisation %in% c("NA", "N/A", "") ~ "Unknown",
  TRUE ~ "Other"
)

annCol$`Type of intervention` <- case_when(
  grepl("^DPC", annCol$`Type of intervention`, ignore.case = TRUE) ~ "DPC",
  grepl("^DPT", annCol$`Type of intervention`, ignore.case = TRUE) ~ "DPT",
  grepl("^SPG", annCol$`Type of intervention`, ignore.case = TRUE) ~ "SPG",
  is.na(annCol$`Type of intervention`) | annCol$`Type of intervention` %in% c("NA", "N/A", "") ~ "Unknown",
  TRUE ~ "Other"
)
annCol[which(annCol$`Type of intervention` == "Unknown"),"Type of intervention"] <- "N/A"
annCol[which(annCol$Localisation == "Unknown"),"Localisation"] <- "N/A"

# correct fibrosis
annCol[which(annCol$FIBROSIS == "0%"), "FIBROSIS"] <- "0"
annCol[which(annCol$FIBROSIS == "2%"), "FIBROSIS"] <- "2"
annCol[which(annCol$FIBROSIS == "0.02"), "FIBROSIS"] <- "2"

## COLLISON
template.collison <- read_excel(file.path(data.path,"Base de données pour Xiaofan.xlsx"), sheet = 4) %>% as.data.frame()
colnames(template.collison) <- c("class","probe")
ntp.collison <- ntp(t(scale(t(as.matrix(tpms[,sam.adk])),scale = TRUE, center = TRUE)),
                    templates = template.collison,
                    seed = 20000112,
                    doPlot = TRUE)
write.table(ntp.collison, file = file.path(res.path,"ntp.collison.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

## MOFFIT TUMORAL
template.moffit1 <- read_excel(file.path(data.path,"Base de données pour Xiaofan.xlsx"), sheet = 5, col_names = NA) %>% as.data.frame()
colnames(template.moffit1) <- c("class","probe")
ntp.moffit1 <- ntp(t(scale(t(as.matrix(tpms[,sam.adk])),scale = TRUE, center = TRUE)),
                   templates = template.moffit1,
                   seed = 20000112,
                   doPlot = TRUE)
write.table(ntp.moffit1, file = file.path(res.path,"ntp.moffit1.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

## MOFFIT STROMA
template.moffit2 <- read_excel(file.path(data.path,"Base de données pour Xiaofan.xlsx"), sheet = 6, col_names = NA) %>% as.data.frame()
colnames(template.moffit2) <- c("class","probe")
ntp.moffit2 <- ntp(t(scale(t(as.matrix(tpms[,sam.adk])),scale = TRUE, center = TRUE)),
                   templates = template.moffit2,
                   seed = 20000112,
                   doPlot = TRUE)
write.table(ntp.moffit2, file = file.path(res.path,"ntp.moffit2.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

## BAILEY
template.bailey <- read_excel(file.path(data.path,"Base de données pour Xiaofan.xlsx"), sheet = 7, col_names = NA) %>% as.data.frame()
colnames(template.bailey) <- c("probe","class")
ntp.bailey <- ntp(t(scale(t(as.matrix(tpms[,sam.adk])),scale = TRUE, center = TRUE)),
                  templates = template.bailey,
                  seed = 20000112,
                  doPlot = TRUE)
write.table(ntp.bailey, file = file.path(res.path,"ntp.bailey.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)

# set previous classifications
annCol$COLLISON <- "N/A"
annCol[sam.adk,"COLLISON"] <- as.character(ntp.collison$prediction)

annCol$MOFFIT_TUMORAL <- "N/A"
annCol[sam.adk,"MOFFIT_TUMORAL"] <- as.character(ntp.moffit1$prediction)

annCol$MOFFIT_STROMA <- "N/A"
annCol[sam.adk,"MOFFIT_STROMA"] <- as.character(ntp.moffit2$prediction)

annCol$BAILEY <- "N/A"
annCol[sam.adk,"BAILEY"] <- as.character(ntp.bailey$prediction)

# calculate tumor purity for all samples
indata <- tpms
write.table(indata,file = file.path(res.path,"pdac_hugo.txt"),sep = "\t",row.names = T,col.names = NA,quote = F)
filterCommonGenes(input.f=file.path(res.path, "pdac_hugo.txt") , output.f=file.path(res.path,"pdac_hugo_ESTIMATE.txt"), id="GeneSymbol")
estimateScore(file.path(res.path,"pdac_hugo_ESTIMATE.txt"), file.path(res.path,"pdac_hugo_estimate_score.txt"), platform="affymetrix")
est.pdac <- read.table(file = file.path(res.path,"pdac_hugo_estimate_score.txt"),header = T,row.names = NULL,check.names = F,stringsAsFactors = F,sep = "\t")
rownames(est.pdac) <- est.pdac[,2]; colnames(est.pdac) <- est.pdac[1,]; est.pdac <- est.pdac[-1,c(-1,-2)];
est.pdac <- sapply(est.pdac, as.numeric); rownames(est.pdac) <- c("StromalScore","ImmuneScore","ESTIMATEScore","TumorPurity"); est.pdac.backup = as.data.frame(est.pdac); colnames(est.pdac.backup) <- colnames(indata)
est.pdac <- as.data.frame(t(est.pdac))
rownames(est.pdac) <- colnames(tpms)

annCol$immune_score <- est.pdac[rownames(annCol), "ImmuneScore"]
annCol$stromal_score <- est.pdac[rownames(annCol), "StromalScore"]
annCol$tumor_purity <- est.pdac[rownames(annCol), "TumorPurity"]
###############################################################
annTrackScale <- function(x, halfwidth = 3) {
  x <- scale(x)
  x[x > halfwidth] <- halfwidth
  x[x < -halfwidth] <- -halfwidth
  return(as.numeric(x))
}
######################################################""
tmp <- annCol$immune_score; names(tmp) <- rownames(annCol)
annCol$immune_score <- annTrackScale(tmp, halfwidth = 3)

tmp <- annCol$stromal_score; names(tmp) <- rownames(annCol)
annCol$stromal_score <- annTrackScale(tmp, halfwidth = 3)

# set colors
annColors <- list()
annColors[["Localisation tumorale"]] <- c("Tissu sain" = "grey95","Tumeur centre" = "black","Tumeur périphérique" = brown)
annColors[["Sex"]] <- c("Female" = red, "Male" = blue, "N/A" = "grey95")
annColors[["Stage"]] <- c("I" = darkblue, "II" = green, "III" = red, "IV" = cherry,"N/A" = "grey95")
annColors[["ADK"]] <- c("1" = "#653C90", "0" = "#A2C9D9", "N/A" = "grey95")
annColors[["NEODJUVANT CHEMOTHERAPY"]] <- c("1" = "#653C90", "0" = "#A2C9D9", "N/A" = "grey95")
annColors[["Traitement"]] <- c("FOLFIRINOX" = "#E41A1C", "FOLFOX" = "#377EB8", "None" = "grey95" )
annColors[["COLLISON"]] <- c("Classical_Collisson" = npg[1], "Exocrine" = npg[2] ,"QM" = npg[3], "N/A" = "grey95")
annColors[["MOFFIT_TUMORAL"]] <- c("Basal_Like" = jco[1], "Classical_Moffitt" = jco[2], "N/A" = "grey95")
annColors[["MOFFIT_STROMA"]] <- c("Activated" = cherry, "Normal" = nake, "N/A" = "grey95")
annColors[["BAILEY"]] <- c("ADEX" = red, "Immunogenic" = blue, "Progenitor" = yellow, "Squamous" = green, "N/A" = "grey95")
annColors[["Differentiation"]] <- c("Bien différencié" = darkblue, "Modérément différencié" = green, "Peu différencié" = red, "Indifférencié" = cherry, "N/A" = "grey95")
annColors[["Localisation"]] <- c("Body" = darkblue, "Head" = green, "Tail" = red, "N/A" = "grey95")
annColors[["Type of intervention"]] <- c("DPC" = npg[1], "DPT" = npg[2] ,"SPG" = npg[3], "N/A" = "grey95")
#annColors[["INFILTRATION DUODELE"]] <- c("1" = "#653C90", "0" = "#A2C9D9", "N/A" = "grey95")
annColors[["FIBROSIS"]] <- c("0" = darkblue, "1" = green, "2" = red, "3" = cherry,"N/A" = "grey95")
annColors[["IMMUNE INFILTRATION"]] <- c("1" = "#653C90", "0" = "#A2C9D9", "N/A" = "grey95")
#annColors[["Health tissue"]] <- c("ADM" = mycol[2], "ADM and PaNin" = mycol[4] ,"Chronic pancreatitis" = mycol[6], "Duodenum" = mycol[8], "Healthy pancreas" = mycol[10], "TIPMP and PaNin" = mycol[12], "N/A" = "grey95")
annColors[["immune_score"]] <- bluered(64)
annColors[["stromal_score"]] <- bluered(64)
annColors[["tumor_purity"]] <- NMF:::ccRamp(c("grey90","black"),64)
annColors[["CT"]] <- NMF:::ccRamp(c("grey90","black"),64)
annColors[["Subtype"]] <- c("C1" = jco[1],"C2" = jco[2], "C3" = jco[4])

#annColors[["CYP_group"]] <- c("High_CYP" = jco[1],"Low_CYP" = jco[2], "N/A" = "grey95")

# unsupervised clustering of all 159 samples including ADK and normal -----------------
indata <- tpms[Mids,rownames(annCol)]
indata <- indata[apply(indata,1,function(x) sum(x > 1) > 0.1*ncol(indata)),]
var.sel <- apply(indata, 1, mad)
#var.sel <- names(var.sel[var.sel > quantile(var.sel, probs = seq(0,1,0.25))[4]])
var.sel <- names(var.sel[var.sel > quantile(var.sel, probs = seq(0,1,0.1))[10]])
indata <- indata[var.sel,]
indata <- t(scale(t(indata), center = T, scale = T))
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D2")
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward.D")
plotdata <- standarize.fun(indata, halfwidth = 3)
hm1 <- pheatmap(plotdata,
                use_raster = T,
                border_color = NA,
                show_rownames = F,
                show_colnames = T,
                cluster_rows = hcg,
                cluster_cols = hcs,
                annotation_col = annCol[,c("Localisation tumorale","NEODJUVANT CHEMOTHERAPY","Traitement","Sex","Stage","ADK","COLLISON","MOFFIT_TUMORAL","MOFFIT_STROMA","BAILEY","Differentiation","Localisation","Type of intervention","FIBROSIS","IMMUNE INFILTRATION","immune_score","stromal_score","tumor_purity","CT")],
                #annotation_row = annRow.regulon[rownames(regulon_activity),"family",drop = F],
                annotation_colors = annColors,
                #cellwidth = 400/ncol(plotdata),
                cellwidth = 8,
                cellheight = 400/nrow(plotdata),
                name = "Expr.",
                color = NMF:::ccRamp(heatmap.BlWtRd2,64))
pdf(file = file.path(fig.path,"unsupervised heatmap of high variable genes in 159 samples of marie for qc.pdf"), width = 25,height = 12)
draw(hm1, heatmap_legend_side = "left", annotation_legend_side = "left")
invisible(dev.off())

# unsupervised clustering of all ADK samples without normal -----------------
annCol.ADK <- subset(annCol, ADK == "1")
indata <- tpms[Mids,rownames(annCol.ADK)]
indata <- indata[apply(indata,1,function(x) sum(x > 1) > 0.1*ncol(indata)),]
var.sel <- apply(indata, 1, mad)
#var.sel <- names(var.sel[var.sel > quantile(var.sel, probs = seq(0,1,0.25))[4]])
var.sel <- names(var.sel[var.sel > quantile(var.sel, probs = seq(0,1,0.1))[10]])
indata <- indata[var.sel,]
indata <- t(scale(t(indata), center = T, scale = T))
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D")
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward.D")
plotdata <- standarize.fun(indata, halfwidth = 3)
hm1 <- pheatmap(plotdata,
                use_raster = F,
                border_color = NA,
                show_rownames = F,
                show_colnames = F,
                cluster_rows = hcg,
                cluster_cols = hcs,
                annotation_col = annCol.ADK[,c("Localisation tumorale","NEODJUVANT CHEMOTHERAPY","Traitement","Sex","Stage","ADK","COLLISON","MOFFIT_TUMORAL","MOFFIT_STROMA","BAILEY","Differentiation","Localisation","Type of intervention","FIBROSIS","IMMUNE INFILTRATION","immune_score","stromal_score","tumor_purity","CT")],
                #annotation_row = annRow.regulon[rownames(regulon_activity),"family",drop = F],
                annotation_colors = annColors,
                cellwidth = 400/ncol(plotdata),
                cellheight = 400/nrow(plotdata),
                name = "Expr.",
                color = NMF:::ccRamp(heatmap.BlWtRd2,64))
pdf(file = file.path(fig.path,"unsupervised heatmap of high variable genes in 125 ADK samples of marie.pdf"), width = 15,height = 12)
draw(hm1, heatmap_legend_side = "left", annotation_legend_side = "left")
invisible(dev.off())

# unsupervised clustering of all ADK samples without normal in tumor center -----------------
annCol.ADK.TC <- subset(annCol, ADK == "1" & `Localisation tumorale` == "Tumeur centre")
indata <- tpms[Mids,rownames(annCol.ADK.TC)]
indata <- indata[apply(indata,1,function(x) sum(x > 1) > 0.1*ncol(indata)),]
var.sel <- apply(indata, 1, mad)
#var.sel <- names(var.sel[var.sel > quantile(var.sel, probs = seq(0,1,0.25))[4]])
var.sel <- names(var.sel[var.sel > quantile(var.sel, probs = seq(0,1,0.1))[10]])
indata <- indata[var.sel,]
indata <- t(scale(t(indata), center = T, scale = T))
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D")
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward.D")
plotdata <- standarize.fun(indata, halfwidth = 3)
hm1 <- pheatmap(plotdata,
                use_raster = F,
                border_color = NA,
                show_rownames = F,
                show_colnames = T,
                cluster_rows = hcg,
                cluster_cols = hcs,
                annotation_col = annCol.ADK.TC[,c("Localisation tumorale","NEODJUVANT CHEMOTHERAPY","Traitement","Sex","Stage","ADK","COLLISON","MOFFIT_TUMORAL","MOFFIT_STROMA","BAILEY","Differentiation","Localisation","Type of intervention","FIBROSIS","IMMUNE INFILTRATION","immune_score","stromal_score","tumor_purity","CT")],
                #annotation_row = annRow.regulon[rownames(regulon_activity),"family",drop = F],
                annotation_colors = annColors,
                cellwidth = 10,
                cellheight = 400/nrow(plotdata),
                name = "Expr.",
                color = NMF:::ccRamp(heatmap.BlWtRd2,64))
pdf(file = file.path(fig.path,"unsupervised heatmap of high variable genes in 75 ADK center samples of marie.pdf"), width = 25,height = 12)
draw(hm1, heatmap_legend_side = "left", annotation_legend_side = "left")
invisible(dev.off())

# unsupervised clustering of all ADK samples without normal in tumor periphery -----------------
annCol.ADK.TP <- subset(annCol, ADK == "1" & `Localisation tumorale` == "Tumeur périphérique")
indata <- tpms[Mids,rownames(annCol.ADK.TP)]
indata <- indata[apply(indata,1,function(x) sum(x > 1) > 0.1*ncol(indata)),]
var.sel <- apply(indata, 1, mad)
#var.sel <- names(var.sel[var.sel > quantile(var.sel, probs = seq(0,1,0.25))[4]])
var.sel <- names(var.sel[var.sel > quantile(var.sel, probs = seq(0,1,0.1))[10]])
indata <- indata[var.sel,]
indata <- t(scale(t(indata), center = T, scale = T))
hcs <- hclust(distanceMatrix(as.matrix(indata), "euclidean"), "ward.D")
hcg <- hclust(distanceMatrix(as.matrix(t(indata)), "euclidean"), "ward.D")
plotdata <- standarize.fun(indata, halfwidth = 3)
hm1 <- pheatmap(plotdata,
                use_raster = F,
                border_color = NA,
                show_rownames = F,
                show_colnames = T,
                cluster_rows = hcg,
                cluster_cols = hcs,
                annotation_col = annCol.ADK.TP[,c("Localisation tumorale","NEODJUVANT CHEMOTHERAPY","Traitement","Sex","Stage","ADK","COLLISON","MOFFIT_TUMORAL","MOFFIT_STROMA","BAILEY","Differentiation","Localisation","Type of intervention","FIBROSIS","IMMUNE INFILTRATION","immune_score","stromal_score","tumor_purity","CT")],
                #annotation_row = annRow.regulon[rownames(regulon_activity),"family",drop = F],
                annotation_colors = annColors,
                cellwidth = 10,
                cellheight = 400/nrow(plotdata),
                name = "Expr.",
                color = NMF:::ccRamp(heatmap.BlWtRd2,64))
pdf(file = file.path(fig.path,"unsupervised heatmap of high variable genes in 50 ADK periphery samples of marie.pdf"), width = 15,height = 12)
draw(hm1, heatmap_legend_side = "left", annotation_legend_side = "left")
invisible(dev.off())

# define survival data
surv.data <- sinfo[,c("Code tube Genomeast","DATE OF INTERVENTION","DATE OF DEATH","DATE OF LAST NEWS","REASON","PROGRESSION M+","DATE")]

# Overall Survival
# Step 1: copy the column to work on
death_raw <- surv.data$`DATE OF DEATH`

# Step 2: identify numeric-like Excel dates and convert them
death_date <- suppressWarnings(as.numeric(death_raw))
death_date <- ifelse(!is.na(death_date),
                     as.character(as.Date(death_date, origin = "1899-12-30")),
                     as.character(death_raw))

# Step 3: remove blanks and other non-standard values
death_date[death_date %in% c("", "NA", "N/A")] <- NA

# Step 4: parse everything that looks like a date into POSIXct
surv.data$DATE_OF_DEATH <- as.POSIXct(death_date, tz = "UTC", tryFormats = c(
  "%Y-%m-%d",
  "%Y/%m/%d",
  "%d/%m/%Y",
  "%m/%d/%Y"
))

# Step 5: same for DATE OF INTERVENTION (if not already done)
surv.data$DATE_OF_INTERVENTION <- as.POSIXct(surv.data$`DATE OF INTERVENTION`,
                                             tz = "UTC",
                                             tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%Y-%m-%d %Z"))

# Step 6: compute OS time in months
surv.data$OS.time <- as.numeric(difftime(surv.data$DATE_OF_DEATH,
                                         surv.data$DATE_OF_INTERVENTION,
                                         units = "days")) / 30
surv.data$OS <- ifelse(!is.na(surv.data$DATE_OF_DEATH), 1, 0)

# progression free survival
convert_mixed_date <- function(x) {
  # remove trailing "UTC" and spaces
  x <- trimws(gsub(" UTC", "", x, fixed = TRUE))
  
  # convert numeric Excel serials to ISO dates
  suppressWarnings({
    num <- suppressWarnings(as.numeric(x))
    x <- ifelse(!is.na(num),
                as.character(as.Date(num, origin = "1899-12-30")),
                as.character(x))
  })
  
  # replace blanks / text NAs
  x[x %in% c("", "NA", "N/A")] <- NA
  
  # parse into POSIXct with flexible formats
  as.POSIXct(x, tz = "UTC", tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%d/%m/%Y"))
}

surv.data$DATE_OF_INTERVENTION <- convert_mixed_date(surv.data$`DATE OF INTERVENTION`)
surv.data$DATE                 <- convert_mixed_date(surv.data$`DATE`)
surv.data$DATE_OF_DEATH        <- convert_mixed_date(surv.data$`DATE OF DEATH`)
surv.data$DATE_OF_LAST_NEWS    <- convert_mixed_date(surv.data$`DATE OF LAST NEWS`)

# initialize output columns
surv.data$PFS.time <- NA_real_
surv.data$PFS      <- NA_integer_

# compute durations (in months)
prog_time  <- as.numeric(difftime(surv.data$DATE, surv.data$DATE_OF_INTERVENTION, units = "days")) / 30
death_time <- as.numeric(difftime(surv.data$DATE_OF_DEATH, surv.data$DATE_OF_INTERVENTION, units = "days")) / 30
last_time  <- as.numeric(difftime(surv.data$DATE_OF_LAST_NEWS, surv.data$DATE_OF_INTERVENTION, units = "days")) / 30

# determine the earliest available event (progression or death)
event_time <- pmin(prog_time, death_time, na.rm = TRUE)

# if both are NA, use last follow-up
event_time[is.infinite(event_time)] <- NA  # pmin returns Inf when both NA
event_time[is.na(event_time) & !is.na(last_time)] <- last_time[is.na(event_time) & !is.na(last_time)]

# assign event indicator: 1 = event (progression or death), 0 = censored
event_status <- ifelse(!is.na(prog_time) | !is.na(death_time), 1,
                       ifelse(!is.na(last_time), 0, NA))

# store results
surv.data$PFS.time <- event_time
surv.data$PFS      <- event_status

surv.data[which(surv.data$PFS.time == 0),"PFS"] <- NA # we dont know the progression actually
surv.data[which(surv.data$PFS.time == 0),"PFS.time"] <- NA
rownames(surv.data) <- surv.data$`Code tube Genomeast`

write.table(surv.data, file = file.path(res.path,"survival data of 335 samples.txt"),sep = "\t",row.names = F,col.names = T,quote = F)


# save image
save.image(file = file.path(workdir,"PDAC.RData"))
#######################################################
#################### TEST ####################################

# ==============================
# 1. Séparer TC / TP
# ==============================
TC <- annCol.ADK[annCol.ADK$`Localisation tumorale` == "Tumeur centre", ]
TP <- annCol.ADK[annCol.ADK$`Localisation tumorale` == "Tumeur périphérique", ]

# ==============================
# 2. VARIABLES
# ==============================
vars <- c("Stage","Differentiation","BAILEY","MOFFIT_TUMORAL",
          "MOFFIT_STROMA","COLLISON","Sex","Traitement",
          "FIBROSIS","IMMUNE INFILTRATION",
          "Localisation","Type of intervention")

cont_vars <- c("immune_score","stromal_score","tumor_purity","CT")

# ==============================
# 3. TESTS CATEGORIELS
# ==============================

# TC
results_cat_TC <- do.call(rbind, lapply(vars, function(v) {
  
  df <- TC[TC[[v]] != "N/A", ]
  tab <- table(df$Cluster, df[[v]])
  
  if (all(dim(tab) > 1)) {
    p <- fisher.test(tab)$p.value
    
    data.frame(
      variable = v,
      p.value = round(p, 3),
      Location = "TC"
    )
  }
}))

# TP
results_cat_TP <- do.call(rbind, lapply(vars, function(v) {
  
  df <- TP[TP[[v]] != "N/A", ]
  tab <- table(df$Cluster, df[[v]])
  
  if (all(dim(tab) > 1)) {
    p <- fisher.test(tab)$p.value
    
    data.frame(
      variable = v,
      p.value = round(p, 3),
      Location = "TP"
    )
  }
}))

# ==============================
# 4. TESTS CONTINUS
# ==============================

# TC
results_cont_TC <- do.call(rbind, lapply(cont_vars, function(v) {
  
  df <- TC[!is.na(TC[[v]]), ]
  
  p <- kruskal.test(df[[v]] ~ df$Cluster)$p.value
  
  data.frame(
    variable = v,
    p.value = round(p, 3),
    Location = "TC"
  )
}))

# TP
results_cont_TP <- do.call(rbind, lapply(cont_vars, function(v) {
  
  df <- TP[!is.na(TP[[v]]), ]
  
  p <- kruskal.test(df[[v]] ~ df$Cluster)$p.value
  
  data.frame(
    variable = v,
    p.value = round(p, 3),
    Location = "TP"
  )
}))

# ==============================
# 5. COMBINER TOUS LES RESULTATS
# ==============================

results_all <- rbind(
  results_cat_TC,
  results_cat_TP,
  results_cont_TC,
  results_cont_TP
)

results_all$p.value <- ifelse(results_all$p.value == 0, "<0.001", results_all$p.value)

# ==============================
# 6. SAUVEGARDE
# ==============================

write.table(results_all,
            file = file.path(res.path, "STAT_RESULTS_TC_TP.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# ==============================
# 7. AFFICHAGE
# ==============================

results_all
##################################################
library(openxlsx)

wb <- createWorkbook()

# ==============================
# TC TABLES
# ==============================
for (v in vars) {
  
  tab <- table(TC$Cluster, TC[[v]])
  tab_df <- as.data.frame.matrix(tab)
  
  # 👉 ajouter colonne Cluster
  tab_df$Cluster <- rownames(tab_df)
  tab_df <- tab_df[, c("Cluster", setdiff(colnames(tab_df), "Cluster"))]
  
  addWorksheet(wb, paste0("TC_", v))
  writeData(wb, sheet = paste0("TC_", v), tab_df)
}

# ==============================
# TP TABLES
# ==============================
for (v in vars) {
  
  tab <- table(TP$Cluster, TP[[v]])
  tab_df <- as.data.frame.matrix(tab)
  
  # 👉 ajouter colonne Cluster
  tab_df$Cluster <- rownames(tab_df)
  tab_df <- tab_df[, c("Cluster", setdiff(colnames(tab_df), "Cluster"))]
  
  addWorksheet(wb, paste0("TP_", v))
  writeData(wb, sheet = paste0("TP_", v), tab_df)
}

# ==============================
# P-VALUES
# ==============================

# optionnel : rendre joli
results_all$p.value <- as.character(results_all$p.value)

addWorksheet(wb, "P_values")
writeData(wb, "P_values", results_all)

# ==============================
# SAVE
# ==============================
saveWorkbook(wb,
             file = file.path(res.path, "STAT_RESULTS.xlsx"),
             overwrite = TRUE)

########################kappa######################"
# ==============================
# 1. Sélection des infos utiles
# ==============================
df <- annCol.ADK[, c("Patient", "Cluster", "Localisation tumorale")]

# ==============================
# 2. Format large (TC vs TP)
# ==============================
df_wide <- df %>%
  pivot_wider(
    names_from = `Localisation tumorale`,
    values_from = Cluster
  )

# ==============================
# 3. Renommer colonnes
# ==============================
colnames(df_wide)[colnames(df_wide) == "Tumeur centre"] <- "TC"
colnames(df_wide)[colnames(df_wide) == "Tumeur périphérique"] <- "TP"

# ==============================
# 4. Nettoyage (cas multiples)
# ==============================
df_wide$TC <- sapply(df_wide$TC, function(x) if(length(x) > 0) x[1] else NA)
df_wide$TP <- sapply(df_wide$TP, function(x) if(length(x) > 0) x[1] else NA)

# ==============================
# 5. Supprimer NA
# ==============================
df_wide <- df_wide %>%
  filter(!is.na(TC) & !is.na(TP))

# ==============================
# 6. IMPORTANT → format facteur
# ==============================
df_wide$TC <- factor(df_wide$TC)
df_wide$TP <- factor(df_wide$TP)

# ==============================
# 7. Calcul du Cohen's kappa
# ==============================
kappa_result <- kappa2(df_wide[, c("TC", "TP")])

print(kappa_result)

# ==============================
# 8. Tableau de concordance
# ==============================
table_TC_TP <- table(df_wide$TC, df_wide$TP)
print(table_TC_TP)
print(kappa_result)

################SURVIE################"""
# ==============================
# SURVIE (TC uniquement)
# ==============================

# garder seulement TC
annCol.TC <- annCol.ADK[annCol.ADK$`Localisation tumorale` == "Tumeur centre", ]

# créer dataframe survie
df.surv <- surv.data
df.surv$Cluster <- annCol.TC[rownames(df.surv), "Cluster"]

# enlever NA
df.surv <- df.surv[!is.na(df.surv$Cluster), ]

# ==============================
# OS
# ==============================

fit.OS <- survfit(Surv(OS.time, OS) ~ Cluster, data = df.surv)

p.OS <- 1 - pchisq(survdiff(Surv(OS.time, OS) ~ Cluster, data = df.surv)$chisq, 2)

ggsurvplot(fit.OS,
           pval = TRUE,
           title = paste0("OS (TC only, p = ", 
                          ifelse(p.OS < 0.001, "<0.001", round(p.OS,3)), ")"))

# ==============================
# PFS
# ==============================

fit.PFS <- survfit(Surv(PFS.time, PFS) ~ Cluster, data = df.surv)

p.PFS <- 1 - pchisq(survdiff(Surv(PFS.time, PFS) ~ Cluster, data = df.surv)$chisq, 2)

ggsurvplot(fit.PFS,
           pval = TRUE,
           title = paste0("PFS (TC only, p = ", 
                          ifelse(p.PFS < 0.001, "<0.001", round(p.PFS,3)), ")"))
###############################DEG######################################
###############################
###############################
# DEG ANALYSIS (LIMMA)
###############################

# ==============================
# DEG ANALYSIS (CLUSTERS)
# ==============================

library(limma)
options(scipen = 999)

# ==============================
# 1. METADATA
# ==============================
pd_cluster <- data.frame(
  Samples = rownames(annCol.ADK),
  Group = annCol.ADK$Cluster
)

# renommer clusters (C1, C2, C3)
pd_cluster$Group <- factor(pd_cluster$Group)
levels(pd_cluster$Group) <- paste0("C", levels(pd_cluster$Group))

# ==============================
# 2. DESIGN MATRIX
# ==============================
design_cluster <- model.matrix(~0 + Group, data = pd_cluster)
colnames(design_cluster) <- levels(pd_cluster$Group)

# ==============================
# 3. EXPRESSION MATRIX
# ==============================
gset_cluster <- tpms[Mids, pd_cluster$Samples]

# ==============================
# 4. MODELE
# ==============================
fit_cluster <- lmFit(gset_cluster, design_cluster)

# ==============================
# 5. CONTRASTS
# ==============================
contrast_cluster <- makeContrasts(
  C1vsC2 = C1 - C2,
  C1vsC3 = C1 - C3,
  C2vsC3 = C2 - C3,
  levels = design_cluster
)

fit_cluster2 <- contrasts.fit(fit_cluster, contrast_cluster)
fit_cluster2 <- eBayes(fit_cluster2)

# ==============================
# 6. RESULTATS
# ==============================
res_C1vsC2 <- topTable(fit_cluster2, coef = "C1vsC2", number = Inf)
res_C1vsC3 <- topTable(fit_cluster2, coef = "C1vsC3", number = Inf)
res_C2vsC3 <- topTable(fit_cluster2, coef = "C2vsC3", number = Inf)

# ==============================
# 7. FORMAT
# ==============================
format_res <- function(res) {
  res$P.Value   <- signif(res$P.Value, 3)
  res$adj.P.Val <- signif(res$adj.P.Val, 3)
  res$logFC     <- round(res$logFC, 3)
  return(res)
}

res_C1vsC2 <- format_res(res_C1vsC2)
res_C1vsC3 <- format_res(res_C1vsC3)
res_C2vsC3 <- format_res(res_C2vsC3)

# ==============================
# 8. FILTRE FDR
# ==============================
res_C1vsC2_sig <- res_C1vsC2[res_C1vsC2$adj.P.Val < 0.05, ]
res_C1vsC3_sig <- res_C1vsC3[res_C1vsC3$adj.P.Val < 0.05, ]
res_C2vsC3_sig <- res_C2vsC3[res_C2vsC3$adj.P.Val < 0.05, ]

# ==============================
# 9. SAVE
# ==============================
write.table(res_C1vsC2, file = file.path(res.path,"DEG_C1vsC2_ALL.txt"), sep="\t", quote=F)
write.table(res_C1vsC3, file = file.path(res.path,"DEG_C1vsC3_ALL.txt"), sep="\t", quote=F)
write.table(res_C2vsC3, file = file.path(res.path,"DEG_C2vsC3_ALL.txt"), sep="\t", quote=F)

write.table(res_C1vsC2_sig, file = file.path(res.path,"DEG_C1vsC2_FDR.txt"), sep="\t", quote=F)
write.table(res_C1vsC3_sig, file = file.path(res.path,"DEG_C1vsC3_FDR.txt"), sep="\t", quote=F)
write.table(res_C2vsC3_sig, file = file.path(res.path,"DEG_C2vsC3_FDR.txt"), sep="\t", quote=F)
###################################################
####################DEG TPvsTC##########################
###############################
###############################
# DEG TC vs TP (PAIRED)
###############################

library(limma)
options(scipen = 999)

# ==============================
# 1. GARDER PATIENTS AVEC TC + TP
# ==============================

patients_TC <- annCol.ADK$Patient[annCol.ADK$`Localisation tumorale` == "Tumeur centre"]
patients_TP <- annCol.ADK$Patient[annCol.ADK$`Localisation tumorale` == "Tumeur périphérique"]

patients_valid <- intersect(patients_TC, patients_TP)

annCol.paired <- annCol.ADK[annCol.ADK$Patient %in% patients_valid, ]

# ==============================
# 2. METADATA
# ==============================

pd <- data.frame(
  Samples = rownames(annCol.paired),
  Patient = annCol.paired$Patient,
  Location = annCol.paired$`Localisation tumorale`
)

# garder seulement TC/TP
pd <- pd[pd$Location %in% c("Tumeur centre", "Tumeur périphérique"), ]

# ==============================
# 3. EXPRESSION
# ==============================

gset <- tpms[Mids, pd$Samples]

# enlever gènes non variables
keep <- apply(gset, 1, var) > 0
gset <- gset[keep, ]

# ==============================
# 4. DESIGN (IMPORTANT)
# ==============================

pd$Location <- factor(pd$Location)
pd$Patient  <- factor(pd$Patient)

design <- model.matrix(~ Patient + Location, data = pd)

# ==============================
# 5. MODELE
# ==============================

fit <- lmFit(gset, design)

fit2 <- eBayes(fit)

# ==============================
# 6. RESULTATS
# ==============================

res_TCvsTP <- topTable(fit2, coef = "LocationTumeur périphérique", number = Inf)

# format
res_TCvsTP$logFC <- round(res_TCvsTP$logFC, 3)
res_TCvsTP$P.Value <- signif(res_TCvsTP$P.Value, 3)
res_TCvsTP$adj.P.Val <- signif(res_TCvsTP$adj.P.Val, 3)

# ==============================
# 7. FILTRE
# ==============================

res_TCvsTP_sig <- res_TCvsTP[res_TCvsTP$adj.P.Val < 0.05, ]

#####################################################################
###PATIENT PAR CLUSTER###
library(dplyr)
library(tidyr)

df <- annCol.ADK[, c("Patient", "Cluster", "Localisation tumorale")]

df_wide <- df %>%
  pivot_wider(
    names_from = `Localisation tumorale`,
    values_from = Cluster
  )

# renommer colonnes
colnames(df_wide)[colnames(df_wide) == "Tumeur centre"] <- "TC"
colnames(df_wide)[colnames(df_wide) == "Tumeur périphérique"] <- "TP"

# par patient
df_wide$TC <- sapply(df_wide$TC, function(x) if(length(x) > 0) x[1] else NA)
df_wide$TP <- sapply(df_wide$TP, function(x) if(length(x) > 0) x[1] else NA)

#TC + TP
df_wide <- df_wide %>%
  filter(!is.na(TC) & !is.na(TP))

# convertir en character
df_wide$TC <- as.character(df_wide$TC)
df_wide$TP <- as.character(df_wide$TP)

df_wide
###############################
geneList <- res_C1vsC3$logFC
names(geneList) <- rownames(res_C1vsC3)

geneList <- geneList[!is.na(geneList)]
geneList <- geneList[!duplicated(names(geneList))]
geneList <- sort(geneList, decreasing = TRUE)
########################
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

run_gsea <- function(res, name){
  
  # ==============================
  # 1. PREPARER geneList
  # ==============================
  geneList <- res$logFC
  names(geneList) <- rownames(res)
  
  geneList <- geneList[!is.na(geneList)]
  geneList <- geneList[!duplicated(names(geneList))]
  geneList <- sort(geneList, decreasing = TRUE)
  
  # ==============================
  # 2. GSEA
  # ==============================
  gsea <- GSEA(
    geneList = geneList,
    TERM2GENE = MSigDB.HMARK,
    pvalueCutoff = 0.05
  )
  
  gsea_res <- as.data.frame(gsea)
  gsea_sig <- gsea_res[gsea_res$p.adjust < 0.05, ]
  
  # ==============================
  # 3. SAVE TABLE
  # ==============================
  write.table(gsea_res,
              file = file.path(res.path, paste0("GSEA_", name, ".txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # ==============================
  # 4. DOTPLOT (résumé global)
  # ==============================
  pdf(file = file.path(fig.path, paste0("GSEA_dotplot_", name, ".pdf")),
      width = 10, height = 8)
  
  print(
    dotplot(gsea, showCategory = 15) +
      ggtitle(paste("GSEA", name))
  )
  
  dev.off()
  
  # ==============================
  # 5. BARPLOT 
  # ==============================
  top <- gsea_sig[order(gsea_sig$NES), ]
  top <- rbind(head(top, 10), tail(top, 10))
  
  pdf(file = file.path(fig.path, paste0("GSEA_barplot_", name, ".pdf")),
      width = 10, height = 8)
  
  print(
    ggplot(top, aes(x = NES, y = reorder(Description, NES))) +
      geom_col(aes(fill = NES)) +
      theme_minimal() +
      labs(title = paste("Top pathways", name),
           x = "NES",
           y = "")
  )
  
  dev.off()
  
  # ==============================
  # 6. ENRICHMENT PLOT
  # ==============================
  if(nrow(gsea_sig) > 0){
    
    pathway <- gsea_sig$ID[1]
    
    pdf(file = file.path(fig.path, paste0("GSEA_enrichment_", name, ".pdf")),
        width = 8, height = 6)
    
    print(gseaplot2(gsea, geneSetID = pathway,
                    title = pathway))
    
    dev.off()
  }
  
  return(gsea)
}
gsea_C1vsC2 <- run_gsea(res_C1vsC2, "C1vsC2")
gsea_C1vsC3 <- run_gsea(res_C1vsC3, "C1vsC3")
gsea_C2vsC3 <- run_gsea(res_C2vsC3, "C2vsC3")
gsea_TCvsTP <- run_gsea(res_TCvsTP, "TCvsTP")
###############
#####################
# COX MODEL
###############################
df.surv$Stage <- annCol.ADK[rownames(df.surv), "Stage"]

# remplacer N/A par NA
df.surv$Stage[df.surv$Stage == "N/A"] <- NA

# convertir en facteur
df.surv$Stage <- factor(df.surv$Stage)
df.cox <- df.surv[!is.na(df.surv$Stage), ]
library(survival)

cox <- coxph(Surv(OS.time, OS) ~ Cluster + Stage, data = df.cox)

summary(cox)
library(survminer)

ggforest(cox, data = df.cox)
ggsurvplot(fit.OS, pval = TRUE)
ggadjustedcurves(cox, data = df.cox, variable = "Cluster")

############################################
###########################################
library(dplyr)
library(tidyr)
library(irr)
library(pheatmap)

# 1. Sélection des colonnes utiles
df <- annCol.ADK[, c("Patient", "Cluster", "Localisation tumorale")]

# 2. Transformation en format large
df_wide <- df %>%
  pivot_wider(
    names_from = `Localisation tumorale`,
    values_from = Cluster
  )

# 3. Vérifier les noms des colonnes
colnames(df_wide)

# 4. Renommer proprement (IMPORTANT : pas à l’aveugle)
colnames(df_wide)[colnames(df_wide) == "Tumeur centre"] <- "TC"
colnames(df_wide)[colnames(df_wide) == "Tumeur périphérique"] <- "TFI"

# 5. Gérer les duplications (listes → prendre 1 valeur)
df_wide$TC  <- sapply(df_wide$TC, function(x) if(length(x) > 0) x[1] else NA)
df_wide$TFI <- sapply(df_wide$TFI, function(x) if(length(x) > 0) x[1] else NA)

# 6. Supprimer les valeurs manquantes
df_wide <- df_wide %>%
  filter(!is.na(TC) & !is.na(TFI))

# 7. Convertir en character (format propre pour kappa)
df_wide$TC  <- as.character(df_wide$TC)
df_wide$TFI <- as.character(df_wide$TFI)

# 8. Vérification
head(df_wide)
table(df_wide$TC, df_wide$TFI)

# 9. Calcul du Cohen’s kappa
kappa_result <- kappa2(df_wide[, c("TC", "TFI")])
print(kappa_result)

# 10. Visualisation
mat <- table(df_wide$TC, df_wide$TFI)

pheatmap(prop.table(mat, 1),
         display_numbers = TRUE,
         main = "Concordance TC vs TFI")


################################
####################
fit.OS <- survfit(Surv(OS.time, OS) ~ 1, data = surv.data)

ggsurvplot(
  fit.OS,
  data = surv.data,
  risk.table = TRUE,
  title = "Overall Survival",
  xlab = "Time (months)",
  ylab = "Survival probability"
)
summary(fit.OS)$table
fit.PFS <- survfit(Surv(PFS.time, PFS) ~ 1, data = surv.data)

ggsurvplot(
  fit.PFS,
  data = surv.data,
  risk.table = TRUE,
  title = "Progression Free Survival",
  xlab = "Time (months)",
  ylab = "Survival probability"
)

ggsurvplot(
  fit.OS,
  data = surv.data,
  risk.table = TRUE,
  surv.median.line = "hv",
  title = "Overall Survival",
  xlab = "Time (months)",
  ylab = "Survival probability"
)



##################################################"
library(pheatmap)

# 1. sélectionner top gènes
top_genes <- rownames(res_C1vsC2_sig)[1:50]  # top 50

# 2. matrice expression
mat <- tpms[top_genes, pd$Samples]

# 3. normalisation (très important)
mat_scaled <- t(scale(t(mat)))

# 4. annotation des colonnes (clusters)
annotation_col <- data.frame(
  Cluster = pd$Group
)
rownames(annotation_col) <- pd$Samples

# 5. heatmap
pheatmap(mat_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 DEG - C1 vs C2")
pdf(file = file.path(fig.path, "Heatmap_DEG_C1vsC2.pdf"),
    width = 8, height = 10)

pheatmap(mat_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         main = "Top 50 DEG - C1 vs C2")

dev.off()






