library(tidyverse)
library(magrittr)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(cowplot)
library(ggplotify)

######### ######### ######### ######### ######### ######### ######### 
######### read in files and prepare datasets #### ######### ######### 
######### ######### ######### ######### ######### ######### ######### 

setwd("W:/corrinne/Cytonuclear/arachis/DEanalysis")


##### ##### ##### ##### ##### ##### 
##### establish functions ### ##### 
##### ##### ##### ##### ##### ##### 

# turn off scientific notation and set seed
options(scipen = 999)
set.seed(8675309)

#### function to merge tsv files from kallisto, report counts #### 
countmerge = function(mypath){
    filenames=list.files(path=mypath,pattern = "*.tsv", recursive=T, full.names=TRUE)
    datalist = lapply(filenames, function(x){read.table(file=x,header=T,colClasses=c(NA,"NULL","NULL",NA,"NULL"))})
    Reduce(function(x,y) {merge(x,y,by="target_id")}, datalist) 
}


#### function to merge tsv files from kallisto, report tpm #### 
tpmmerge = function(mypath){
    filenames=list.files(path=mypath,pattern = "*.tsv", recursive=T, full.names=TRUE)
    datalist = lapply(filenames, function(x){read.table(file=x,header=T,colClasses=c(NA,"NULL","NULL","NULL",NA))})
    Reduce(function(x,y) {merge(x,y,by="target_id")}, datalist) 
}

#### function to perform DESeq ####
ddsMake <- function (countdata, metadata) {
    ddsTmp <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design=~species, tidy = TRUE)
    ddsTmp <- ddsTmp[rowSums(counts(ddsTmp))/nrow(metadata) >= 1,]
    ddsTmp <- estimateSizeFactors(ddsTmp)
    ddsTmp <- DESeq(ddsTmp)
    ddsTmp
}

#### function to make PCA/pheatmap images ####
PCApheat <- function (dds_obj,subgenome,cpmt_y_n) {
    vsd <- vst(dds_obj, blind=FALSE)
    colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
    sampleDists <- dist(t(assay(vsd)))
    sampleDistMatrix <- as.matrix(sampleDists)
    rownames(sampleDistMatrix) <- paste(vsd$species, vsd$rep, sep="-")
    colnames(sampleDistMatrix) <- paste(vsd$species, vsd$rep, sep="-")
    	
    cpmtAll <- intersect(rownames(vsd),cpmt)
    CMvsd <- vsd[cpmtAll,]
    CMsampleDists <- dist(t(assay(CMvsd)))
    CMsampleDistMatrix <- as.matrix(CMsampleDists)
    rownames(CMsampleDistMatrix) <- paste(CMvsd$species, CMvsd$rep, sep="-")
    colnames(CMsampleDistMatrix) <- paste(CMvsd$species, CMvsd$rep, sep="-")
	
    a <- plotPCA(vsd, intgroup="species")	
    c <- as.grob(pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists,clustering_distance_cols=sampleDists,col=colors,silent=T))

    if(cpmt_y_n=="yes") {
        b <- plotPCA(vsd[cpmtAll,], intgroup="species")
        d <- as.grob(pheatmap(CMsampleDistMatrix,clustering_distance_rows=CMsampleDists, clustering_distance_cols=CMsampleDists,col=colors,silent=T))
        myplot <- plot_grid(a,b,c,d)
        myheight=8
	}
	
	if(cpmt_y_n=="no") {
        myplot <- plot_grid(a,c)
        myheight=5
	} 
    
    ggsave(paste0("PCA_pheatmap.",subgenome,".jpg"), device="jpeg", units="in", height=myheight, width=11, plot=myplot)
}

### function to generate contingency tables
### gives the object sigResults_CyMiraCat_direction_ct
### table is row 1= test category, row 2=non-targeted category
### column 1 = number DE, column 2= number not DE (nDE)
contingencyTable <- function (sigResults, origdds, CyMiraCat, direction) {
    if (direction=="both") {
        NTDE <- nrow(sigResults[row.names(sigResults) %in% NT,])
        DE <- nrow(sigResults[row.names(sigResults) %in% CyMiraCat,])
    }

    if (direction=="up") {
        NTDE <- nrow(sigResults[row.names(sigResults) %in% NT & sigResults[2] >0,])
        DE <- nrow(sigResults[row.names(sigResults) %in% CyMiraCat & sigResults[2] >0,])
    }

    if (direction=="down") {
        NTDE <- nrow(sigResults[row.names(sigResults) %in% NT & sigResults$log2FoldChange <0,])
        DE <- nrow(sigResults[row.names(sigResults) %in% CyMiraCat & sigResults$log2FoldChange <0,])
    }

    NTnDE <- length(intersect(row.names(origdds),NT))- NTDE
    nDE <- length(intersect(row.names(origdds),CyMiraCat))- DE

    mat <- cbind(c(DE,NTDE),c(nDE,NTnDE))
    mat
}






##### ##### ##### ##### ##### ##### 
##### read in files ### ##### #####
##### ##### ##### ##### ##### ##### 

### use the above functions to aggregate files into dataframes
countdata <- countmerge("W:/corrinne/Cytonuclear/arachis/DEanalysis/counts")
tpmdata <- tpmmerge("W:/corrinne/Cytonuclear/arachis/DEanalysis/counts")

# Ahyp_1 and Aipadur_5 have poor mapping; remove
countdata <- countdata %>% select(!contains("Ahyp_1")) %>% select(!contains("Aipadur_5"))
tpmdata <- tpmdata %>% select(!contains("Ahyp_1")) %>% select(!contains("Aipadur_5"))

write.table(countdata, file="count.arachis.tbl", quote=F, row.names=F, sep="\t")
write.table(tpmdata, file="tpm.arachis.tbl", quote=F, row.names=F, sep="\t")

### generate metadata framework
metadata <- data.frame(id=names(countdata[,-1]), data.frame(id=names(countdata[,-1])) %>% separate(id, c("species", "rep", NA)))
countdata[,-1] <- round(countdata[,-1],0)

### define cp/mt genes
cp <- countdata$target_id[grepl("cp_",countdata$target_id)]
mt <- countdata$target_id[grepl("mt_",countdata$target_id)]
cpmt <- c(cp,mt)

### generate df and metadata for subgenomes
### Aipa is dad for Ah, mom for Ai
### Adur is mom for Ah, dad for Ai
category <- read.table("gene.pairs", col.names=c("category","dad","mom")) %>% 
     filter(dad %in% countdata$target_id) %>% filter(mom %in% countdata$target_id)
category$both <- paste0(category$dad,"_",category$mom)

fullcategory <- read.table("arachis.expanded.category.longForm", col.names=c("category","CyMira","CyMiraSub","gene","subgenome")) %>% 
     filter(gene %in% countdata$target_id) %>% select(gene,subgenome,category,CyMira,CyMiraSub) %>%
     add_row(gene=cpmt, subgenome="NotApp", category="NotApp",CyMira="NotApp", CyMiraSub="NotApp")

Acountdata <- countdata %>% filter(target_id %in% category$mom) %>% select(!contains("Aipa_"))
Dcountdata <- countdata %>% filter(target_id %in% category$dad) %>% select(contains("target") | contains("Ahyp") | contains("Aipa"))

Ametadata <- filter(metadata, !grepl("Aipa$",species))
Dmetadata <- filter(metadata, !grepl("Adur",species)) 

### get CyMIRA categories and make vectors with names
DI <- scan("Dual-targeted_Interacting.list",what="character")
DN <- scan("Dual-targeted_Non-interacting.list",what="character")
MI <- scan("Mitochondria-targeted_Interacting.list",what="character")
MN <- scan("Mitochondria-targeted_Non-interacting.list",what="character")
NT <- scan("Not-organelle-targeted.list",what="character")  
PI <- scan("Plastid-targeted_Interacting.list",what="character")
PN <- scan("Plastid-targeted_Non-interacting.list",what="character")
catList <- c("DI","DN","MI","MN","PI","PN")
catListcpmt <- c(catList,"cp","mt")


### generate other lists
directionList <- c("both","up","down")
momSubgenomes <- c("AhypD","AipadurD")
dadSubgenomes <- c("AhypI","AipadurI")
sumGenomes <- c("AhypSD","AipadurSD","AhypSI","AipadurSI")
sided <- c('two.sided','greater','less')

#### get the count sum of homoeolog pairs ####
homoeoPairs <- category[,2:3]
Scountdata <- countdata[countdata$target_id %in% cpmt,]
Scountdata[,-1] <- round(Scountdata[,-1],0)

for (i in 1:nrow(homoeoPairs)) {
    Scountdata[nrow(Scountdata)+1,] <- rbind(
        cbind(paste0(homoeoPairs[i,1],"_",homoeoPairs[i,2]),
        countdata[countdata$target_id==homoeoPairs[i,1], -1] + 
        countdata[countdata$target_id==homoeoPairs[i,2], -1]))
}

#### get a ratio of counts maternal/paternal, in that order ####
### paternal must be first column or edit below order
### edit the columns that are selected from countdata for polyploid only
Rcountdata <- countdata[1:2,]
Rcountdata <- Rcountdata[-(1:2),c(1,11:19)] 

for (i in 1:nrow(homoeoPairs)) {
    Rcountdata[nrow(Rcountdata)+1,] <- rbind(
        cbind(paste0(homoeoPairs[i,1],"_",homoeoPairs[i,2]),
        countdata[countdata$target_id==homoeoPairs[i,2], c(11:19)]/ 
        countdata[countdata$target_id==homoeoPairs[i,1], c(11:19)]))
}



##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 
######  can check results for homoeologs partitioned, all samples ##### #####
##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### ##### 

dds <- ddsMake(countdata,metadata)
PCApheat(dds,"both","yes")

resultsNames(dds)
res <- results(dds)
summary(res) 

rld <- rlog(dds)
lmTable <- as.data.frame(assay(rld)) %>% rownames_to_column(var = "gene") %>% full_join(fullcategory)
write.table(lmTable,file="arachis.lm.tbl", quote=F, sep="\t", row.names=F)


##### ##### ##### ##### ##### ##### 
##### Maternal (A) comparison ##### 
##### ##### ##### ##### ##### ##### 

Adds <- ddsMake(Acountdata,Ametadata)
PCApheat(Adds,"A","no")

resultsNames(Adds)
Ares <- results(Adds)
summary(Ares) 

### grab DE results, Atha relative to Asue, no cp/mt
AhypDres <- results(Adds, contrast=c("species","Ahyp","Adur"), alpha=0.05)
AhypDres <- AhypDres[order(AhypDres$padj),]
AhypDsig <- subset(AhypDres, padj<0.05)
write.table(AhypDsig,file="DEGsigResults.AhypvAdur.tbl", quote=F, sep="\t")

AipadurDres <- results(Adds, contrast=c("species","Aipadur","Adur"), alpha=0.05)
AipadurDres <- AipadurDres[order(AipadurDres$padj),]
AipadurDsig <- subset(AipadurDres, padj<0.05)
write.table(AipadurDsig,file="DEGsigResults.AipadurvAdur.tbl", quote=F, sep="\t")


### just a couple more checks, for fun
resLFCAhd <- lfcShrink(Adds,coef="species_Ahyp_vs_Adur", type="apeglm")
plotMA(resLFCAhd, ylim = c(-5, 5))

resLFCAid <- lfcShrink(Adds,coef="species_Aipadur_vs_Adur", type="apeglm")
plotMA(resLFCAid, ylim = c(-5, 5))


Avsd <- vst(Adds, blind=FALSE)
topVarGenesA <- head(order(rowVars(assay(Avsd)), decreasing = TRUE), 20)
matA  <- assay(Avsd)[ topVarGenesA, ]
matA  <- matA - rowMeans(matA)
pheatmap(t(matA))



##### ##### ##### ##### ##### ##### 
##### Paternal (D) comparison ##### 
##### ##### ##### ##### ##### ##### 

Ddds <- ddsMake(Dcountdata,Dmetadata)
PCApheat(Ddds,"D","no")

resultsNames(Ddds)
Dres <- results(Ddds)
summary(Dres) 

### grab DE results, Dt relative to D5, no cp/mt
AhypIres <- results(Ddds, contrast=c("species","Ahyp","Aipa"), alpha=0.05)
AhypIres <- AhypIres[order(AhypIres$padj),]
AhypIsig <- subset(AhypIres, padj<0.05)
write.table(AhypIsig,file="DEGsigResults.AhypvAipa.tbl", quote=F, sep="\t")

AipadurIres <- results(Ddds, contrast=c("species","Aipadur","Aipa"), alpha=0.05)
AipadurIres <- AipadurIres[order(AipadurIres$padj),]
AipadurIsig <- subset(AipadurIres, padj<0.05)
write.table(AipadurIsig,file="DEGsigResults.AipadurvAipa.tbl", quote=F, sep="\t")

### just a couple more checks, for fun
resLFCAhi <- lfcShrink(Ddds,coef="species_Aipa_vs_Ahyp", type="apeglm")
plotMA(resLFCAhi, ylim = c(-5, 5))

resLFCAii <- lfcShrink(Ddds,coef="species_Aipadur_vs_Ahyp", type="apeglm")
plotMA(resLFCAii, ylim = c(-5, 5))

Dvsd <- vst(Ddds, blind=FALSE)
topVarGenesD <- head(order(rowVars(assay(Dvsd)), decreasing = TRUE), 20)
matD  <- assay(Dvsd)[ topVarGenesD, ]
matD  <- matD - rowMeans(matD)
pheatmap(t(matD))


##### ##### ##### ##### ##### ##### #####
##### diploid-polyploid comparison  ##### 
##### ##### ##### ##### ##### ##### ##### 

Sdds <- ddsMake(Scountdata,metadata)
PCApheat(Sdds,"summed","yes")

resultsNames(Sdds)
Sres <- results(Sdds)
summary(Sres)

### grab DE results, Asue relative to Atha/are

AhypSDres <- results(Sdds, contrast=c("species","Ahyp","Adur"), alpha=0.05)
AhypSDres <- AhypSDres[order(AhypSDres$padj),]
AhypSDsig <- subset(AhypSDres, padj<0.05)
write.table(AhypSDsig,file="DEGsigResults.AhypSvAdurS.tbl", quote=F, sep="\t")


AipadurSDres <- results(Sdds, contrast=c("species","Aipadur","Adur"), alpha=0.05)
AipadurSDres <- AipadurSDres[order(AipadurSDres$padj),]
AipadurSDsig <- subset(AipadurSDres, padj<0.05)
write.table(AipadurSDsig,file="DEGsigResults.AipadurSvAdurS.tbl", quote=F, sep="\t")


AhypSIres <- results(Sdds, contrast=c("species","Ahyp","Aipa"), alpha=0.05)
AhypSIres <- AhypSIres[order(AhypSIres$padj),]
AhypSIsig <- subset(AhypSIres, padj<0.05)
write.table(AhypSIsig,file="DEGsigResults.AhypSvAipaS.tbl", quote=F, sep="\t")


AipadurSIres <- results(Sdds, contrast=c("species","Aipadur","Aipa"), alpha=0.05)
AipadurSIres <- AipadurSIres[order(AipadurSIres$padj),]
AipadurSIsig <- subset(AipadurSIres, padj<0.05)
write.table(AipadurSIsig,file="DEGsigResults.AipadurSvAipaS.tbl", quote=F, sep="\t")


### summary(AhypSDres)
#out of 10753 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 614, 5.7%
#LFC < 0 (down)     : 681, 6.3%
#outliers [1]       : 121, 1.1%
#low counts [2]     : 0, 0%
#(mean count < 1)

### summary(AhypSIres)
#out of 10753 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 1884, 18%
#LFC < 0 (down)     : 1725, 16%
#outliers [1]       : 121, 1.1%
#low counts [2]     : 0, 0%
#(mean count < 1)

### summary(AipadurSIres)
#out of 10753 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 1944, 18%
#LFC < 0 (down)     : 1617, 15%
#outliers [1]       : 121, 1.1%
#low counts [2]     : 0, 0%
#(mean count < 1)

### summary(AipadurSDres)
#out of 10753 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 113, 1.1%
#LFC < 0 (down)     : 91, 0.85%
#outliers [1]       : 121, 1.1%
#low counts [2]     : 828, 7.7%


############################################
##### get expression relationships #########
#####      vis-a-vis dominance     #########
############################################

Dipres <- results(Sdds, contrast=c("species","Adur","Aipa"), alpha=0.05)

cpmtcat <- rbind(
	category %>% select(category,both) %>% dplyr::rename(gene=both),
	as.data.frame(cbind(gene=cpmt,category=gsub("_.*","",cpmt)))
)
	
is.gt.Ah <- data.frame(gene=character(),PvM=character(),PvD=character(),MvD=character(), stringsAsFactors=F)

for (gene in Scountdata$target_id) {
	is.gt.Ah[nrow(is.gt.Ah)+1,"gene"] <- gene
	is.gt.Ah[is.gt.Ah$gene==gene,"MvD"] <- as.data.frame(Dipres)[gene,] %>% 
   	 mutate(MvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'mom',
               padj < 0.05 & log2FoldChange < 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(MvD)

	is.gt.Ah[is.gt.Ah$gene==gene,"PvM"] <- as.data.frame(AhypSDres)[gene,] %>% 
   	 mutate(PvM = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'mom',
               TRUE ~ 'error' ) ) %>% select(PvM)

	is.gt.Ah[is.gt.Ah$gene==gene,"PvD"] <- as.data.frame(AhypSIres)[gene,] %>% 
   	 mutate(PvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(PvD)
}

is.gt.Ah <- is.gt.Ah %>% full_join(.,cpmtcat)


dominanceCatsAh <- plyr::count(is.gt.Ah, vars = c('PvM','PvD','MvD')) %>% 
  	mutate(type = case_when(
		MvD == "dad" & PvM == "pp" & PvD == "dad" ~ "category01", 
		MvD == "mom" & PvM == "mom" & PvD == "pp" ~ "category12",
		MvD == "dad" & PvM == "pp" & PvD == "eq" ~ "category02",
		MvD == "mom" & PvM == "mom" & PvD == "eq" ~ "category11",
		MvD == "mom" & PvM == "eq" & PvD == "pp" ~ "category04",
		MvD == "dad" & PvM == "eq" & PvD == "dad" ~ "category09",		
		MvD == "dad" & PvM == "mom" & PvD == "dad" ~ "category03", 
		MvD == "eq" & PvM == "mom" & PvD == "dad" ~ "category07",
		MvD == "mom" & PvM == "mom" & PvD == "dad" ~ "category10",
		MvD == "dad" & PvM == "pp" & PvD == "pp" ~ "category05",
		MvD == "mom" & PvM == "pp" & PvD == "pp" ~ "category06",
		MvD == "eq" & PvM == "pp" & PvD == "pp" ~ "category08",
		MvD == "eq" & PvM == "eq" & PvD == "eq" ~ "all_equal",
		MvD == "error" & PvM == "error" & PvD == "error" ~ "no_info",
            TRUE ~ 'undet' ) )



dominanceCatsCytoAh <- plyr::count(is.gt.Ah, vars = c('PvM','PvD','MvD','category')) %>% 
  	mutate(type = case_when(
		MvD == "dad" & PvM == "pp" & PvD == "dad" ~ "category01", 
		MvD == "mom" & PvM == "mom" & PvD == "pp" ~ "category12",
		MvD == "dad" & PvM == "pp" & PvD == "eq" ~ "category02",
		MvD == "mom" & PvM == "mom" & PvD == "eq" ~ "category11",
		MvD == "mom" & PvM == "eq" & PvD == "pp" ~ "category04",
		MvD == "dad" & PvM == "eq" & PvD == "dad" ~ "category09",		
		MvD == "dad" & PvM == "mom" & PvD == "dad" ~ "category03", 
		MvD == "eq" & PvM == "mom" & PvD == "dad" ~ "category07",
		MvD == "mom" & PvM == "mom" & PvD == "dad" ~ "category10",
		MvD == "dad" & PvM == "pp" & PvD == "pp" ~ "category05",
		MvD == "mom" & PvM == "pp" & PvD == "pp" ~ "category06",
		MvD == "eq" & PvM == "pp" & PvD == "pp" ~ "category08",
		MvD == "eq" & PvM == "eq" & PvD == "eq" ~ "all_equal",
		MvD == "error" & PvM == "error" & PvD == "error" ~ "no_info",
            TRUE ~ 'undet' ) )


dominanceTableAh <- bind_rows(dominanceCatsAh,dominanceCatsCytoAh) %>% 
	select(-c(PvM,PvD,MvD)) %>% 
	filter(type != "undet") %>%
	filter(type != "no_info") %>%
	pivot_wider(names_from=type, values_from=freq) %>%
	tidyr::replace_na(list(category="all")) %>%
	replace(is.na(.), 0) %>%
	select(category,category01, category12, category02, category11, category04, category09, category03, category07, category10, category05, category06, category08, all_equal) 


write.table(dominanceTableAh,file="AhExpressionDominance.tbl", quote=F, row.names=F, sep="\t")



	
is.gt.Aid <- data.frame(gene=character(),PvM=character(),PvD=character(),MvD=character(), stringsAsFactors=F)

for (gene in Scountdata$target_id) {
	is.gt.Aid[nrow(is.gt.Aid)+1,"gene"] <- gene
	is.gt.Aid[is.gt.Aid$gene==gene,"MvD"] <- as.data.frame(Dipres)[gene,] %>% 
   	 mutate(MvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange < 0 ~ 'mom',
               padj < 0.05 & log2FoldChange > 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(MvD)

	is.gt.Aid[is.gt.Aid$gene==gene,"PvM"] <- as.data.frame(AipadurSIres)[gene,] %>% 
   	 mutate(PvM = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'mom',
               TRUE ~ 'error' ) ) %>% select(PvM)

	is.gt.Aid[is.gt.Aid$gene==gene,"PvD"] <- as.data.frame(AipadurSDres)[gene,] %>% 
   	 mutate(PvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(PvD)
}

is.gt.Aid <- is.gt.Aid %>% full_join(.,cpmtcat)


dominanceCatsAid <- plyr::count(is.gt.Aid, vars = c('PvM','PvD','MvD')) %>% 
  	mutate(type = case_when(
		MvD == "dad" & PvM == "pp" & PvD == "dad" ~ "category01", 
		MvD == "mom" & PvM == "mom" & PvD == "pp" ~ "category12",
		MvD == "dad" & PvM == "pp" & PvD == "eq" ~ "category02",
		MvD == "mom" & PvM == "mom" & PvD == "eq" ~ "category11",
		MvD == "mom" & PvM == "eq" & PvD == "pp" ~ "category04",
		MvD == "dad" & PvM == "eq" & PvD == "dad" ~ "category09",		
		MvD == "dad" & PvM == "mom" & PvD == "dad" ~ "category03", 
		MvD == "eq" & PvM == "mom" & PvD == "dad" ~ "category07",
		MvD == "mom" & PvM == "mom" & PvD == "dad" ~ "category10",
		MvD == "dad" & PvM == "pp" & PvD == "pp" ~ "category05",
		MvD == "mom" & PvM == "pp" & PvD == "pp" ~ "category06",
		MvD == "eq" & PvM == "pp" & PvD == "pp" ~ "category08",
		MvD == "eq" & PvM == "eq" & PvD == "eq" ~ "all_equal",
		MvD == "error" & PvM == "error" & PvD == "error" ~ "no_info",
            TRUE ~ 'undet' ) )



dominanceCatsCytoAid <- plyr::count(is.gt.Aid, vars = c('PvM','PvD','MvD','category')) %>% 
  	mutate(type = case_when(
		MvD == "dad" & PvM == "pp" & PvD == "dad" ~ "category01", 
		MvD == "mom" & PvM == "mom" & PvD == "pp" ~ "category12",
		MvD == "dad" & PvM == "pp" & PvD == "eq" ~ "category02",
		MvD == "mom" & PvM == "mom" & PvD == "eq" ~ "category11",
		MvD == "mom" & PvM == "eq" & PvD == "pp" ~ "category04",
		MvD == "dad" & PvM == "eq" & PvD == "dad" ~ "category09",		
		MvD == "dad" & PvM == "mom" & PvD == "dad" ~ "category03", 
		MvD == "eq" & PvM == "mom" & PvD == "dad" ~ "category07",
		MvD == "mom" & PvM == "mom" & PvD == "dad" ~ "category10",
		MvD == "dad" & PvM == "pp" & PvD == "pp" ~ "category05",
		MvD == "mom" & PvM == "pp" & PvD == "pp" ~ "category06",
		MvD == "eq" & PvM == "pp" & PvD == "pp" ~ "category08",
		MvD == "eq" & PvM == "eq" & PvD == "eq" ~ "all_equal",
		MvD == "error" & PvM == "error" & PvD == "error" ~ "no_info",
            TRUE ~ 'undet' ) )


dominanceTableAid <- bind_rows(dominanceCatsAid,dominanceCatsCytoAid) %>% 
	select(-c(PvM,PvD,MvD)) %>% 
	filter(type != "undet") %>%
	filter(type != "no_info") %>%
	pivot_wider(names_from=type, values_from=freq) %>%
	tidyr::replace_na(list(category="all")) %>%
	replace(is.na(.), 0) %>%
	add_column(category12 = NA) %>%
	select(category,category01, category12, category02, category11, category04, category09, category03, category07, category10, category05, category06, category08, all_equal) 


write.table(dominanceTableAid,file="AidExpressionDominance_revised.tbl", quote=F, row.names=F, sep="\t")



### get the significant DE genes each category
## also parsed by up/down regulation


AhypD <- read.table("DEGsigResults.AhypvAdur.tbl")
AipadurD <- read.table("DEGsigResults.AipadurvAdur.tbl")

AhypI <- read.table("DEGsigResults.AhypvAipa.tbl")
AipadurI <- read.table("DEGsigResults.AipadurvAipa.tbl")

AhypSD <- read.table("DEGsigResults.AhypSvAdurS.tbl")
AipadurSD <- read.table("DEGsigResults.AipadurSvAdurS.tbl")

AhypSI <- read.table("DEGsigResults.AhypSvAipaS.tbl")
AipadurSI <- read.table("DEGsigResults.AipadurSvAipaS.tbl")



## do fisher's exact for Mom subgenome
mom <- data.frame(species=character(),category=character(),direction=character(),sided=character(),DE=integer(), nDE=integer(), pvalueDE=double(), stringsAsFactors=F)
for (A in momSubgenomes) {
    for (c in catList) {
        for (d in directionList) {
            x <- contingencyTable(get(A),Adds,get(c),d)            
            for (s in sided) {
                ft <- fisher.test(x, alternative=s)
                mom[nrow(mom)+1,] <- rbind(A,c,d,s,x[1,1],x[1,2],ft$p.value)        
                if (c=="PI" & s=="less") { mom[nrow(mom)+1,] <- rbind(A,"NT",d,"NA",x[2,1],x[2,2],"NA") }
}}}}
write.table(mom,file="AdurvAhorAi.fisherexact.tbl", quote=F, row.names=F, sep="\t")

 
           
## make contingency tables for Dad subgenome
dad <- data.frame(species=character(),category=character(),direction=character(),sided=character(),DE=integer(), nDE=integer(), pvalueDE=double(), stringsAsFactors=F)
for (A in dadSubgenomes) {
    for (c in catList) {
        for (d in directionList) {
            x <- contingencyTable(get(A),Ddds,get(c),d)            
            for (s in sided) {
                ft <- fisher.test(x, alternative=s)
                dad[nrow(dad)+1,] <- rbind(A,c,d,s,x[1,1],x[1,2],ft$p.value)        
                if (c=="PI" & s=="less") { dad[nrow(dad)+1,] <- rbind(A,"NT",d,"NA",x[2,1],x[2,2],"NA") }
}}}}
write.table(dad,file="AipaAhorAi.fisherexact.tbl", quote=F, row.names=F, sep="\t")




## make contingency tables for summed polyploid
summed <- data.frame(species=character(),category=character(),direction=character(),sided=character(),DE=integer(), nDE=integer(), pvalueDE=double(), stringsAsFactors=F)
for (A in sumGenomes) {
    for (c in catListcpmt) {
        for (d in directionList) {
            x <- contingencyTable(get(A),Sdds,get(c),d)            
            for (s in sided) {
                ft <- fisher.test(x, alternative=s)
                summed[nrow(summed)+1,] <- rbind(A,c,d,s,x[1,1],x[1,2],ft$p.value)        
                if (c=="PI" & s=="less") { summed[nrow(summed)+1,] <- rbind(A,"NT",d,"NA",x[2,1],x[2,2],"NA") }
}}}}
write.table(summed,file="AsuevAareorAtha.fisherexact.tbl", quote=F, row.names=F, sep="\t")


pubtable <- rbind(mom,dad,summed) %>% group_by(species, category, direction) %>% 
          slice_min(pvalueDE,n=1) %>% 
          mutate(sided=if_else(pvalueDE <= 0.05, paste0(sided,"_SIGNIF"), paste0(sided,"_notSig"))) %>% 
          select(-pvalueDE)

write.table(pubtable,file="publicationFisherArachis.tbl", row.names=F, col.names=T, sep="\t", quote=F)

############################
##################################
######################################
# add Kruskal-Wallis tests??? maybe not because categories/genes are not independent?
######################################
##################################
############################

categoryP <- dplyr::rename(category,target_id=both) %>% select(target_id,category) 

ratioSum <- Rcountdata %>% rowwise() %>% 
     transmute(target_id,  
     Ahypmean = mean(c(Ahyp_2,Ahyp_3,Ahyp_4,Ahyp_5), na.rm=TRUE),
     Ahypmedian = median(c(Ahyp_2,Ahyp_3,Ahyp_4,Ahyp_5), na.rm=TRUE),
     Aidumean = mean(c(Aipadur_1,Aipadur_2,Aipadur_3,Aipadur_4), na.rm=TRUE),
     Aidumedian = median(c(Aipadur_1,Aipadur_2,Aipadur_3,Aipadur_4), na.rm=TRUE),
     AhypNAcnt= sum(is.na(c(Ahyp_2,Ahyp_3,Ahyp_4,Ahyp_5))),
     AiduNAcnt= sum(is.na(c(Aipadur_1,Aipadur_2,Aipadur_3,Aipadur_4)))) %>%
     left_join(., categoryP, by='target_id') %>% 
     gather("tmp", "value", starts_with('A')) %>% 
     separate(tmp, into=c('species', 'type' ), sep=c(4))

ratioMean <- ratioSum %>% 
#    filter(category != 'Not-organelle-targeted') %>% 
    filter(!is.na(value)) %>% 
    filter(value != 'Inf') %>% 
    filter(value != '-Inf') %>%
    filter(type == 'mean')

ratioMedian <- ratioSum %>% 
#    filter(category != 'Not-organelle-targeted') %>% 
    filter(!is.na(value)) %>% 
    filter(value != 'Inf') %>% 
    filter(value != '-Inf') %>%
    filter(type == 'median')

Mn <- ggplot(ratioMean, aes(x=category,y=value,fill=species)) +  
    geom_hline(yintercept=1) + 
    geom_violin(
    aes(alpha=0.2), trim = FALSE,
    position = position_dodge(0.9) 
    ) +
    geom_boxplot(
    aes(alpha=0.3), width = 0.15,
    position = position_dodge(0.9)
    ) +
    scale_fill_manual(values=c("deepskyblue3", "purple")) +
    ylab("Mean value Mom/Dad") + scale_y_log10(limits=c(0.1,10))

Md <- ggplot(ratioMedian, aes(x=category,y=value,fill=species)) + 
    geom_hline(yintercept=1) + 
    geom_violin(
    aes(alpha=0.2), trim = FALSE,
    position = position_dodge(0.9) 
    ) +
    geom_boxplot(
    aes(alpha=0.3), width = 0.15,
    position = position_dodge(0.9)
    ) +
    scale_fill_manual(values=c("deepskyblue3", "purple")) +
    ylab("Median value Mom/Dad")+ scale_y_log10(limits=c(0.1,10))



MnMd <- plot_grid(Mn,Md,nrow=2)
ggsave(paste0("Mean_median_mat2pat_ratio.jpg"), device="jpeg", units="in", height=8, width=22, plot=MnMd)



##### ##### ##### ##### ##### ##### 
##### Subgenome comparison ## ##### 
##### ##### ##### ##### ##### ##### 

PcountdataH <- category %>% 
	select(dad,mom,both) %>% 
	dplyr::rename(target_id=dad) %>%
     	left_join(., Dcountdata, by='target_id') %>% 
	select(c(mom,both,starts_with('Ahyp'))) %>%
	rename_at(vars(starts_with('Ahyp')), ~ str_replace(., 'Ahyp', 'AhypDad')) %>% 
	dplyr::rename(target_id=mom) %>%
     	left_join(., Acountdata, by='target_id') %>%
	select(c(both,starts_with('Ahyp'))) %>%
	rename_at(vars(starts_with('Ahyp_')), ~ str_replace(., 'Ahyp', 'AhypMom')) %>% 
	column_to_rownames(var = "both") 

PmetadataH <- data.frame(id=names(PcountdataH), data.frame(id=names(PcountdataH)) %>% separate(id, c("species", "rep", NA)))

PddsH <- DESeqDataSetFromMatrix(countData = round(PcountdataH,0), colData=PmetadataH, design = ~species) 

libTotalH <- rep(colSums(countdata[,7:10]),2) %>%
	set_names(names(PcountdataH))
libSizeH <- libTotalH/mean(libTotalH)

sizeFactors(PddsH) <- libSizeH
PddsH <- DESeq(PddsH)

PresH <- results(DESeq(PddsH),contrast=c("species","AhypDad","AhypMom"),alpha=0.05)

PsigH <- subset(PresH[order(PresH$log2FoldChange),], padj<0.05)

write.table(PsigH,file="DEGsigResults.HPatvMat.tbl", quote=F, sep="\t")



PcountdataAid <- category %>% 
	select(dad,mom,both) %>% 
	dplyr::rename(target_id=dad) %>%
     	left_join(., Dcountdata, by='target_id') %>% 
	select(c(mom,both,starts_with('Aipadur'))) %>% 
	rename_at(vars(starts_with('Aipadur')), ~ str_replace(., 'Aipadur', 'AidMom')) %>% 
	dplyr::rename(target_id=mom) %>%
     	left_join(., Acountdata, by='target_id') %>%
	select(c(both,starts_with('Aipadur'),starts_with('Aid'))) %>%
	rename_at(vars(starts_with('Aipadur_')), ~ str_replace(., 'Aipadur', 'AidDad')) %>% 
	column_to_rownames(var = "both") 	

PmetadataAid <- data.frame(id=names(PcountdataAid), data.frame(id=names(PcountdataAid)) %>% separate(id, c("species", "rep", NA)))

PddsAid <- DESeqDataSetFromMatrix(countData = round(PcountdataAid,0), colData=PmetadataAid, design = ~species) 

libTotalAid <- rep(colSums(countdata[,11:14]),2) %>%
	set_names(names(PcountdataAid))
libSizeAid <- libTotalAid/mean(libTotalAid)

sizeFactors(PddsAid) <- libSizeAid
PddsAid <- DESeq(PddsAid)

PresAid <- results(DESeq(PddsAid),contrast=c("species","AidDad","AidMom"),alpha=0.05)

PsigAid <- subset(PresAid[order(PresAid$log2FoldChange),], padj<0.05)

write.table(PsigAid,file="DEGsigResults.AidPatvMat.tbl", quote=F, sep="\t")



Bcountdata <- category %>% 
	select(dad,mom,both) %>% 
	dplyr::rename(target_id=dad) %>%
     	left_join(., Dcountdata, by='target_id') %>% 
	select(c(mom,both,starts_with('Aipa_'))) %>%
	dplyr::rename(target_id=mom) %>%
     	left_join(., Acountdata, by='target_id') %>%
	select(c(both,starts_with(c('Aipa_','Adur')))) %>% 
	dplyr::rename(target_id=both)	

Bmetadata <- data.frame(id=names(Bcountdata[,-1]), data.frame(id=names(Bcountdata[,-1])) %>% separate(id, c("species", "rep", NA)))

Bdds <- ddsMake(Bcountdata,Bmetadata)
Bres <- results(Bdds,contrast=c("species","Aipa","Adur"),alpha=0.05)
Bsig <- subset(Bres[order(Bres$log2FoldChange),], padj<0.05)
write.table(Bsig,file="DEGsigResults.AipavAdur.tbl", quote=F, sep="\t")



options(scipen=10)

ResCompare <- as.data.frame(PresH) %>% 
	rownames_to_column(var='both') %>%
	dplyr::arrange(desc(abs(log2FoldChange))) %>%
	mutate(HypRank=1:nrow(PresH)) %>%
     	right_join(., category, by='both') %>% 
	dplyr::rename(target_id=both, HypPatMatFC=log2FoldChange, HypPatMatpadj=padj) %>%
     	left_join(., (as.data.frame(PresAid) %>% rownames_to_column(var='target_id') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(AidRank=1:nrow(PresAid))), by='target_id') %>% 
	dplyr::rename(AidPatMatFC=log2FoldChange, AidPatMatpadj=padj) %>%
     	left_join(., (as.data.frame(AhypDres) %>% rownames_to_column(var='mom') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(HypMomRank=1:nrow(AhypDres))), by='mom') %>% 
	dplyr::rename(HypMatDurFC=log2FoldChange, HypMatDurpadj=padj) %>%
     	left_join(., (as.data.frame(AipadurDres) %>% rownames_to_column(var='mom') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(AidDadRank=1:nrow(AipadurDres))), by='mom') %>% 
	dplyr::rename(AidPatDurFC=log2FoldChange, AidPatDurpadj=padj) %>%
     	left_join(., (as.data.frame(AhypIres) %>% rownames_to_column(var='dad') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(HypDadRank=1:nrow(AhypIres))), by='dad') %>% 
	dplyr::rename(HypPatIpaFC=log2FoldChange, HypPatIpapadj=padj) %>% 
     	left_join(., (as.data.frame(AipadurIres) %>% rownames_to_column(var='dad') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(AidMomRank=1:nrow(AipadurIres))), by='dad') %>% 
	dplyr::rename(AidMatIpaFC=log2FoldChange, AidMatIpapadj=padj) %>%
     	left_join(., (as.data.frame(Bres) %>% rownames_to_column(var='target_id') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(DipRank=1:nrow(Bres))), by='target_id') %>% 
	dplyr::rename(IpaDurFC=log2FoldChange,IpaDurpadj=padj) %>% 
	select(matches('target_id|category|Hyp|Aid|Dur|Rank'))

write.table(ResCompare,file="DEGResults.ResCompare.arachis.tbl", quote=F, sep="\t", row.names=F)






## make contingency tables for homoeolog comparison
PsigHdf <- as.data.frame(PsigH)

homoeoCompH <- data.frame(category=character(),sided=character(), FC=character(0),
	MatMoreCat=integer(), PatMoreCat=integer(), 
	pvalueDE=double(), stringsAsFactors=F)
for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(PsigHdf[row.names(PsigHdf) %in% get(c) & PsigHdf$log2FoldChange <0,])
		PatMoreCat <- nrow(PsigHdf[row.names(PsigHdf) %in% get(c) & PsigHdf$log2FoldChange >0,])
		MatMoreNT <- nrow(PsigHdf[row.names(PsigHdf) %in% NT & PsigHdf$log2FoldChange <0,])
		PatMoreNT <- nrow(PsigHdf[row.names(PsigHdf) %in% NT & PsigHdf$log2FoldChange >0,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoCompH[nrow(homoeoCompH)+1,] <- rbind(c,s,"any",mat[1,1],mat[2,1],ft$p.value) } }
homoeoCompH[nrow(homoeoCompH)+1,] <- rbind("NT","NA","any",mat[1,2],mat[2,2],"NA") 


for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(PsigHdf[row.names(PsigHdf) %in% get(c) & PsigHdf$log2FoldChange <(-2),])
		PatMoreCat <- nrow(PsigHdf[row.names(PsigHdf) %in% get(c) & PsigHdf$log2FoldChange >2,])
		MatMoreNT <- nrow(PsigHdf[row.names(PsigHdf) %in% NT & PsigHdf$log2FoldChange <(-2),])
		PatMoreNT <- nrow(PsigHdf[row.names(PsigHdf) %in% NT & PsigHdf$log2FoldChange >2,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoCompH[nrow(homoeoCompH)+1,] <- rbind(c,s,"twofold",mat[1,1],mat[2,1],ft$p.value) }}
homoeoCompH[nrow(homoeoCompH)+1,] <- rbind("NT","NA","twofold",mat[1,2],mat[2,2],"NA") 


HCpubtableH <- homoeoCompH %>% 
          mutate(sided=if_else(pvalueDE <= 0.05, paste0(sided,"_SIGNIF"), paste0(sided,"_notSig"))) %>% 
          select(-pvalueDE)

write.table(HCpubtableH,file="publicationFisherHomoeolog.Ahyp.tbl", row.names=F, col.names=T, sep="\t", quote=F)




PsigAiddf <- as.data.frame(PsigAid)

homoeoCompAid <- data.frame(category=character(),sided=character(), FC=character(0),
	MatMoreCat=integer(), PatMoreCat=integer(), 
	pvalueDE=double(), stringsAsFactors=F)
for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(PsigAiddf[row.names(PsigAiddf) %in% get(c) & PsigAiddf$log2FoldChange <0,])
		PatMoreCat <- nrow(PsigAiddf[row.names(PsigAiddf) %in% get(c) & PsigAiddf$log2FoldChange >0,])
		MatMoreNT <- nrow(PsigAiddf[row.names(PsigAiddf) %in% NT & PsigAiddf$log2FoldChange <0,])
		PatMoreNT <- nrow(PsigAiddf[row.names(PsigAiddf) %in% NT & PsigAiddf$log2FoldChange >0,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoCompAid[nrow(homoeoCompAid)+1,] <- rbind(c,s,"any",mat[1,1],mat[2,1],ft$p.value) } }
homoeoCompAid[nrow(homoeoCompAid)+1,] <- rbind("NT","NA","any",mat[1,2],mat[2,2],"NA") 


for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(PsigAiddf[row.names(PsigAiddf) %in% get(c) & PsigAiddf$log2FoldChange <(-2),])
		PatMoreCat <- nrow(PsigAiddf[row.names(PsigAiddf) %in% get(c) & PsigAiddf$log2FoldChange >2,])
		MatMoreNT <- nrow(PsigAiddf[row.names(PsigAiddf) %in% NT & PsigAiddf$log2FoldChange <(-2),])
		PatMoreNT <- nrow(PsigAiddf[row.names(PsigAiddf) %in% NT & PsigAiddf$log2FoldChange >2,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoCompAid[nrow(homoeoCompAid)+1,] <- rbind(c,s,"twofold",mat[1,1],mat[2,1],ft$p.value) }}
homoeoCompAid[nrow(homoeoCompAid)+1,] <- rbind("NT","NA","twofold",mat[1,2],mat[2,2],"NA") 


HCpubtableAid <- homoeoCompAid %>% 
          mutate(sided=if_else(pvalueDE <= 0.05, paste0(sided,"_SIGNIF"), paste0(sided,"_notSig"))) %>% 
          select(-pvalueDE)

write.table(HCpubtableAid,file="publicationFisherHomoeolog.Aid.tbl", row.names=F, col.names=T, sep="\t", quote=F)







ls()
rm(list = ls()[grep("_", ls())])

