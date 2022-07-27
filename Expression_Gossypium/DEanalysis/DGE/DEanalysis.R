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

setwd("W:/corrinne/Cytonuclear/cotton/DEanalysis")


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
countdata <- countmerge("W:/corrinne/Cytonuclear/cotton/DEanalysis/counts")
tpmdata <- tpmmerge("W:/corrinne/Cytonuclear/cotton/DEanalysis/counts")

write.table(countdata, file="count.cotton.tbl", quote=F, sep="\t", row.names=F)
write.table(tpmdata, file="tpm.cotton.tbl", quote=F, row.names=F, sep="\t")

### generate metadata framework
metadata <- data.frame(id=names(countdata[,-1]), data.frame(id=names(countdata[,-1])) %>% separate(id, c("species", "rep", NA)))
countdata[,-1] <- round(countdata[,-1],0)

### define cp/mt genes
cpmt <- countdata$target_id[!grepl("Gohir.",countdata$target_id)]
cp <- cpmt[grepl("cp_",cpmt)]
mt <- cpmt[!grepl("cp_",cpmt)]

### generate df and metadata for subgenomes
category <- read.table("gene.pairs", col.names=c("category","dad","mom")) %>% 
     filter(dad %in% countdata$target_id) %>% filter(mom %in% countdata$target_id)
category$both <- paste0(category$dad,"_",category$mom)


dadCategory <- category %>%
	select(category,dad) %>%
	mutate(subgenome="dad") %>%
	dplyr::rename(gene=dad) 

momCategory <- category %>%
	select(category,mom) %>%
	mutate(subgenome="mom") %>%
	dplyr::rename(gene=mom) 

fullCategory <- rbind(dadCategory,momCategory) %>% 
	filter(gene %in% countdata$target_id) %>%
	add_row(gene=cpmt, subgenome="NotApp", category="NotApp")


#fullcategory <- read.table("cotton.expanded.category.longForm", col.names=c("category","CyMira","CyMiraSub","gene","subgenome")) %>% 
#     filter(gene %in% countdata$target_id) %>% select(gene,subgenome,category,CyMira,CyMiraSub) %>%
#     add_row(gene=cpmt, subgenome="NotApp", category="NotApp",CyMira="NotApp", CyMiraSub="NotApp")


Acountdata <- countdata %>% filter(target_id %in% category$mom) %>% select(!contains("D5"))
Dcountdata <- countdata %>% filter(target_id %in% category$dad) %>% select(!contains("A2"))

Ametadata <- filter(metadata, !grepl("D5",species))
Dmetadata <- filter(metadata, !grepl("A2",species)) 

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
momSubgenomes <- c("AD1A","AD2A")
dadSubgenomes <- c("AD1D","AD2D")
sumGenomes <- c("AD1SA","AD2SA","AD1SD","AD2SD")
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
Rcountdata <- Rcountdata[-(1:2),c(1,7:16)] 

for (i in 1:nrow(homoeoPairs)) {
    Rcountdata[nrow(Rcountdata)+1,] <- rbind(
        cbind(paste0(homoeoPairs[i,1],"_",homoeoPairs[i,2]),
        countdata[countdata$target_id==homoeoPairs[i,2], 7:16]/ 
        countdata[countdata$target_id==homoeoPairs[i,1], 7:16]))
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
lmTable <- as.data.frame(assay(rld)) %>% rownames_to_column(var = "gene") %>% full_join(fullCategory) %>% drop_na()
write.table(lmTable,file="cotton.lm.tbl", quote=F, sep="\t")

##### ##### ##### ##### ##### ##### ##### #####
##### testing what is wrong with tpm ##### ####
##### ##### ##### ##### ##### ##### ##### ##### 

# homoeolog.test <- c(category$mom, category$dad)
# temp.test <- countdata %>% filter(target_id %in% homoeolog.test)
# temp.test$genome <- gsub("Gohir[.]","",temp.test$target_id)
# temp.test$genome <- gsub("[01][0-9]G.*","",temp.test$genome)

# Z <- colSums(temp.test[temp.test$genome=="A",2:21])/(colSums(temp.test[,2:21]))
# Y <- colSums(temp.test[temp.test$genome=="D",2:21])/(colSums(temp.test[,2:21]))

# temp.test.all <- countdata %>% filter(!target_id %in% cpmt)
# temp.test.all$genome <- gsub("Gohir[.]","",temp.test.all$target_id)
# temp.test.all$genome <- gsub("[01][0-9]G.*","",temp.test.all$genome)

# X <- colSums(temp.test.all[temp.test.all$genome=="A",2:21])/(colSums(temp.test.all[,2:21]))
# W <- colSums(temp.test.all[temp.test.all$genome=="D",2:21])/(colSums(temp.test.all[,2:21]))


# temp.test.tpm <- tpmdata %>% filter(target_id %in% homoeolog.test)
# temp.test.tpm$genome <- gsub("Gohir[.]","",temp.test.tpm$target_id)
# temp.test.tpm$genome <- gsub("[01][0-9]G.*","",temp.test.tpm$genome)

# A <- colSums(temp.test.tpm[temp.test.tpm$genome=="A",2:21])/(colSums(temp.test.tpm[,2:21]))
# B <- colSums(temp.test.tpm[temp.test.tpm$genome=="D",2:21])/(colSums(temp.test.tpm[,2:21]))

# temp.test.tpm.all <- tpmdata %>% filter(!target_id %in% cpmt)
# temp.test.tpm.all$genome <- gsub("Gohir[.]","",temp.test.tpm.all$target_id)
# temp.test.tpm.all$genome <- gsub("[01][0-9]G.*","",temp.test.tpm.all$genome)

# C <- colSums(temp.test.tpm.all[temp.test.tpm.all$genome=="A",2:21])/(colSums(temp.test.tpm.all[,2:21]))
# D <- colSums(temp.test.tpm.all[temp.test.tpm.all$genome=="D",2:21])/(colSums(temp.test.tpm.all[,2:21]))

# rbind(Z,Y,X,W,A,B,C,D)


##### ##### ##### ##### ##### ##### 
##### Maternal (A) comparison ##### 
##### ##### ##### ##### ##### ##### 

Adds <- ddsMake(Acountdata,Ametadata)
PCApheat(Adds,"A","no")

resultsNames(Adds)
Ares <- results(Adds)
summary(Ares) 

### grab DE results, At relative to D5, no cp/mt
AD1resA <- results(Adds, contrast=c("species","AD1","A2"), alpha=0.05)
AD2resA <- results(Adds, contrast=c("species","AD2","A2"), alpha=0.05)

AD1resA <- AD1resA[order(AD1resA$padj),]
AD2resA <- AD2resA[order(AD2resA$padj),]

AD1sigA <- subset(AD1resA, padj<0.05)
AD2sigA <- subset(AD2resA, padj<0.05)

write.table(AD1sigA,file="DEGsigResults.AD1AvA.tbl", quote=F, sep="\t")
write.table(AD2sigA,file="DEGsigResults.AD2AvA.tbl", quote=F, sep="\t")

### just a couple more checks, for fun
resLFCAD1 <- lfcShrink(Adds,coef="species_AD1_vs_A2", type="apeglm")
plotMA(resLFCAD1, ylim = c(-5, 5))

resLFCAD2 <- lfcShrink(Adds,coef="species_AD2_vs_A2", type="apeglm")
plotMA(resLFCAD2, ylim = c(-5, 5))

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
AD1resD <- results(Ddds, contrast=c("species","AD1","D5"), alpha=0.05)
AD2resD <- results(Ddds, contrast=c("species","AD2","D5"), alpha=0.05)

AD1resD <- AD1resD[order(AD1resD$padj),]
AD2resD <- AD2resD[order(AD2resD$padj),]

AD1sigD <- subset(AD1resD, padj<0.05)
AD2sigD <- subset(AD2resD, padj<0.05)

write.table(AD1sigD,file="DEGsigResults.AD1DvD.tbl", quote=F, sep="\t")
write.table(AD2sigD,file="DEGsigResults.AD2DvD.tbl", quote=F, sep="\t")

### just a couple more checks, for fun
resLFCAD1D <- lfcShrink(Ddds,coef="species_D5_vs_AD1", type="apeglm")
plotMA(resLFCAD1D, ylim = c(-5, 5))

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

### grab DE results, Dt relative to D5, no cp/mt

AD1resSA <- results(Sdds, contrast=c("species","AD1","A2"), alpha=0.05)
AD2resSA <- results(Sdds, contrast=c("species","AD2","A2"), alpha=0.05)
AD1resSD <- results(Sdds, contrast=c("species","AD1","D5"), alpha=0.05)
AD2resSD <- results(Sdds, contrast=c("species","AD2","D5"), alpha=0.05)

AD1resSA <- AD1resSA[order(AD1resSA$padj),]
AD2resSA <- AD2resSA[order(AD2resSA$padj),]
AD1resSD <- AD1resSD[order(AD1resSD$padj),]
AD2resSD <- AD2resSD[order(AD2resSD$padj),]

AD1sigSA <- subset(AD1resSA, padj<0.05)
AD2sigSA <- subset(AD2resSA, padj<0.05)
AD1sigSD <- subset(AD1resSD, padj<0.05)
AD2sigSD <- subset(AD2resSD, padj<0.05)

write.table(AD1sigSA,file="DEGsigResults.AD1SAvA.tbl", quote=F, sep="\t")
write.table(AD2sigSA,file="DEGsigResults.AD2SAvA.tbl", quote=F, sep="\t")
write.table(AD1sigSD,file="DEGsigResults.AD1SDvD.tbl", quote=F, sep="\t")
write.table(AD2sigSD,file="DEGsigResults.AD2SDvD.tbl", quote=F, sep="\t")



############################################
##### get expression relationships #########
#####      vis-a-vis dominance     #########
############################################

A2D5res <- results(Sdds, contrast=c("species","A2","D5"), alpha=0.05)

cpmtcat <- rbind(
	category %>% select(category,both) %>% dplyr::rename(gene=both),
	as.data.frame(cbind(gene=cpmt,category=gsub("_.*","",cpmt)))
)
	
is.gt.AD1 <- data.frame(gene=character(),PvM=character(),PvD=character(),MvD=character(), stringsAsFactors=F)

for (gene in Scountdata$target_id) {
	is.gt.AD1[nrow(is.gt.AD1)+1,"gene"] <- gene
	is.gt.AD1[is.gt.AD1$gene==gene,"MvD"] <- as.data.frame(A2D5res)[gene,] %>% 
   	 mutate(MvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'mom',
               padj < 0.05 & log2FoldChange < 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(MvD)

	is.gt.AD1[is.gt.AD1$gene==gene,"PvM"] <- as.data.frame(AD1resSA)[gene,] %>% 
   	 mutate(PvM = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'mom',
               TRUE ~ 'error' ) ) %>% select(PvM)

	is.gt.AD1[is.gt.AD1$gene==gene,"PvD"] <- as.data.frame(AD1resSD)[gene,] %>% 
   	 mutate(PvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(PvD)
}

is.gt.AD1 <- is.gt.AD1 %>% full_join(.,cpmtcat)


dominanceCatsAD1 <- plyr::count(is.gt.AD1, vars = c('PvM','PvD','MvD')) %>% 
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



dominanceCatsCytoAD1 <- plyr::count(is.gt.AD1, vars = c('PvM','PvD','MvD','category')) %>% 
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


dominanceTableAD1 <- bind_rows(dominanceCatsAD1,dominanceCatsCytoAD1) %>% 
	select(-c(PvM,PvD,MvD)) %>% 
	filter(type != "undet") %>%
	filter(type != "no_info") %>%
	pivot_wider(names_from=type, values_from=freq) %>%
	tidyr::replace_na(list(category="all")) %>%
	replace(is.na(.), 0) %>%
	select(category,category01, category12, category02, category11, category04, category09, category03, category07, category10, category05, category06, category08, all_equal) 



write.table(dominanceTableAD1,file="AD1ExpressionDominance.tbl", quote=F, row.names=F, sep="\t")


	
is.gt.AD2 <- data.frame(gene=character(),PvM=character(),PvD=character(),MvD=character(), stringsAsFactors=F)

for (gene in Scountdata$target_id) {
	is.gt.AD2[nrow(is.gt.AD2)+1,"gene"] <- gene
	is.gt.AD2[is.gt.AD2$gene==gene,"MvD"] <- as.data.frame(A2D5res)[gene,] %>% 
   	 mutate(MvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'mom',
               padj < 0.05 & log2FoldChange < 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(MvD)

	is.gt.AD2[is.gt.AD2$gene==gene,"PvM"] <- as.data.frame(AD2resSA)[gene,] %>% 
   	 mutate(PvM = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'mom',
               TRUE ~ 'error' ) ) %>% select(PvM)

	is.gt.AD2[is.gt.AD2$gene==gene,"PvD"] <- as.data.frame(AD2resSD)[gene,] %>% 
   	 mutate(PvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(PvD)
}

is.gt.AD2 <- is.gt.AD2 %>% full_join(.,cpmtcat)


dominanceCatsAD2 <- plyr::count(is.gt.AD2, vars = c('PvM','PvD','MvD')) %>% 
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



dominanceCatsCytoAD2 <- plyr::count(is.gt.AD2, vars = c('PvM','PvD','MvD','category')) %>% 
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


dominanceTableAD2 <- bind_rows(dominanceCatsAD2,dominanceCatsCytoAD2) %>% 
	select(-c(PvM,PvD,MvD)) %>% 
	filter(type != "undet") %>%
	filter(type != "no_info") %>%
	pivot_wider(names_from=type, values_from=freq) %>%
	tidyr::replace_na(list(category="all")) %>%
	replace(is.na(.), 0) %>%
	select(category,category01, category12, category02, category11, category04, category09, category03, category07, category10, category05, category06, category08, all_equal) 


write.table(dominanceTableAD2,file="AD2ExpressionDominance.tbl", quote=F, row.names=F, sep="\t")



### get the significant DE genes each category
## also parsed by up/down regulation

AD1A <- read.table("DEGsigResults.AD1AvA.tbl")
AD2A <- read.table("DEGsigResults.AD2AvA.tbl")
AD1D <- read.table("DEGsigResults.AD1DvD.tbl")
AD2D <- read.table("DEGsigResults.AD2DvD.tbl")

AD1SA <- read.table("DEGsigResults.AD1SAvA.tbl")
AD2SA <- read.table("DEGsigResults.AD2SAvA.tbl")
AD1SD <- read.table("DEGsigResults.AD1SDvD.tbl")
AD2SD <- read.table("DEGsigResults.AD2SDvD.tbl")


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
write.table(mom,file="ADsAvA.fisherexact.tbl", quote=F, row.names=F, sep="\t")

 
           
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
write.table(dad,file="ADsDvD.fisherexact.tbl", quote=F, row.names=F, sep="\t")




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
write.table(summed,file="ADSvAorD.fisherexact.tbl", quote=F, row.names=F, sep="\t")


pubtable <- rbind(mom,dad,summed) %>% group_by(species, category, direction) %>% 
          slice_min(pvalueDE,n=1) %>% 
          mutate(sided=if_else(pvalueDE <= 0.05, paste0(sided,"_SIGNIF"), paste0(sided,"_notSig"))) %>% 
          select(-pvalueDE)

write.table(pubtable,file="publicationFisherCotton.tbl", row.names=F, col.names=T, sep="\t", quote=F)
`
############################
##################################
######################################
add Kruskal-Wallis tests??? maybe not because categories/genes are not independent?
######################################
##################################
############################

categoryP <- dplyr::rename(category,target_id=both) %>% select(target_id,category) 

ratioSum <- Rcountdata %>% rowwise() %>% 
     transmute(target_id,  
     AD1mean = mean(c(AD1_1,AD1_2,AD1_3,AD1_4,AD1_5), na.rm=TRUE),
     AD1median = median(c(AD1_1,AD1_2,AD1_3,AD1_4,AD1_5), na.rm=TRUE),
     AD2mean = mean(c(AD2_1,AD2_2,AD2_3,AD2_4,AD2_5), na.rm=TRUE),
     AD2median = median(c(AD2_1,AD2_2,AD2_3,AD2_4,AD2_5), na.rm=TRUE),
     AD1NAcnt= sum(is.na(c(AD1_1,AD1_2,AD1_3,AD1_4,AD1_5))),
     AD2NAcnt= sum(is.na(c(AD2_1,AD2_2,AD2_3,AD2_4,AD2_5)))) %>%
     left_join(., categoryP, by='target_id') %>% 
     gather("tmp", "value", starts_with('AD')) %>% 
     separate(tmp, into=c('species', 'type' ), sep=c(3))

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

PcountdataAD1 <- category %>% 
	select(dad,mom,both) %>% 
	dplyr::rename(target_id=dad) %>%
     	left_join(., Dcountdata, by='target_id') %>% 
	select(c(mom,both,starts_with('AD1'))) %>%
	rename_at(vars(starts_with('AD1')), ~ str_replace(., 'AD1', 'AD1Dad')) %>% 
	dplyr::rename(target_id=mom) %>%
     	left_join(., Acountdata, by='target_id') %>%
	select(c(both,starts_with('AD1'))) %>%
	rename_at(vars(starts_with('AD1_')), ~ str_replace(., 'AD1', 'AD1Mom')) %>% 
	column_to_rownames(var = "both") 	

PmetadataAD1 <- data.frame(id=names(PcountdataAD1), data.frame(id=names(PcountdataAD1)) %>% separate(id, c("species", "rep", NA)))

PddsAD1 <- DESeqDataSetFromMatrix(countData = round(PcountdataAD1,0), colData=PmetadataAD1, design = ~species) 

libTotalAD1 <- rep(colSums(countdata[,7:11]),2) %>%
	set_names(names(PcountdataAD1))
libSizeAD1 <- libTotalAD1/mean(libTotalAD1)

sizeFactors(PddsAD1) <- libSizeAD1
PddsAD1 <- DESeq(PddsAD1)

PresAD1 <- results(DESeq(PddsAD1),contrast=c("species","AD1Dad","AD1Mom"),alpha=0.05)

PsigAD1 <- subset(PresAD1[order(PresAD1$log2FoldChange),], padj<0.05)

write.table(PsigAD1,file="DEGsigResults.AD1PatvMat.tbl", quote=F, sep="\t")



PcountdataAD2 <- category %>% 
	select(dad,mom,both) %>% 
	dplyr::rename(target_id=dad) %>%
     	left_join(., Dcountdata, by='target_id') %>% 
	select(c(mom,both,starts_with('AD2'))) %>%
	rename_at(vars(starts_with('AD2')), ~ str_replace(., 'AD2', 'AD2Dad')) %>% 
	dplyr::rename(target_id=mom) %>%
     	left_join(., Acountdata, by='target_id') %>%
	select(c(both,starts_with('AD2'))) %>%
	rename_at(vars(starts_with('AD2_')), ~ str_replace(., 'AD2', 'AD2Mom')) %>% 
	column_to_rownames(var = "both") 	

PmetadataAD2 <- data.frame(id=names(PcountdataAD2), data.frame(id=names(PcountdataAD2)) %>% separate(id, c("species", "rep", NA)))

PddsAD2 <- DESeqDataSetFromMatrix(countData = round(PcountdataAD2,0), colData=PmetadataAD2, design = ~species) 

libTotalAD2 <- rep(colSums(countdata[,12:16]),2) %>%
	set_names(names(PcountdataAD2))
libSizeAD2 <- libTotalAD2/mean(libTotalAD2)

sizeFactors(PddsAD2) <- libSizeAD2
PddsAD2 <- DESeq(PddsAD2)

PresAD2 <- results(DESeq(PddsAD2),contrast=c("species","AD2Dad","AD2Mom"),alpha=0.05)

PsigAD2 <- subset(PresAD2[order(PresAD2$log2FoldChange),], padj<0.05)

write.table(PsigAD2,file="DEGsigResults.AD2PatvMat.tbl", quote=F, sep="\t")



Bcountdata <- category %>% 
	select(dad,mom,both) %>% 
	dplyr::rename(target_id=dad) %>%
     	left_join(., Dcountdata, by='target_id') %>% 
	select(c(mom,both,starts_with('D5'))) %>%
	dplyr::rename(target_id=mom) %>%
     	left_join(., Acountdata, by='target_id') %>%
	select(c(both,starts_with(c('D5','A2')))) %>% 
	dplyr::rename(target_id=both)	

Bmetadata <- data.frame(id=names(Bcountdata[,-1]), data.frame(id=names(Bcountdata[,-1])) %>% separate(id, c("species", "rep", NA)))

Bdds <- ddsMake(Bcountdata,Bmetadata)
Bres <- results(Bdds,contrast=c("species","D5","A2"),alpha=0.05)
Bsig <- subset(Bres[order(Bres$log2FoldChange),], padj<0.05)
write.table(Bsig,file="DEGsigResults.D5vA2.tbl", quote=F, sep="\t")


# AssigA
# 1429

# AssigD
# 1977

# Psig
# 1777

# Bsig
# 3103

options(scipen=10)

ResCompare <- as.data.frame(PresAD1) %>% 
	rownames_to_column(var='both') %>%
	dplyr::arrange(desc(abs(log2FoldChange))) %>%
	mutate(AD1Rank=1:nrow(PresAD1)) %>%
     	right_join(., category, by='both') %>% 
	dplyr::rename(target_id=both, AD1PatMatFC=log2FoldChange, AD1PatMatpadj=padj) %>%
     	left_join(., (as.data.frame(PresAD2) %>% rownames_to_column(var='target_id') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(AD2Rank=1:nrow(PresAD2))), by='target_id') %>% 
	dplyr::rename(AD2PatMatFC=log2FoldChange, AD2PatMatpadj=padj) %>%
     	left_join(., (as.data.frame(AD1resA) %>% rownames_to_column(var='mom') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(AD1MomRank=1:nrow(AD1resA))), by='mom') %>% 
	dplyr::rename(AD1MatMomFC=log2FoldChange, AD1MatMompadj=padj) %>%
     	left_join(., (as.data.frame(AD2resA) %>% rownames_to_column(var='mom') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(AD2MomRank=1:nrow(AD2resA))), by='mom') %>% 
	dplyr::rename(AD2MatMomFC=log2FoldChange, AD2MatMompadj=padj) %>%
     	left_join(., (as.data.frame(AD1resD) %>% rownames_to_column(var='dad') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(AD1DadRank=1:nrow(AD1resD))), by='dad') %>% 
	dplyr::rename(AD1PatDadFC=log2FoldChange, AD1PatDadpadj=padj) %>% 
     	left_join(., (as.data.frame(AD2resD) %>% rownames_to_column(var='dad') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(AD2DadRank=1:nrow(AD2resD))), by='dad') %>% 
	dplyr::rename(AD2PatDadFC=log2FoldChange, AD2PatDadpadj=padj) %>%
     	left_join(., (as.data.frame(Bres) %>% rownames_to_column(var='target_id') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(DipRank=1:nrow(Bres))), by='target_id') %>% 
	dplyr::rename(DadMomFC=log2FoldChange,DadMompadj=padj) %>% 
	select(matches('target_id|category|AD1|AD2|Rank|DadMom'))



#	filter(PolyRank!='NA' | MomRank!='NA' | DadRank!='NA' | DipRank!='NA') %>%
#	filter(PolyRank<177 | MomRank<142 | DadRank<197 | DipRank<310)




write.table(ResCompare,file="DEGResults.ResCompare.cotton.tbl", quote=F, sep="\t", row.names=F)









## make contingency tables for homoeolog comparison
PsigAD1df <- as.data.frame(PsigAD1)

homoeoCompAD1 <- data.frame(category=character(),sided=character(), FC=character(0),
	MatMoreCat=integer(), PatMoreCat=integer(), 
	pvalueDE=double(), stringsAsFactors=F)
for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(PsigAD1df[row.names(PsigAD1df) %in% get(c) & PsigAD1df$log2FoldChange <0,])
		PatMoreCat <- nrow(PsigAD1df[row.names(PsigAD1df) %in% get(c) & PsigAD1df$log2FoldChange >0,])
		MatMoreNT <- nrow(PsigAD1df[row.names(PsigAD1df) %in% NT & PsigAD1df$log2FoldChange <0,])
		PatMoreNT <- nrow(PsigAD1df[row.names(PsigAD1df) %in% NT & PsigAD1df$log2FoldChange >0,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoCompAD1[nrow(homoeoCompAD1)+1,] <- rbind(c,s,"any",mat[1,1],mat[2,1],ft$p.value) } }
homoeoCompAD1[nrow(homoeoCompAD1)+1,] <- rbind("NT","NA","any",mat[1,2],mat[2,2],"NA") 


for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(PsigAD1df[row.names(PsigAD1df) %in% get(c) & PsigAD1df$log2FoldChange <(-2),])
		PatMoreCat <- nrow(PsigAD1df[row.names(PsigAD1df) %in% get(c) & PsigAD1df$log2FoldChange >2,])
		MatMoreNT <- nrow(PsigAD1df[row.names(PsigAD1df) %in% NT & PsigAD1df$log2FoldChange <(-2),])
		PatMoreNT <- nrow(PsigAD1df[row.names(PsigAD1df) %in% NT & PsigAD1df$log2FoldChange >2,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoCompAD1[nrow(homoeoCompAD1)+1,] <- rbind(c,s,"twofold",mat[1,1],mat[2,1],ft$p.value) }}
homoeoCompAD1[nrow(homoeoCompAD1)+1,] <- rbind("NT","NA","twofold",mat[1,2],mat[2,2],"NA") 


HCpubtableAD1 <- homoeoCompAD1 %>% 
          mutate(sided=if_else(pvalueDE <= 0.05, paste0(sided,"_SIGNIF"), paste0(sided,"_notSig"))) %>% 
          select(-pvalueDE)

write.table(HCpubtableAD1,file="publicationFisherHomoeolog.AD1.tbl", row.names=F, col.names=T, sep="\t", quote=F)




PsigAD2df <- as.data.frame(PsigAD2)

homoeoCompAD2 <- data.frame(category=character(),sided=character(), FC=character(0),
	MatMoreCat=integer(), PatMoreCat=integer(), 
	pvalueDE=double(), stringsAsFactors=F)
for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(PsigAD2df[row.names(PsigAD2df) %in% get(c) & PsigAD2df$log2FoldChange <0,])
		PatMoreCat <- nrow(PsigAD2df[row.names(PsigAD2df) %in% get(c) & PsigAD2df$log2FoldChange >0,])
		MatMoreNT <- nrow(PsigAD2df[row.names(PsigAD2df) %in% NT & PsigAD2df$log2FoldChange <0,])
		PatMoreNT <- nrow(PsigAD2df[row.names(PsigAD2df) %in% NT & PsigAD2df$log2FoldChange >0,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoCompAD2[nrow(homoeoCompAD2)+1,] <- rbind(c,s,"any",mat[1,1],mat[2,1],ft$p.value) } }
homoeoCompAD2[nrow(homoeoCompAD2)+1,] <- rbind("NT","NA","any",mat[1,2],mat[2,2],"NA") 


for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(PsigAD2df[row.names(PsigAD2df) %in% get(c) & PsigAD2df$log2FoldChange <(-2),])
		PatMoreCat <- nrow(PsigAD2df[row.names(PsigAD2df) %in% get(c) & PsigAD2df$log2FoldChange >2,])
		MatMoreNT <- nrow(PsigAD2df[row.names(PsigAD2df) %in% NT & PsigAD2df$log2FoldChange <(-2),])
		PatMoreNT <- nrow(PsigAD2df[row.names(PsigAD2df) %in% NT & PsigAD2df$log2FoldChange >2,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoCompAD2[nrow(homoeoCompAD2)+1,] <- rbind(c,s,"twofold",mat[1,1],mat[2,1],ft$p.value) }}
homoeoCompAD2[nrow(homoeoCompAD2)+1,] <- rbind("NT","NA","twofold",mat[1,2],mat[2,2],"NA") 


HCpubtableAD2 <- homoeoCompAD2 %>% 
          mutate(sided=if_else(pvalueDE <= 0.05, paste0(sided,"_SIGNIF"), paste0(sided,"_notSig"))) %>% 
          select(-pvalueDE)

write.table(HCpubtableAD2,file="publicationFisherHomoeolog.AD2.tbl", row.names=F, col.names=T, sep="\t", quote=F)



ls()
rm(list = ls()[grep("_", ls())])

