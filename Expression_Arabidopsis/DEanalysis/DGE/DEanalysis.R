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

setwd("W:/corrinne/Cytonuclear/arabidopsis/DEanalysis12kV2")


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
countdata <- countmerge("W:/corrinne/Cytonuclear/arabidopsis/DEanalysis/counts")
tpmdata <- tpmmerge("W:/corrinne/Cytonuclear/arabidopsis/DEanalysis/counts")

write.table(countdata, file="count.arabidopsis.tbl", quote=F, row.names=F, sep="\t")
write.table(tpmdata, file="tpm.arabidopsis.tbl", quote=F, row.names=F, sep="\t")

### generate metadata framework
metadata <- data.frame(id=names(countdata[,-1]), data.frame(id=names(countdata[,-1])) %>% separate(id, c("species", "rep", NA)))
countdata[,-1] <- round(countdata[,-1],0)

### define cp/mt genes
cpmt <- countdata$target_id[!grepl("AsA.",countdata$target_id)]
cp <- cpmt[grepl("cp_",cpmt)]
mt <- cpmt[!grepl("cp_",cpmt)]

### generate df and metadata for subgenomes
### AsAa is dad
### AsAt is mom
category <- read.table("gene.pairs", col.names=c("category","dad","mom")) %>% 
     filter(dad %in% countdata$target_id) %>% filter(mom %in% countdata$target_id)
category$both <- paste0(category$dad,"_",category$mom)

fullcategory <- read.table("arabidopsis.expanded.category.longForm", col.names=c("category","CyMira","CyMiraSub","gene","subgenome")) %>% 
     filter(gene %in% countdata$target_id) %>% select(gene,subgenome,category,CyMira,CyMiraSub) %>%
     add_row(gene=cpmt, subgenome="NotApp", category="NotApp",CyMira="NotApp", CyMiraSub="NotApp")

Acountdata <- countdata %>% filter(target_id %in% category$mom) %>% select(!contains("Aare"))
Dcountdata <- countdata %>% filter(target_id %in% category$dad) %>% select(!contains("Atha"))

Ametadata <- filter(metadata, !grepl("Aare",species))
Dmetadata <- filter(metadata, !grepl("Atha",species)) 

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
momSubgenomes <- c("AsAt")
dadSubgenomes <- c("AsAa")
sumGenomes <- c("AsST","AsSA")
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
Rcountdata <- Rcountdata[-(1:2),c(1,9:13)] 

for (i in 1:nrow(homoeoPairs)) {
    Rcountdata[nrow(Rcountdata)+1,] <- rbind(
        cbind(paste0(homoeoPairs[i,1],"_",homoeoPairs[i,2]),
        countdata[countdata$target_id==homoeoPairs[i,2], 9:13]/ 
        countdata[countdata$target_id==homoeoPairs[i,1], 9:13]))
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
write.table(lmTable,file="arabidopsis.lm.tbl", quote=F, sep="\t", row.names=F)


##### ##### ##### ##### ##### ##### 
##### Maternal (A) comparison ##### 
##### ##### ##### ##### ##### ##### 

Adds <- ddsMake(Acountdata,Ametadata)
PCApheat(Adds,"A","no")

resultsNames(Adds)
Ares <- results(Adds)
summary(Ares) 

### grab DE results, Atha relative to Asue, no cp/mt
AsresA <- results(Adds, contrast=c("species","Asue","Atha"), alpha=0.05)

AsresA <- AsresA[order(AsresA$padj),]

AssigA <- subset(AsresA, padj<0.05)

write.table(AssigA,file="DEGsigResults.AsuevAtha.tsv", quote=F, sep="\t")
write.table(AsresA,file="DEGAllResults.AsuevAtha.tsv", quote=F, sep="\t")



##### ##### ##### ##### ##### ##### 
##### Paternal (D) comparison ##### 
##### ##### ##### ##### ##### ##### 

Ddds <- ddsMake(Dcountdata,Dmetadata)
PCApheat(Ddds,"D","no")

resultsNames(Ddds)
Dres <- results(Ddds)
summary(Dres) 

### grab DE results, Dt relative to D5, no cp/mt
AsresD <- results(Ddds, contrast=c("species","Asue","Aare"), alpha=0.05)

AsresD <- AsresD[order(AsresD$padj),]

AssigD <- subset(AsresD, padj<0.05)

write.table(AssigD,file="DEGsigResults.AsuevAare.tsv", quote=F, sep="\t")
write.table(AsresD,file="DEGAllResults.AsuevAare.tsv", quote=F, sep="\t")



##### ##### ##### ##### ##### ##### #####
##### diploid-polyploid comparison  ##### 
##### ##### ##### ##### ##### ##### ##### 

Sdds <- ddsMake(Scountdata,metadata)
PCApheat(Sdds,"summed","yes")

resultsNames(Sdds)
Sres <- results(Sdds)
summary(Sres)

### grab DE results, Asue relative to Atha/are

AsresSA <- results(Sdds, contrast=c("species","Asue","Atha"), alpha=0.05)
AsresSD <- results(Sdds, contrast=c("species","Asue","Aare"), alpha=0.05)

AsresSA <- AsresSA[order(AsresSA$padj),]
AsresSD <- AsresSD[order(AsresSD$padj),]

AssigSA <- subset(AsresSA, padj<0.05)
AssigSD <- subset(AsresSD, padj<0.05)

write.table(AssigSA,file="DEGsigResults.AsueSvAtha.tsv", quote=F, sep="\t")
write.table(AssigSD,file="DEGsigResults.AsueSvAare.tsv", quote=F, sep="\t")

write.table(AsresSA,file="DEGAllResults.AsueSvAtha.tsv", quote=F, sep="\t")
write.table(AsresSD,file="DEGAllResults.AsueSvAare.tsv", quote=F, sep="\t")


### get the significant DE genes each category
## also parsed by up/down regulation

AsAt <- read.table("DEGsigResults.AsuevAtha.tsv")
AsAa <- read.table("DEGsigResults.AsuevAare.tsv")

AsST <- read.table("DEGsigResults.AsueSvAtha.tsv")
AsSA <- read.table("DEGsigResults.AsueSvAare.tsv")

############################################
##### get expression relationships #########
#####      vis-a-vis dominance     #########
############################################

AsresTA <- results(Sdds, contrast=c("species","Atha","Aare"), alpha=0.05)

cpmtcat <- rbind(
	category %>% select(category,both) %>% dplyr::rename(gene=both),
	as.data.frame(cbind(gene=cpmt,category=gsub("_.*","",cpmt)))
)
	
is.gt <- data.frame(gene=character(),PvM=character(),PvD=character(),MvD=character(), stringsAsFactors=F)

for (gene in Scountdata$target_id) {
	is.gt[nrow(is.gt)+1,"gene"] <- gene
	is.gt[is.gt$gene==gene,"MvD"] <- as.data.frame(AsresTA)[gene,] %>% 
   	 mutate(MvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'mom',
               padj < 0.05 & log2FoldChange < 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(MvD)

	is.gt[is.gt$gene==gene,"PvM"] <- as.data.frame(AsresSA)[gene,] %>% 
   	 mutate(PvM = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'mom',
               TRUE ~ 'error' ) ) %>% select(PvM)

	is.gt[is.gt$gene==gene,"PvD"] <- as.data.frame(AsresSD)[gene,] %>% 
   	 mutate(PvD = case_when(
       	   padj >= 0.05 ~ 'eq',
               padj < 0.05 & log2FoldChange > 0 ~ 'pp',
               padj < 0.05 & log2FoldChange < 0 ~ 'dad',
               TRUE ~ 'error' ) ) %>% select(PvD)
}

is.gt <- is.gt %>% full_join(.,cpmtcat)


dominanceCats <- plyr::count(is.gt, vars = c('PvM','PvD','MvD')) %>% 
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



dominanceCatsCyto <- plyr::count(is.gt, vars = c('PvM','PvD','MvD','category')) %>% 
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


dominanceTable <- bind_rows(dominanceCats,dominanceCatsCyto) %>% 
	select(-c(PvM,PvD,MvD)) %>% 
	filter(type != "undet") %>%
	filter(type != "no_info") %>%
	pivot_wider(names_from=type, values_from=freq) %>%
	tidyr::replace_na(list(category="all")) %>%
	replace(is.na(.), 0) %>%
	select(category,category01, category12, category02, category11, category04, category09, category03, category07, category10, category05, category06, category08, all_equal) 

write.table(dominanceTable,file="ArabidopsisExpressionDominance.tsv", quote=F, row.names=F, sep="\t")




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
write.table(mom,file="AsuevAtha.fisherexact.tsv", quote=F, row.names=F, sep="\t")

 
           
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
write.table(dad,file="AsuevAare.fisherexact.tsv", quote=F, row.names=F, sep="\t")




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
write.table(summed,file="AsuevAareorAtha.fisherexact.tsv", quote=F, row.names=F, sep="\t")


pubtable <- rbind(mom,dad,summed) %>% group_by(species, category, direction) %>% 
          slice_min(pvalueDE,n=1) %>% 
          mutate(sided=if_else(pvalueDE <= 0.05, paste0(sided,"_SIGNIF"), paste0(sided,"_notSig"))) %>% 
          select(-pvalueDE)

write.table(pubtable,file="publicationFisherArabidopsis.tsv", row.names=F, col.names=T, sep="\t", quote=F)


############################

categoryP <- dplyr::rename(category,target_id=both) %>% select(target_id,category) 

ratioSum <- Rcountdata %>% rowwise() %>% 
     transmute(target_id,  
     Asuemean = mean(c(Asue_1,Asue_2,Asue_3,Asue_4,Asue_5), na.rm=TRUE),
     Asuemedian = median(c(Asue_1,Asue_2,Asue_3,Asue_4,Asue_5), na.rm=TRUE),
     AsueNAcnt= sum(is.na(c(Asue_1,Asue_2,Asue_3,Asue_4,Asue_5)))) %>%
     left_join(., categoryP, by='target_id') %>% 
     gather("tmp", "value", starts_with('As')) %>% 
     separate(tmp, into=c('species', 'type' ), sep=c(4))

ratioMean <- ratioSum %>% 
    filter(!is.na(value)) %>% 
    filter(value != 'Inf') %>% 
    filter(value != '-Inf') %>%
    filter(type == 'mean')

ratioMedian <- ratioSum %>% 
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

Pcountdata <- category %>% 
	select(dad,mom,both) %>% 
	dplyr::rename(target_id=dad) %>%
     	left_join(., Dcountdata, by='target_id') %>% 
	select(c(mom,both,starts_with('Asue'))) %>%
	rename_at(vars(starts_with('Asue')), ~ str_replace(., 'Asue', 'AsueDad')) %>% 
	dplyr::rename(target_id=mom) %>%
     	left_join(., Acountdata, by='target_id') %>%
	select(c(both,starts_with('Asue'))) %>%
	rename_at(vars(starts_with('Asue_')), ~ str_replace(., 'Asue', 'AsueMom')) %>% 
	column_to_rownames(var = "both") 	

Pmetadata <- data.frame(id=names(Pcountdata), data.frame(id=names(Pcountdata)) %>% separate(id, c("species", "rep", NA)))

Pdds <- DESeqDataSetFromMatrix(countData = round(Pcountdata,0), colData=Pmetadata, design = ~species) 

libTotal <- rep(colSums(countdata[,9:13]),2) %>%
	set_names(names(Pcountdata))
libSize <- libTotal/mean(libTotal)

sizeFactors(Pdds) <- libSize
Pdds <- DESeq(Pdds)

Pres <- results(DESeq(Pdds),contrast=c("species","AsueDad","AsueMom"),alpha=0.05)

Psig <- subset(Pres[order(Pres$log2FoldChange),], padj<0.05)

write.table(Psig,file="DEGsigResults.AsueDadvMom.tsv", quote=F, sep="\t")





Bcountdata <- category %>% 
	select(dad,mom,both) %>% 
	dplyr::rename(target_id=dad) %>%
     	left_join(., Dcountdata, by='target_id') %>% 
	select(c(mom,both,starts_with('Aare'))) %>%
	dplyr::rename(target_id=mom) %>%
     	left_join(., Acountdata, by='target_id') %>%
	select(c(both,starts_with(c('Aare','Atha')))) %>% 
	dplyr::rename(target_id=both)	

Bmetadata <- data.frame(id=names(Bcountdata[,-1]), data.frame(id=names(Bcountdata[,-1])) %>% separate(id, c("species", "rep", NA)))

Bdds <- ddsMake(Bcountdata,Bmetadata)
Bres <- results(Bdds,contrast=c("species","Aare","Atha"),alpha=0.05)
Bsig <- subset(Bres[order(Bres$log2FoldChange),], padj<0.05)
write.table(Bsig,file="DEGsigResults.AarevAtha.tsv", quote=F, sep="\t")


# AssigA
# 2856

# AssigD
# 3895

# Psig
# 3549

# Bsig
# 6317

SigCompare <- as.data.frame(Psig) %>% 
	rownames_to_column(var='both') %>%
	dplyr::arrange(desc(abs(log2FoldChange))) %>%
	mutate(PolyRank=1:nrow(Psig)) %>%
     	right_join(., category, by='both') %>% 
	select(c(both,log2FoldChange,category,dad,mom,PolyRank)) %>%
	dplyr::rename(target_id=both, PatMatFC=log2FoldChange) %>%
     	left_join(., (as.data.frame(AssigA) %>% rownames_to_column(var='mom') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(MomRank=1:nrow(AssigA))), by='mom') %>% 
	dplyr::rename(MatMomFC=log2FoldChange) %>%
     	left_join(., (as.data.frame(AssigD) %>% rownames_to_column(var='dad') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(DadRank=1:nrow(AssigD))), by='dad') %>% 
	dplyr::rename(PatDadFC=log2FoldChange) %>%
     	left_join(., (as.data.frame(Bsig) %>% rownames_to_column(var='target_id') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(DipRank=1:nrow(Bsig))), by='target_id') %>% 
	dplyr::rename(DadMomFC=log2FoldChange) %>%
	select(c(target_id,category,PatMatFC,DadMomFC,MatMomFC,PatDadFC,PolyRank,DipRank,MomRank,DadRank)) %>%
	filter(PolyRank!='NA' | MomRank!='NA' | DadRank!='NA' | DipRank!='NA') %>%
	filter(PolyRank<177 | MomRank<142 | DadRank<197 | DipRank<310)




write.table(SigCompare,file="DEGsigResults.SigCompare.tsv", quote=F, sep="\t", row.names=F)



options(scipen=10)

ResCompare <- as.data.frame(Pres) %>% 
	rownames_to_column(var='both') %>%
	dplyr::arrange(desc(abs(log2FoldChange))) %>%
	mutate(PolyRank=1:nrow(Pres)) %>%
     	right_join(., category, by='both') %>% 
	dplyr::rename(target_id=both, PatMatFC=log2FoldChange, PatMatpadj=padj) %>%
     	left_join(., (as.data.frame(AsresA) %>% rownames_to_column(var='mom') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(MatMomRank=1:nrow(AsresA))), by='mom') %>% 
	dplyr::rename(MatMomFC=log2FoldChange, MatMompadj=padj) %>%
     	left_join(., (as.data.frame(AsresD) %>% rownames_to_column(var='dad') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(PatDadRank=1:nrow(AsresD))), by='dad') %>% 
	dplyr::rename(PatDadFC=log2FoldChange, PatDadpadj=padj) %>% 
     	left_join(., (as.data.frame(Bres) %>% rownames_to_column(var='target_id') %>% dplyr::arrange(desc(abs(log2FoldChange))) %>% mutate(DipRank=1:nrow(Bres))), by='target_id') %>% 
	dplyr::rename(DadMomFC=log2FoldChange,DadMompadj=padj) %>% 
	select(matches('target_id|category|Mat|Pat|DadMom|Rank|padj'))

write.table(ResCompare,file="DEGResults.ResCompare.arabidopsis.tsv", quote=F, sep="\t", row.names=F)



## make contingency tables for homoeolog comparison
Psigdf <- as.data.frame(Psig)

homoeoComp <- data.frame(category=character(),sided=character(), FC=character(0),
	MatMoreCat=integer(), PatMoreCat=integer(), 
	pvalueDE=double(), stringsAsFactors=F)
for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(Psig[row.names(Psigdf) %in% get(c) & Psigdf$log2FoldChange <0,])
		PatMoreCat <- nrow(Psig[row.names(Psigdf) %in% get(c) & Psigdf$log2FoldChange >0,])
		MatMoreNT <- nrow(Psig[row.names(Psigdf) %in% NT & Psigdf$log2FoldChange <0,])
		PatMoreNT <- nrow(Psig[row.names(Psigdf) %in% NT & Psigdf$log2FoldChange >0,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoComp[nrow(homoeoComp)+1,] <- rbind(c,s,"any",mat[1,1],mat[2,1],ft$p.value) } }
homoeoComp[nrow(homoeoComp)+1,] <- rbind("NT","NA","any",mat[1,2],mat[2,2],"NA") 


for (c in catList) { 
	for (s in sided) {
		MatMoreCat <- nrow(Psig[row.names(Psigdf) %in% get(c) & Psigdf$log2FoldChange <(-2),])
		PatMoreCat <- nrow(Psig[row.names(Psigdf) %in% get(c) & Psigdf$log2FoldChange >2,])
		MatMoreNT <- nrow(Psig[row.names(Psigdf) %in% NT & Psigdf$log2FoldChange <(-2),])
		PatMoreNT <- nrow(Psig[row.names(Psigdf) %in% NT & Psigdf$log2FoldChange >2,])
		mat <- cbind(c(MatMoreCat,PatMoreCat),c(MatMoreNT,PatMoreNT))
		ft <- fisher.test(mat,alternative=s)
		homoeoComp[nrow(homoeoComp)+1,] <- rbind(c,s,"twofold",mat[1,1],mat[2,1],ft$p.value) }}
homoeoComp[nrow(homoeoComp)+1,] <- rbind("NT","NA","twofold",mat[1,2],mat[2,2],"NA") 





HCpubtable <- homoeoComp %>% 
          mutate(sided=if_else(pvalueDE <= 0.05, paste0(sided,"_SIGNIF"), paste0(sided,"_notSig"))) %>% 
          select(-pvalueDE)

write.table(HCpubtable ,file="publicationFisherHomoeolog.Arabidopsis.tsv", row.names=F, col.names=T, sep="\t", quote=F)















### should we use lfcShrink???
PresLFC <- lfcShrink(Pdds,coef="species_AsueMom_vs_AsueDad", type="apeglm")
plotMA(PresLFC , ylim = c(-5, 5))

Pvsd <- vst(Pdds, blind=FALSE)
topVarGenesP <- head(order(rowVars(assay(Pvsd)), decreasing = TRUE), 20)
matP  <- assay(Pvsd)[ topVarGenesP, ]
matP  <- matP - rowMeans(matP)
pheatmap(t(matP))



Log fold change shrinkage for visualization and ranking
Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. To shrink the LFC, we pass the dds object to the function lfcShrink. Below we specify to use the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves on the previous estimator.

ls()
rm(list = ls()[grep("_", ls())])

