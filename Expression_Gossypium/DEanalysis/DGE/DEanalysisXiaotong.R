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
setDTthreads(12) # setting the data.table threads
options(datatable.fread.datatable=FALSE)


#### function to perform DESeq ####
ddsMake <- function (countdata, metadata) {
    ddsTmp <- DESeqDataSetFromMatrix(countData=countdata, colData=metadata, design=~species, tidy = TRUE)
    ddsTmp <- ddsTmp[rowSums(counts(ddsTmp))/nrow(metadata) >= 1,]
    ddsTmp <- estimateSizeFactors(ddsTmp)
    ddsTmp <- DESeq(ddsTmp)
    ddsTmp
}


##### ##### ##### ##### ##### ##### 
##### read in files ### ##### #####
##### ##### ##### ##### ##### ##### 

### use the above functions to aggregate files into dataframes
countdata <- countmerge("W:/corrinne/Cytonuclear/cotton/DEanalysis/counts")
tpmdata <- tpmmerge("W:/corrinne/Cytonuclear/cotton/DEanalysis/counts")

write.table(countdata, file="count.cotton.tbl", quote=F, sep="\t", row.names=F)
write.table(tpmdata, file="tpm.cotton.tbl", quote=F, row.names=F, sep="\t")

countdata <- fread("count.cotton.tbl", sep="\t")
tpmdata <- fread("tpm.cotton.tbl", sep="\t")

### generate metadata framework
metadata <- data.frame(id=names(countdata[,-1]), data.frame(id=names(countdata[,-1])) %>% separate(id, c("species", "rep", NA)))
countdata[,-1] <- round(countdata[,-1],0)


### generate df and metadata for subgenomes
category <- read.table("gene.pairs", col.names=c("category","dad","mom")) %>% 
     filter(dad %in% countdata$target_id) %>% filter(mom %in% countdata$target_id)
category$both <- paste0(category$dad,"_",category$mom)


Acountdata <- countdata %>% filter(grepl("Gohir.A",target_id)) %>% select(!contains("D5"))
Dcountdata <- countdata %>% filter(grepl("Gohir.D",target_id)) %>% select(!contains("A2"))

Ametadata <- filter(metadata, !grepl("D5",species))
Dmetadata <- filter(metadata, !grepl("A2",species)) 



#### get the count sum of homoeolog pairs ####
homoeoPairs <- category[,2:3]
cpmt <- NULL
Scountdata <- countdata[countdata$target_id %in% cpmt,]

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

resultsNames(dds)
res <- results(dds)
summary(res) 


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




##### ##### ##### ##### ##### ##### 
##### Paternal (D) comparison ##### 
##### ##### ##### ##### ##### ##### 

Ddds <- ddsMake(Dcountdata,Dmetadata)

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









ls()
rm(list = ls()[grep("_", ls())])

