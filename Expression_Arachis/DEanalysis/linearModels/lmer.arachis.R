setwd("Y:/corrinne/Cytonuclear/arachis/DEanalysis/linearModels")

library(tidyverse)
library(lme4)
library(car)
library(sjPlot)
library(rstatix)
library(emmeans)
library(data.table)
library(gridExtra)

##### global options #####
emm_options(lmerTest.limit = 100000)
setDTthreads(10)
options(datatable.fread.datatable=FALSE)


# read and filter the table
rlogfull <- fread("arachis.lm.tbl", header=T, sep="\t") %>% 
	filter(str_detect(gene,"arahy")) %>%
	filter(!is.na(subgenome)) %>%
	filter(!is.na(Aipadur_1)) %>% 
	rename(Subgenome=subgenome) %>%
	rename(Targeting=category)%>%
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-"))%>%
	select(!contains("CyMira")) # put NOT first alphabetically

rlog <- rlogfull %>% 
	select(!contains(c("Adur_","Aipa_"),ignore.case=F))

# read in and number the pairs then extract into two 2-column df
pairs <- fread("arachis.gene.pairs", col.names=c("category","dad","mom"), sep="\t") %>% rowid_to_column("genePairID")
one <- pairs[, c(1, 3)]
two <- pairs[, c(1, 4)]
names(one) <- c("genePairID", "gene")
names(two) <- c("genePairID", "gene")


# join each df with the rlog df to assign genePairID to gene
data <- bind_rows(one,two) %>%
	right_join(.,rlog,by='gene') %>%
	pivot_longer(cols=starts_with("A"), names_to="plantID", values_to="rlog")


# separate Adur and Aipadur, and fix Aipadur
# parents are backwards in Aipadur, so we need to reverse designations
Ahdata <- data %>% filter(str_detect(plantID,"Ahyp"))
Aiddata <- data %>% filter(str_detect(plantID,"Aipadur")) %>%
	mutate(Subgenome = str_replace(Subgenome,"mom","pat")) %>%
	mutate(Subgenome = str_replace(Subgenome,"dad","mom")) %>% 
	mutate(Subgenome = str_replace(Subgenome,"pat","dad"))


##### run model 1, Ah #####

## run mixed model
Ahm1 <- lmer(rlog ~ Subgenome * Targeting + (1 | genePairID) + (1 | plantID), data=Ahdata, REML=T)
summary(Ahm1)

## make a table of fixed effects
Ahm1.fixedEffects <- summary(Ahm1)$coefficients
write.table(Ahm1.fixedEffects, file="Ahm1.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Ahm1, type=2)
Ahm1.anova <- anova_summary(Anova(Ahm1, type=3))
fwrite(Ahm1.anova, file="Ahm1.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## graph fixed effects with CI
Ahm1.emmp <- emmip(Ahm1, Subgenome ~ Targeting, plotit=F, CIs=T)
Ahm1.fixedPlot <- emmip_ggplot(Ahm1.emmp, style = "factor", dodge = 0.1, 
  facetlab = "label_context", CIarg = list(lwd = 2, alpha = 0.5),
  PIarg = list(lwd = 1.25, alpha = 0.33)) + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("Ahm1.fixedEffects.jpg", device="jpeg", units="in", height=8, width=12, plot=Ahm1.fixedPlot)

## perform contrasts and generate table
Ahm1.emm <- emmeans(Ahm1, pairwise ~ Subgenome*Targeting, lmer.df='satterthwaite',adjust='Holm')
Ahm1.contrasts <- as.data.frame(Ahm1.emm$contrasts)

Ahm1.tbl <- Ahm1.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Parent1","Category1","Parent2","Category2"),sep=" ") %>%
	filter( ((Parent1 == Parent2 & Category1 != Category2) | (Parent1 != Parent2 & Category1 == Category2)) & p.value<0.05)

fwrite(Ahm1.tbl, file="Ahm1.SigContrasts.tsv", quote=F, row.names=F, sep="\t")



##### run model 1, Aid #####

## run mixed model
Aidm1 <- lmer(rlog ~ Subgenome * Targeting + (1 | genePairID) + (1 | plantID), data=Aiddata, REML=T)
# note: model fails to converge
summary(Aidm1)

## make a table of fixed effects
Aidm1.fixedEffects <- summary(Aidm1)$coefficients
write.table(Aidm1.fixedEffects, file="Aidm1.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Aidm1, type=2)
Aidm1.anova <- anova_summary(Anova(Aidm1, type=3))
fwrite(Aidm1.anova, file="Aidm1.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## graph fixed effects with CI
Aidm1.emmp <- emmip(Aidm1, Subgenome ~ Targeting, plotit=F, CIs=T)
Aidm1.fixedPlot <- emmip_ggplot(Aidm1.emmp, style = "factor", dodge = 0.1, 
  facetlab = "label_context", CIarg = list(lwd = 2, alpha = 0.5),
  PIarg = list(lwd = 1.25, alpha = 0.33)) + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("Aidm1.fixedEffects.jpg", device="jpeg", units="in", height=8, width=12, plot=Aidm1.fixedPlot)

## perform contrasts and generate table
Aidm1.emm <- emmeans(Aidm1, pairwise ~ Subgenome*Targeting, lmer.df='satterthwaite',adjust='Holm')
Aidm1.contrasts <- as.data.frame(Aidm1.emm$contrasts)

Aidm1.tbl <- Aidm1.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Parent1","Category1","Parent2","Category2"),sep=" ") %>%
	filter( ((Parent1 == Parent2 & Category1 != Category2) | (Parent1 != Parent2 & Category1 == Category2)) & p.value<0.05)

fwrite(Aidm1.tbl, file="Aidm1.SigContrasts.tsv", quote=F, row.names=F, sep="\t")




#####
# genePairID is not a good random predictor... too many df and not enough independent observations
# note the Inf degrees freedom
# make new model rlogRatio ~ Targeting + (1 | plantID)
# or other new model delta_rlog ~ Targeting + (1 | plantID)
#####

## make a df that has diploids combined
rlogmomAh <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","Adur"),ignore.case=F)) %>% 
	setNames(gsub("Adur", "Dip", names(.)))

rlogdadAh <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","Aipa_"))) %>% 
	setNames(gsub("Aipa_", "Dip_", names(.)))

rlogdipAh <- bind_rows(rlogmomAh,rlogdadAh) %>% 
	right_join(.,rlog,by='gene') %>%
	select(!(contains(c("subgenome","Aipadur")))) 



# parents are backwards in Aipadur, so we need to reverse designations
rlogmomAid <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","Aipa_"))) %>% 
	setNames(gsub("Aipa_", "Dip_", names(.)))

rlogdadAid <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","Adur"),ignore.case=F)) %>% 
	setNames(gsub("Adur", "Dip", names(.)))

rlogdipAid <- bind_rows(rlogmomAid,rlogdadAid) %>% 
	right_join(.,rlog,by='gene') %>%
	select(!(contains(c("subgenome","Ahyp")))) 




##### count difference #####
# for uniformity, always mom - dad or mom/dad
# compared to momDip - dadDip or momDip/dadDip
homoeoPairsAh <- pairs[,3:4] %>% filter(dad %in% rlog$gene & mom %in% rlog$gene)
homoeoPairsAid <- pairs[,3:4] %>% filter(dad %in% rlog$gene & mom %in% rlog$gene) %>% 
	relocate(dad, .after=mom) %>%
	rename(mom = dad, dad = mom)
	


### make df for Ah
## make the delta
DrlogAh <- rlogdipAh[1:2,]
DrlogAh <- DrlogAh[-(1:2),]

for (i in 1:nrow(homoeoPairsAh)) {
    DrlogAh[nrow(DrlogAh)+1,] <- rbind(
        cbind(gene=paste0(homoeoPairsAh[i,2],"_",homoeoPairsAh[i,1]),
        rlogdipAh[rlogdipAh$gene==homoeoPairsAh[i,2], -c(1,11)] - 
        rlogdipAh[rlogdipAh$gene==homoeoPairsAh[i,1], -c(1,11)],
        rlogdipAh[rlogdipAh$gene==homoeoPairsAh[i,1], c(11)]))
}

DdataAh <- DrlogAh %>% 
	select(!contains(c("Dip_1"))) %>%
	unite("Plant2", c(Dip_2,Ahyp_2), sep="X") %>% 
	unite("Plant3", c(Dip_3,Ahyp_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,Ahyp_4), sep="X") %>% 
	unite("Plant5", c(Dip_5,Ahyp_5), sep="X") %>% 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipDelta","tetDelta"), sep="X") %>%
	mutate(across(ends_with("Delta"), as.double))
                                                                 
fwrite(DdataAh, file="Ah.delta.tsv", quote=F, row.names=F, sep="\t")



## make the ratio for Ah
RrlogAh <- rlogdipAh[1:2,]
RrlogAh <- RrlogAh[-(1:2),]

for (i in 1:nrow(homoeoPairsAh)) {
    RrlogAh[nrow(RrlogAh)+1,] <- rbind(
        cbind(gene=paste0(homoeoPairsAh[i,2],"_",homoeoPairsAh[i,1]),
        rlogdipAh[rlogdipAh$gene==homoeoPairsAh[i,2], -c(1,11)] / 
        rlogdipAh[rlogdipAh$gene==homoeoPairsAh[i,1], -c(1,11)],
        rlogdipAh[rlogdipAh$gene==homoeoPairsAh[i,1], c(11)]))
}


RdataAh <- RrlogAh %>% 
	select(!contains(c("Dip_1"))) %>%
	unite("Plant2", c(Dip_2,Ahyp_2), sep="X") %>% 
	unite("Plant3", c(Dip_3,Ahyp_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,Ahyp_4), sep="X") %>% 
	unite("Plant5", c(Dip_5,Ahyp_5), sep="X") %>% 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipRatio","tetRatio"), sep="X") %>%
	mutate(across(ends_with("Ratio"), as.double))

fwrite(RdataAh, file="Ah.ratio.tsv", quote=F, row.names=F, sep="\t")



### make df for Aid
## make the delta
DrlogAid <- rlogdipAid[1:2,]
DrlogAid <- DrlogAid[-(1:2),]

for (i in 1:nrow(homoeoPairsAid)) {
    DrlogAid[nrow(DrlogAid)+1,] <- rbind(
        cbind(gene=paste0(homoeoPairsAid[i,2],"_",homoeoPairsAid[i,1]),
        rlogdipAid[rlogdipAid$gene==homoeoPairsAid[i,2], -c(1,11)] - 
        rlogdipAid[rlogdipAid$gene==homoeoPairsAid[i,1], -c(1,11)],
        rlogdipAid[rlogdipAid$gene==homoeoPairsAid[i,1], c(11)]))
}

DdataAid <- DrlogAid %>% 
	select(!contains(c("Dip_5"))) %>%
	unite("Plant1", c(Dip_1,Aipadur_1), sep="X") %>% 
	unite("Plant2", c(Dip_2,Aipadur_2), sep="X") %>% 
	unite("Plant3", c(Dip_3,Aipadur_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,Aipadur_4), sep="X") %>% 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipDelta","tetDelta"), sep="X") %>%
	mutate(across(ends_with("Delta"), as.double))

fwrite(DdataAid, file="Aid.delta.tsv", quote=F, row.names=F, sep="\t")



## make the ratio for Aid
RrlogAid <- rlogdipAid[1:2,]
RrlogAid <- RrlogAid[-(1:2),]

for (i in 1:nrow(homoeoPairsAid)) {
    RrlogAid[nrow(RrlogAid)+1,] <- rbind(
        cbind(gene=paste0(homoeoPairsAid[i,2],"_",homoeoPairsAid[i,1]),
        rlogdipAid[rlogdipAid$gene==homoeoPairsAid[i,2], -c(1,11)] / 
        rlogdipAid[rlogdipAid$gene==homoeoPairsAid[i,1], -c(1,11)],
        rlogdipAid[rlogdipAid$gene==homoeoPairsAid[i,1], c(11)]))
}


RdataAid <- RrlogAid %>% 
	select(!contains(c("Dip_5"))) %>%
	unite("Plant1", c(Dip_1,Aipadur_1), sep="X") %>% 
	unite("Plant2", c(Dip_2,Aipadur_2), sep="X") %>% 
	unite("Plant3", c(Dip_3,Aipadur_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,Aipadur_4), sep="X") %>% 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipRatio","tetRatio"), sep="X") %>%
	mutate(across(ends_with("Ratio"), as.double))


fwrite(RdataAid, file="Aid.ratio.tsv", quote=F, row.names=F, sep="\t")







##### run model 2, Ah #####

### run mixed model, ratio ###
Ahm2R <- lmer(tetRatio ~ Targeting + (1 | plantID), data=RdataAh, REML=T)
summary(Ahm2R)
# boundary (singular) fit: see ?isSingular

Ahm2R <- lm(tetRatio ~ Targeting, data=RdataAh)
summary(Ahm2R)

## save the model plot
jpeg(file="Ahm2.Ratio.modelPlot.jpeg")
plot(Ahm2R)
dev.off()

## make a table of fixed effects
Ahm2R.fixedEffects <- as.data.frame(summary(Ahm2R)$coefficients)
fwrite(Ahm2R.fixedEffects, file="Ahm2.Ratio.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Ahm2R, type=2)
Ahm2R.anova <- anova_summary(Anova(Ahm2R, type=3))
fwrite(Ahm2R.anova, file="Ahm2.Ratio.anovaSummary.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Ahm2R.emm <- emmeans(Ahm2R, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Ahm2R.contrasts <- as.data.frame(Ahm2R.emm$contrasts)

Ahm2R.tbl <- Ahm2R.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") %>%
	filter(p.value<0.05)

fwrite(Ahm2R.tbl, file="Ahm2.Ratio.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


### run mixed model, delta ###
Ahm2D <- lmer(tetDelta ~ Targeting + (1 | plantID), data=DdataAh, REML=T)
summary(Ahm2D)
# boundary (singular) fit: see ?isSingular

Ahm2D <- lm(tetDelta ~ Targeting, data=DdataAh)
summary(Ahm2D)


## save the model plot
jpeg(file="Ahm2.Delta.modelPlot.jpeg")
plot(Ahm2D)
dev.off()

## make a table of fixed effects
Ahm2D.fixedEffects <- as.data.frame(summary(Ahm2D)$coefficients)
fwrite(Ahm2D.fixedEffects, file="Ahm2.Delta.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Ahm2D, type=2)
Ahm2D.anova <- anova_summary(Anova(Ahm2D, type=3))
fwrite(Ahm2D.anova, file="Ahm2.Delta.anovaSummary.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Ahm2D.emm <- emmeans(Ahm2D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Ahm2D.contrasts <- as.data.frame(Ahm2D.emm$contrasts)

Ahm2D.tbl <- Ahm2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ")  %>%
	filter( p.value<0.05)

fwrite(Ahm2D.tbl, file="Ahm2.Delta.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


Ahm2D.all.tbl <- Ahm2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") 

fwrite(Ahm2D.all.tbl, file="Ahm2.Delta.AllContrasts.tsv", quote=F, row.names=F, sep="\t")


##### run model 2, Aid #####

### run mixed model, ratio ###
Aidm2R <- lmer(tetRatio ~ Targeting + (1 | plantID), data=RdataAid, REML=T)
summary(Aidm2R)

Aidm2R <- lm(tetRatio ~ Targeting, data=RdataAid)
summary(Aidm2R)

## save the model plot
jpeg(file="Aidm2.Ratio.modelPlot.jpeg")
plot(Aidm2R)
dev.off()

## make a table of fixed effects
Aidm2R.fixedEffects <- as.data.table(summary(Aidm2R)$coefficients)
fwrite(Aidm2R.fixedEffects, file="Aidm2.Ratio.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Aidm2R, type=2)
Aidm2R.anova <- anova_summary(Anova(Aidm2R, type=3))
fwrite(Aidm2R.anova, file="Aidm2.Ratio.anovaSummary.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Aidm2R.emm <- emmeans(Aidm2R, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Aidm2R.contrasts <- as.data.frame(Aidm2R.emm$contrasts)

Aidm2R.tbl <- Aidm2R.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") %>%
	filter(p.value<0.05)

fwrite(Aidm2R.tbl, file="Aidm2.Ratio.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


### run mixed model, Aid delta ###
Aidm2D <- lmer(tetDelta ~ Targeting + (1 | plantID), data=DdataAid, REML=T)
summary(Aidm2D)

Aidm2D <- lm(tetDelta ~ Targeting, data=DdataAid)
summary(Aidm2D)

## save the model plot
jpeg(file="Aidm2.Delta.modelPlot.jpeg")
plot(Aidm2D)
dev.off()

## make a table of fixed effects
Aidm2D.fixedEffects <- as.data.frame(summary(Aidm2D)$coefficients)
fwrite(Aidm2D.fixedEffects, file="Aidm2.Delta.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Aidm2D, type=2)
Aidm2D.anova <- anova_summary(Anova(Aidm2D, type=3))
fwrite(Aidm2D.anova, file="Aidm2.Delta.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## perform contrasts and generate table
Aidm2D.emm <- emmeans(Aidm2D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Aidm2D.contrasts <- as.data.frame(Aidm2D.emm$contrasts)

Aidm2D.tbl <- Aidm2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") %>%
	filter( p.value<0.05)

fwrite(Aidm2D.tbl, file="Aidm2.Delta.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


Aidm2D.all.tbl <- Aidm2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") 

fwrite(Aidm2D.all.tbl, file="Aidm2.Delta.AllContrasts.tsv", quote=F, row.names=F, sep="\t")



##### run model 3, Ah #####

### run mixed model 3, ratio ###
Ahm3R <- lmer(tetRatio ~ Targeting*dipRatio + (1 | plantID), data=RdataAh, REML=T)
summary(Ahm3R)
#boundary (singular) fit: see ?isSingular

Ahm3R <- lm(tetRatio ~ Targeting*dipRatio, data=RdataAh)
summary(Ahm3R)

## save the model plot
jpeg(file="Ahm3.Ratio.modelPlot.jpeg")
plot(Ahm3R)
dev.off()

## make a table of fixed effects
Ahm3R.fixedEffects <- as.data.frame(summary(Ahm3R)$coefficients)
fwrite(Ahm3R.fixedEffects, file="Ahm3.Ratio.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Ahm3R, type=2)
Ahm3R.anova <- anova_summary(Anova(Ahm3R, type=3))
fwrite(Ahm3R.anova, file="Ahm3.Ratio.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## graph fixed effects
emmip(Ahm3R, dipRatio ~ Targeting, cov.reduce=F)

## perform contrasts and generate table
Ahm3R.emm <- emmeans(Ahm3R, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Ahm3R.contrasts <- as.data.frame(Ahm3R.emm$contrasts)

Ahm3R.tbl <- Ahm3R.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") %>%
	filter( p.value<0.05)

fwrite(Ahm3R.tbl, file="Ahm3.Ratio.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


### run mixed model 3, delta ###
Ahm3D <- lmer(tetDelta ~ Targeting*dipDelta + (1 | plantID), data=DdataAh, REML=T)
summary(Ahm3D)

Ahm3D <- lm(tetDelta ~ Targeting*dipDelta, data=DdataAh)
summary(Ahm3D)



## save the model plot
jpeg(file="Ahm3.Delta.modelPlot.jpeg")
plot(Ahm3D)
dev.off()

## make a table of fixed effects
Ahm3D.fixedEffects <- as.data.frame(summary(Ahm3D)$coefficients)
fwrite(Ahm3D.fixedEffects, file="Ahm3.Delta.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Ahm3D, type=2)
Ahm3D.anova <- anova_summary(Anova(Ahm3D, type=3))
fwrite(Ahm3D.anova, file="Ahm3.Delta.anovaSummary.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Ahm3D.emm <- emmeans(Ahm3D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
# NOTE: Results may be misleading due to involvement in interactions
Ahm3D.contrasts <- as.data.frame(Ahm3D.emm$contrasts)

Ahm3D.tbl <- Ahm3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") %>%
	filter( p.value<0.05)

fwrite(Ahm3D.tbl, file="Ahm3.Delta.SigContrasts.tsv", quote=F, row.names=F, sep="\t")

Ahm3D.all.tbl <- Ahm3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") 

fwrite(Ahm3D.all.tbl, file="Ahm3.Delta.AllContrasts.tsv", quote=F, row.names=F, sep="\t")






##### run model 3, Aid #####

### run mixed model 3, ratio ###
Aidm3R <- lmer(tetRatio ~ Targeting*dipRatio + (1 | plantID), data=RdataAid, REML=T)
summary(Aidm3R)

Aidm3R <- lm(tetRatio ~ Targeting*dipRatio, data=RdataAid)
summary(Aidm3R)

## save the model plot
jpeg(file="Aidm3.Ratio.modelPlot.jpeg")
plot(Aidm3R)
dev.off()

## make a table of fixed effects
Aidm3R.fixedEffects <- as.data.frame(summary(Aidm3R)$coefficients)
fwrite(Aidm3R.fixedEffects, file="Aidm3.Ratio.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Aidm3R, type=2)
Aidm3R.anova <- anova_summary(Anova(Aidm3R, type=3))
fwrite(Aidm3R.anova, file="Aidm3.Ratio.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## perform contrasts and generate table
Aidm3R.emm <- emmeans(Aidm3R, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
# NOTE: Results may be misleading due to involvement in interactions
Aidm3R.contrasts <- as.data.frame(Aidm3R.emm$contrasts)

Aidm3R.tbl <- Aidm3R.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") %>%
	filter( p.value<0.05)

fwrite(Aidm3R.tbl, file="Aidm3.Ratio.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


### run mixed model 3, delta ###
Aidm3D <- lmer(tetDelta ~ Targeting*dipDelta + (1 | plantID), data=DdataAid, REML=T)
summary(Aidm3D)

Aidm3D <- lm(tetDelta ~ Targeting*dipDelta, data=DdataAid)
summary(Aidm3D)


## save the model plot
jpeg(file="Aidm3.Delta.modelPlot.jpeg")
plot(Aidm3D)
dev.off()

## make a table of fixed effects
Aidm3D.fixedEffects <- as.data.frame(summary(Aidm3D)$coefficients)
fwrite(Aidm3D.fixedEffects, file="Aidm3.Delta.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Aidm3D, type=2)
Aidm3D.anova <- anova_summary(Anova(Aidm3D, type=3))
fwrite(Aidm3D.anova, file="Aidm3.Delta.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## perform contrasts and generate table
Aidm3D.emm <- emmeans(Aidm3D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
# NOTE: Results may be misleading due to involvement in interactions
Aidm3D.contrasts <- as.data.frame(Aidm3D.emm$contrasts)

Aidm3D.tbl <- Aidm3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") %>%
	filter( p.value<0.05)

fwrite(Aidm3D.tbl, file="Aidm3.Delta.SigContrasts.tsv", quote=F, row.names=F, sep="\t")



Aidm3D.all.tbl <- Aidm3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") 

fwrite(Aidm3D.all.tbl, file="Aidm3.Delta.AllContrasts.tsv", quote=F, row.names=F, sep="\t")











#####
# it weirds me out that delta and ratio give different results
# combining mom & dad doesn't speak to the subgenomes directly
# let's model tetMom, tetDad separately together for Model 4
#####


SSrlogmomAh <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","Adur","Ahyp"),ignore.case=F)) %>% 
	setNames(gsub("Adur", "dipMom", names(.))) %>% 
	setNames(gsub("Ahyp", "tetMom", names(.))) %>% 
	rename(momGene = gene)

SSrlogdadAh <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","Aipa_","Ahyp"))) %>% 
	setNames(gsub("Aipa", "dipDad", names(.))) %>% 
	setNames(gsub("Ahyp", "tetDad", names(.))) %>% 
	rename(dadGene = gene)


SSrlogAh <- pairs %>% select(!genePairID) %>%
	rename(Targeting = category) %>%
	rename(dadGene = dad) %>%
	rename(momGene = mom)  %>% 
	left_join(.,SSrlogmomAh,by='momGene') %>% 
	left_join(.,SSrlogdadAh,by='dadGene') %>%
	filter(!(is.na(dipMom_1) | is.na(dipDad_1))) %>%
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-"))



SSdataAh <- SSrlogAh %>% 
	select(!contains(c("dipMom_5","dipDad_5"))) %>% 
	unite("PlantA", c(dipMom_1, dipDad_1, tetMom_2, tetDad_2), sep="X") %>% 
	unite("PlantB", c(dipMom_2, dipDad_2, tetMom_3, tetDad_3), sep="X") %>% 
	unite("PlantC", c(dipMom_3, dipDad_3, tetMom_4, tetDad_4), sep="X") %>% 
	unite("PlantD", c(dipMom_4, dipDad_4, tetMom_5, tetDad_5), sep="X") %>% 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipMom","dipDad","tetMom","tetDad"), sep="X") %>% 
	mutate(across(ends_with("Mom"), as.double)) %>% 
	mutate(across(ends_with("Dad"), as.double))




# Aid is opposite parentage, so here are the mental gymnastics
SSrlogmomAid <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","Aipa"),ignore.case=F)) %>% 
	setNames(gsub("Aipa_", "dipMom_", names(.))) %>% 
	setNames(gsub("Aipadur", "tetMom", names(.))) %>% 
	rename(momGene = gene)

SSrlogdadAid <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","Aipadur","Adur"))) %>% 
	setNames(gsub("Aipadur", "tetDad", names(.))) %>% 
	setNames(gsub("Adur", "dipDad", names(.))) %>% 
	rename(dadGene = gene)


SSrlogAid <- pairs %>% select(!genePairID) %>%
	rename(Targeting = category) %>%
	rename(dadGene = mom) %>%
	rename(momGene = dad)  %>% 
	left_join(.,SSrlogmomAid,by='momGene') %>% 
	left_join(.,SSrlogdadAid,by='dadGene') %>%
	filter(!(is.na(dipMom_1) | is.na(dipDad_1))) %>%
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-"))



SSdataAid <- SSrlogAid %>% 
	select(!contains(c("dipMom_5","dipDad_5"))) %>% 
	unite("PlantA", c(dipMom_1, dipDad_1, tetMom_1, tetDad_1), sep="X") %>% 
	unite("PlantB", c(dipMom_2, dipDad_2, tetMom_2, tetDad_2), sep="X") %>% 
	unite("PlantC", c(dipMom_3, dipDad_3, tetMom_3, tetDad_3), sep="X") %>% 
	unite("PlantD", c(dipMom_4, dipDad_4, tetMom_4, tetDad_4), sep="X") %>% 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipMom","dipDad","tetMom","tetDad"), sep="X") %>% 
	mutate(across(ends_with("Mom"), as.double)) %>% 
	mutate(across(ends_with("Dad"), as.double))



fwrite(SSdataAh, file="Ah.SS.tsv", quote=F, row.names=F, sep="\t")
fwrite(SSdataAid, file="Aid.SS.tsv", quote=F, row.names=F, sep="\t")

# check data
plot(SSdataAh$tetDad ~ SSdataAh$tetMom)
plot(SSdataAid$tetDad ~ SSdataAid$tetMom)


SSplotAh <- SSdataAh %>% 
	select(Targeting, tetMom, tetDad) %>%
	mutate(Targeting = str_replace(Targeting,"-targeted_Non-interacting","")) %>%
	mutate(Targeting = str_replace(Targeting,"-targeted_Interacting","")) %>%
	mutate(Targeting = str_replace(Targeting,"-organelle-targeted","")) %>% 
	arrange(Targeting) %>%
	ggplot(aes(x=tetMom, y=tetDad, color=Targeting)) + 
	   geom_point(size=0.9) + 
	   scale_color_manual(values=c('lightgrey','red3','blue4', 'green4')) +
	   geom_abline(intercept = 0, slope = 1) +
	   ggtitle("Gossypium hirsutum") + 
	   theme(plot.title = element_text(face = "italic")) + 
	   labs(x="Maternal Homoeolog Expression (rlog)", y="Paternal Homoeolog Expression (rlog)")

ggsave("AhMomvsDad.jpg",SSplotAh)


SSplotAid <- SSdataAid %>% 
	select(Targeting, tetMom, tetDad) %>%
	mutate(Targeting = str_replace(Targeting,"-targeted_Non-interacting","")) %>%
	mutate(Targeting = str_replace(Targeting,"-targeted_Interacting","")) %>%
	mutate(Targeting = str_replace(Targeting,"-organelle-targeted","")) %>% 
	arrange(Targeting) %>%
	ggplot(aes(x=tetMom, y=tetDad, color=Targeting)) + 
	   geom_point(size=0.9) + 
	   scale_color_manual(values=c('lightgrey','red3','blue4', 'green4')) +
	   geom_abline(intercept = 0, slope = 1) +
	   ggtitle("Arachis ipaensis x Arachis duranensis") + 
	   theme(plot.title = element_text(face = "italic")) + 
	   labs(x="Maternal Homoeolog Expression (rlog)", y="Paternal Homoeolog Expression (rlog)")

ggsave("AidMomvsDad.jpg",SSplotAid)





##### run model 4, Ah #####

Ahm4 <- lm(cbind(tetMom, tetDad) ~ Targeting*dipMom*dipDad, data = SSdataAh)
summary(Ahm4)

momPlotAh <- emmip_ggplot(emmip(Ahm4, Targeting ~ dipMom | subgenome, mult.name = "subgenome", at=list(dipMom=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
dadPlotAh <- emmip_ggplot(emmip(Ahm4, Targeting ~ dipDad | subgenome, mult.name = "subgenome", at=list(dipDad=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))

momdadPlotAh <- grid.arrange(momPlotAh,dadPlotAh,ncol=1)
ggsave("Ahm4.linearPrediction.jpg",momdadPlotAh, height=8, width=8)


# also run the models independently to access results easier
# will give same results
Ahm4Mom <- lm(tetMom ~ Targeting*dipMom*dipDad, data = SSdataAh)
Ahm4Dad <- lm(tetDad ~ Targeting*dipMom*dipDad, data = SSdataAh)


# check residual normality
SSdataAh$predictMom <- predict(Ahm4)[,1]
SSdataAh$predictDad <- predict(Ahm4)[,2]

SSdataAh$residMom <- SSdataAh$tetMom - SSdataAh$predictMom
SSdataAh$residDad <- SSdataAh$tetDad - SSdataAh$predictDad

SSdataAh$zresidMom <- (SSdataAh$residMom - mean(SSdataAh$residMom))/sd(SSdataAh$residMom)
SSdataAh$zresidDad <- (SSdataAh$residDad - mean(SSdataAh$residDad))/sd(SSdataAh$residDad)

histMomAh <- hist(SSdataAh$zresidMom, plot=F)
histDadAh <- hist(SSdataAh$zresidDad, plot=F)

cDad <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
cMom <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

#jpeg(file="Ahm4.modelPlot.normality.jpeg")
plot(histMomAh, col=cMom)
plot(histDadAh, col=cDad, add=T)
#dev.off()


# heteroscedasticity
plot(SSdataAh$zresidMom ~ as.factor(SSdataAh$Targeting), col=cMom)
plot(SSdataAh$zresidDad ~ as.factor(SSdataAh$Targeting), col=cDad, add=T)

plot(SSdataAh$zresidMom ~ SSdataAh$tetDad, col=cMom)
plot(SSdataAh$zresidDad ~ SSdataAh$tetMom, col=cDad)



## make a table of fixed effects
Ahm4.fixedEffects <- as.data.frame(tidy(Ahm4) %>%
	mutate(term= str_replace(term,"_Non-interacting","NI")) %>% 
	mutate(term= str_replace(term,"Targeting","")) %>% 
	mutate(term= str_replace(term,"Dual-targeted","D")) %>% 
	mutate(term= str_replace(term,"Mitochondria-targeted","M")) %>% 
	mutate(term= str_replace(term,"Plastid-targeted","P")) %>% 
	mutate(term= str_replace(term,"_Interacting","I"))
)

write.table(Ahm4.fixedEffects, file="Ahm4.fixedEffects.tsv", quote=F, row.names=F, sep="\t") 

Ahm4.fixedEffects.sig <- Ahm4.fixedEffects %>%
	filter(p.value<0.05)



## run a type III anova and make table of summary
Ahm4.anova <- Anova(Ahm4, type=3)
Ahm4Mom.anova <- anova_summary(Anova(Ahm4Mom, type=3))
Ahm4Dad.anova <- anova_summary(Anova(Ahm4Dad, type=3))

write.table(Ahm4Mom.anova, file="Ahm4.anovaSummary.Mom.tsv", quote=F, row.names=F, sep="\t")
write.table(Ahm4Dad.anova, file="Ahm4.anovaSummary.Dad.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Ahm4Mom.emm <- emmeans(Ahm4Mom, pairwise ~ Targeting*dipMom*dipDad, lmer.df='satterthwaite',adjust='Holm')
Ahm4Mom.contrasts <- as.data.frame(Ahm4Mom.emm$contrasts)
Ahm4Mom.contrasts$model <- "Ahm4Mom"

Ahm4Dad.emm <- emmeans(Ahm4Dad, pairwise ~ Targeting*dipMom*dipDad, lmer.df='satterthwaite',adjust='Holm')
Ahm4Dad.contrasts <- as.data.frame(Ahm4Dad.emm$contrasts)
Ahm4Dad.contrasts$model <- "Ahm4Dad"

Ahm4.emm <- emmeans(Ahm4, pairwise ~ Targeting*dipMom*dipDad, lmer.df='satterthwaite',adjust='Holm')
Ahm4.contrasts <- as.data.frame(Ahm4.emm$contrasts)
Ahm4.contrasts$model <- "Ahm4full"

Ahm4.tbl <- rbind(Ahm4.contrasts,Ahm4Mom.contrasts,Ahm4Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,NA,"Category2",NA,NA),sep=" ") %>%
	filter(p.value<0.05)

write.table(Ahm4.tbl, file="Ahm4.SigContrasts.tsv", quote=F, row.names=F, sep="\t")



##### run model 4, Aid #####

Aidm4 <- lm(cbind(tetMom, tetDad) ~ Targeting*dipMom*dipDad, data = SSdataAid)
summary(Aidm4)

momPlotAid <- emmip_ggplot(emmip(Aidm4, Targeting ~ dipMom | subgenome, mult.name = "subgenome", at=list(dipMom=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
dadPlotAid <- emmip_ggplot(emmip(Aidm4, Targeting ~ dipDad | subgenome, mult.name = "subgenome", at=list(dipDad=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))

momdadPlotAid <- grid.arrange(momPlotAid,dadPlotAid,ncol=1)
ggsave("Aidm4.linearPrediction.jpg",momdadPlotAid, height=8, width=10)


# also run the models independently to access results easier
# will give same results
Aidm4Mom <- lm(tetMom ~ Targeting*dipMom*dipDad, data = SSdataAid)
Aidm4Dad <- lm(tetDad ~ Targeting*dipMom*dipDad, data = SSdataAid)


# check residual normality
SSdataAid$predictMom <- predict(Aidm4)[,1]
SSdataAid$predictDad <- predict(Aidm4)[,2]

SSdataAid$residMom <- SSdataAid$tetMom - SSdataAid$predictMom
SSdataAid$residDad <- SSdataAid$tetDad - SSdataAid$predictDad

SSdataAid$zresidMom <- (SSdataAid$residMom - mean(SSdataAid$residMom))/sd(SSdataAid$residMom)
SSdataAid$zresidDad <- (SSdataAid$residDad - mean(SSdataAid$residDad))/sd(SSdataAid$residDad)

histMomAid <- hist(SSdataAid$zresidMom, plot=F)
histDadAid <- hist(SSdataAid$zresidDad, plot=F)

# residuals are looking weird....

plot(histMomAid, col=cMom)
plot(histDadAid, col=cDad, add=T)



AidHist <- ggplot(data=SSdataAid) +
	geom_histogram(aes(x=zresidMom),alpha=0.7,fill="red",bins=70) + 
	geom_histogram(aes(x=zresidDad),alpha=0.7,fill="blue",bins=70)

ggsave("Aidm4.modelPlot.normality.jpeg",AidHist)






# heteroscedasticity; mom is high, dad is low
plot(SSdataAid$zresidMom ~ as.factor(SSdataAid$Targeting), col=cMom)
plot(SSdataAid$zresidDad ~ as.factor(SSdataAid$Targeting), col=cDad, add=T)

# mom has a weird shape here
plot(SSdataAid$zresidMom ~ SSdataAid$tetDad, col=cMom)
plot(SSdataAid$zresidDad ~ SSdataAid$tetMom, col=cDad)



## make a table of fixed effects
Aidm4.fixedEffects <- as.data.frame(tidy(Aidm4) %>%
	mutate(term= str_replace(term,"_Non-interacting","NI")) %>% 
	mutate(term= str_replace(term,"Targeting","")) %>% 
	mutate(term= str_replace(term,"Dual-targeted","D")) %>% 
	mutate(term= str_replace(term,"Mitochondria-targeted","M")) %>% 
	mutate(term= str_replace(term,"Plastid-targeted","P")) %>% 
	mutate(term= str_replace(term,"_Interacting","I"))
)

write.table(Aidm4.fixedEffects, file="Aidm4.fixedEffects.tsv", quote=F, row.names=F, sep="\t") 

Aidm4.fixedEffects.sig <- Aidm4.fixedEffects %>%
	filter(p.value<0.05)



## run a type III anova and make table of summary
Aidm4.anova <- Anova(Aidm4, type=3)
Aidm4Mom.anova <- anova_summary(Anova(Aidm4Mom, type=3))
Aidm4Dad.anova <- anova_summary(Anova(Aidm4Dad, type=3))

write.table(Aidm4Mom.anova, file="Aidm4.anovaSummary.Mom.tsv", quote=F, row.names=F, sep="\t")
write.table(Aidm4Dad.anova, file="Aidm4.anovaSummary.Dad.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Aidm4Mom.emm <- emmeans(Aidm4Mom, pairwise ~ Targeting*dipMom*dipDad, lmer.df='satterthwaite',adjust='Holm')
Aidm4Mom.contrasts <- as.data.frame(Aidm4Mom.emm$contrasts)
Aidm4Mom.contrasts$model <- "Aidm4Mom"

Aidm4Dad.emm <- emmeans(Aidm4Dad, pairwise ~ Targeting*dipMom*dipDad, lmer.df='satterthwaite',adjust='Holm')
Aidm4Dad.contrasts <- as.data.frame(Aidm4Dad.emm$contrasts)
Aidm4Dad.contrasts$model <- "Aidm4Dad"

Aidm4.emm <- emmeans(Aidm4, pairwise ~ Targeting*dipMom*dipDad, lmer.df='satterthwaite',adjust='Holm')
Aidm4.contrasts <- as.data.frame(Aidm4.emm$contrasts)
Aidm4.contrasts$model <- "Aidm4full"

Aidm4.tbl <- rbind(Aidm4.contrasts,Aidm4Mom.contrasts,Aidm4Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,NA,"Category2",NA,NA),sep=" ") %>%
	filter(p.value<0.05)

write.table(Aidm4.tbl, file="Aidm4.SigContrasts.tsv", quote=F, row.names=F, sep="\t")











##### construct reciprocal df for model 5, Ah #####

deldataAh <- SSdataAh %>% 
	mutate(deltaMomDad = dipMom - dipDad) %>%
	select(c(Targeting,dadGene,momGene,plantID,tetMom,tetDad,deltaMomDad))


##### run model 5, Ah #####

Ahm5 <- lm(cbind(tetMom, tetDad) ~ Targeting*deltaMomDad, data = deldataAh)
summary(Ahm5)

Plotm5Ah <- emmip_ggplot(emmip(Ahm5, Targeting ~ deltaMomDad | subgenome, mult.name = "subgenome", at=list(deltaMomDad=seq(-10,10,by=2)),plotit=F))
ggsave("Ahm5.linearPrediction.jpg",Plotm5Ah,height=5,width=10)
goodPlotAh <- Plotm5Ah + scale_color_discrete(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"))
plot(goodPlotAh)

# also run the models independently to access results easier
# will give same results
Ahm5Mom <- lm(tetMom ~ Targeting*deltaMomDad, data = deldataAh)
Ahm5Dad <- lm(tetDad ~ Targeting*deltaMomDad, data = deldataAh)


# check residual normality
deldataAh$predictMom <- predict(Ahm5)[,1]
deldataAh$predictDad <- predict(Ahm5)[,2]

deldataAh$residMom <- deldataAh$tetMom - deldataAh$predictMom
deldataAh$residDad <- deldataAh$tetDad - deldataAh$predictDad

deldataAh$zresidMom <- (deldataAh$residMom - mean(deldataAh$residMom))/sd(deldataAh$residMom)
deldataAh$zresidDad <- (deldataAh$residDad - mean(deldataAh$residDad))/sd(deldataAh$residDad)

histMom5Ah <- hist(deldataAh$zresidMom, plot=F)
histDad5Ah <- hist(deldataAh$zresidDad, plot=F)

#jpeg(file="Ahm5.modelPlot.normality.jpeg")
plot(histMom5Ah, col=cMom)
plot(histDad5Ah, col=cDad, add=T)
#dev.off()


# heteroscedasticity
plot(deldataAh$zresidMom ~ as.factor(deldataAh$Targeting), col=cMom)
plot(deldataAh$zresidDad ~ as.factor(deldataAh$Targeting), col=cDad, add=T)

plot(deldataAh$zresidMom ~ deldataAh$tetDad, col=cMom)
plot(deldataAh$zresidDad ~ deldataAh$tetMom, col=cDad)



## make a table of fixed effects
Ahm5.fixedEffects <- as.data.frame(tidy(Ahm5) %>%
	mutate(term= str_replace(term,"_Non-interacting","NI")) %>% 
	mutate(term= str_replace(term,"Targeting","")) %>% 
	mutate(term= str_replace(term,"Dual-targeted","D")) %>% 
	mutate(term= str_replace(term,"Mitochondria-targeted","M")) %>% 
	mutate(term= str_replace(term,"Plastid-targeted","P")) %>% 
	mutate(term= str_replace(term,"_Interacting","I"))
)

write.table(Ahm5.fixedEffects, file="Ahm5.fixedEffects.tsv", quote=F, row.names=F, sep="\t") 

Ahm5.fixedEffects.sig <- Ahm5.fixedEffects %>%
	filter(p.value<0.05)



## run a type III anova and make table of summary
Ahm5.anova <- Anova(Ahm5, type=3)
Ahm5Mom.anova <- anova_summary(Anova(Ahm5Mom, type=3))
Ahm5Dad.anova <- anova_summary(Anova(Ahm5Dad, type=3))

write.table(Ahm5Mom.anova, file="Ahm5.anovaSummary.Mom.tsv", quote=F, row.names=F, sep="\t")
write.table(Ahm5Dad.anova, file="Ahm5.anovaSummary.Dad.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Ahm5Mom.emm <- emmeans(Ahm5Mom, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Ahm5Mom.contrasts <- as.data.frame(Ahm5Mom.emm$contrasts)
Ahm5Mom.contrasts$model <- "Ahm5Mom"

Ahm5Dad.emm <- emmeans(Ahm5Dad, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Ahm5Dad.contrasts <- as.data.frame(Ahm5Dad.emm$contrasts)
Ahm5Dad.contrasts$model <- "Ahm5Dad"

Ahm5.emm <- emmeans(Ahm5, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Ahm5.contrasts <- as.data.frame(Ahm5.emm$contrasts)
Ahm5.contrasts$model <- "Ahm5full"

Ahm5.tbl <- rbind(Ahm5.contrasts,Ahm5Mom.contrasts,Ahm5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") %>%
	filter(p.value<0.05)

write.table(Ahm5.tbl, file="Ahm5.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


Ahm5.all.tbl <- rbind(Ahm5.contrasts,Ahm5Mom.contrasts,Ahm5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") 

write.table(Ahm5.all.tbl, file="Ahm5.AllContrasts.tsv", quote=F, row.names=F, sep="\t")



##### construct reciprocal df for model 5, Aid #####

deldataAid <- SSdataAid %>% 
	mutate(deltaMomDad = dipMom - dipDad) %>%
	select(c(Targeting,dadGene,momGene,plantID,tetMom,tetDad,deltaMomDad))


##### run model 5, Aid #####

Aidm5 <- lm(cbind(tetMom, tetDad) ~ Targeting*deltaMomDad, data = deldataAid)
summary(Aidm5)

Plotm5Aid <- emmip_ggplot(emmip(Aidm5, Targeting ~ deltaMomDad | subgenome, mult.name = "subgenome", at=list(deltaMomDad=seq(-10,10,by=2)),plotit=F))
ggsave("Aidm5.linearPrediction.jpg",Plotm5Aid, height=5, width=10)
goodPlotAid <- Plotm5Aid + scale_color_discrete(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"))


# also run the models independently to access results easier
# will give same results
Aidm5Mom <- lm(tetMom ~ Targeting*deltaMomDad, data = deldataAid)
Aidm5Dad <- lm(tetDad ~ Targeting*deltaMomDad, data = deldataAid)


# check residual normality
deldataAid$predictMom <- predict(Aidm5)[,1]
deldataAid$predictDad <- predict(Aidm5)[,2]

deldataAid$residMom <- deldataAid$tetMom - deldataAid$predictMom
deldataAid$residDad <- deldataAid$tetDad - deldataAid$predictDad

deldataAid$zresidMom <- (deldataAid$residMom - mean(deldataAid$residMom))/sd(deldataAid$residMom)
deldataAid$zresidDad <- (deldataAid$residDad - mean(deldataAid$residDad))/sd(deldataAid$residDad)

histMom5Aid <- hist(deldataAid$zresidMom, plot=F)
histDad5Aid <- hist(deldataAid$zresidDad, plot=F)

#residuals look better here...
#jpeg(file="Aidm5.modelPlot.normality.jpeg")
plot(histMom5Aid, col=cMom)
plot(histDad5Aid, col=cDad, add=T)
#dev.off()


# heteroscedasticity, also looks better than m4
plot(deldataAid$zresidMom ~ as.factor(deldataAid$Targeting), col=cMom)
plot(deldataAid$zresidDad ~ as.factor(deldataAid$Targeting), col=cDad, add=T)

plot(deldataAid$zresidMom ~ deldataAid$tetDad, col=cMom)
plot(deldataAid$zresidDad ~ deldataAid$tetMom, col=cDad)



## make a table of fixed effects
Aidm5.fixedEffects <- as.data.frame(tidy(Aidm5) %>%
	mutate(term= str_replace(term,"_Non-interacting","NI")) %>% 
	mutate(term= str_replace(term,"Targeting","")) %>% 
	mutate(term= str_replace(term,"Dual-targeted","D")) %>% 
	mutate(term= str_replace(term,"Mitochondria-targeted","M")) %>% 
	mutate(term= str_replace(term,"Plastid-targeted","P")) %>% 
	mutate(term= str_replace(term,"_Interacting","I"))
)

write.table(Aidm5.fixedEffects, file="Aidm5.fixedEffects.tsv", quote=F, row.names=F, sep="\t") 

Aidm5.fixedEffects.sig <- Aidm5.fixedEffects %>%
	filter(p.value<0.05)



## run a type III anova and make table of summary
Aidm5.anova <- Anova(Aidm5, type=3)
Aidm5Mom.anova <- anova_summary(Anova(Aidm5Mom, type=3))
Aidm5Dad.anova <- anova_summary(Anova(Aidm5Dad, type=3))

write.table(Aidm5Mom.anova, file="Aidm5.anovaSummary.Mom.tsv", quote=F, row.names=F, sep="\t")
write.table(Aidm5Dad.anova, file="Aidm5.anovaSummary.Dad.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Aidm5Mom.emm <- emmeans(Aidm5Mom, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Aidm5Mom.contrasts <- as.data.frame(Aidm5Mom.emm$contrasts)
Aidm5Mom.contrasts$model <- "Aidm5Mom"

Aidm5Dad.emm <- emmeans(Aidm5Dad, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Aidm5Dad.contrasts <- as.data.frame(Aidm5Dad.emm$contrasts)
Aidm5Dad.contrasts$model <- "Aidm5Dad"

Aidm5.emm <- emmeans(Aidm5, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Aidm5.contrasts <- as.data.frame(Aidm5.emm$contrasts)
Aidm5.contrasts$model <- "Aidm5full"

Aidm5.tbl <- rbind(Aidm5.contrasts,Aidm5Mom.contrasts,Aidm5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") %>%
	filter(p.value<0.05)

write.table(Aidm5.tbl, file="Aidm5.SigContrasts.tsv", quote=F, row.names=F, sep="\t")



Aidm5.all.tbl <- rbind(Aidm5.contrasts,Aidm5Mom.contrasts,Aidm5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") %>%
	filter(p.value<0.05)

write.table(Aidm5.all.tbl, file="Aidm5.AllContrasts.tsv", quote=F, row.names=F, sep="\t")




































