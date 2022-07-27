setwd("Y:/corrinne/Cytonuclear/chenopodium/DEanalysis/linearModels")

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
rlogfull <- fread("chenopodium.lm.tbl", header=T, sep="\t") %>% 
	filter(str_detect(gene,"AUR")) %>%
	filter(!is.na(subgenome)) %>%
	filter(!is.na(Cq_1))

rlog <- rlogfull %>% 
	select(!contains(c("Cp","Cs")))
		
# read in and number the pairs then extract into two 2-column df
pairs <- fread("chenopodium.gene.pairs", col.names=c("category","dad","mom"), sep="\t") %>% 
	rowid_to_column("genePairID")
one <- pairs[, c(1, 3)]
two <- pairs[, c(1, 4)]
names(one) <- c("genePairID", "gene")
names(two) <- c("genePairID", "gene")


# join each df with the rlog df to assign genePairID to gene
data <- bind_rows(one,two) %>%
	right_join(.,rlog,by='gene') %>%
	select(!contains("CyMira")) %>% 
	pivot_longer(cols=starts_with("Cq"), names_to="plantID", values_to="rlog") %>%
	rename(Subgenome=subgenome) %>%
	rename(Targeting=category)


##### run model 1, Cq1 #####

## run mixed model
Cqm1 <- lmer(rlog ~ Subgenome * Targeting + (1 | genePairID) + (1 | plantID), data=data, REML=T)
summary(Cqm1)

## make a table of fixed effects
Cqm1.fixedEffects <- summary(Cqm1)$coefficients
write.table(Cqm1.fixedEffects, file="Cqm1.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Cqm1, type=2)
Cqm1.anova <- anova_summary(Anova(Cqm1, type=3))
write.table(Cqm1.anova, file="Cqm1.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## graph fixed effects with CI
Cqm1.emmp <- emmip(Cqm1, Subgenome ~ Targeting, plotit=F, CIs=T)
Cqm1.fixedPlot <- emmip_ggplot(Cqm1.emmp, style = "factor", dodge = 0.1, 
  facetlab = "label_context", CIarg = list(lwd = 2, alpha = 0.5),
  PIarg = list(lwd = 1.25, alpha = 0.33)) + theme(axis.text.x = element_text(angle = 45, hjust=1))
ggsave("Cqm1.fixedEffects.jpg", device="jpeg", units="in", height=8, width=12, plot=Cqm1.fixedPlot)

## perform contrasts and generate table
Cqm1.emm <- emmeans(Cqm1, pairwise ~ Subgenome*Targeting, lmer.df='satterthwaite',adjust='Holm')
Cqm1.contrasts <- as.data.frame(Cqm1.emm$contrasts)

Cqm1.tbl <- Cqm1.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Parent1","Category1","Parent2","Category2"),sep=" ") %>%
	filter( ((Parent1 == Parent2 & Category1 != Category2) | (Parent1 != Parent2 & Category1 == Category2)) & p.value<0.05)

write.table(Cqm1.tbl, file="Cqm1.SigContrasts.tsv", quote=F, row.names=F, sep="\t")



#####
# genePairID is not a good random predictor... too many df and not enough independent observations
# note the Inf degrees freedom
# make new model rlogRatio ~ Targeting + (1 | plantID)
# or other new model delta_rlog ~ Targeting + (1 | plantID)
#####

## make a df that has diploids combined
rlogmom <- rlogfull %>% filter(subgenome=="mom") %>% 
	select(contains(c("gene","Cp"))) %>% 
	setNames(gsub("Cp", "Dip", names(.))) %>%
	rename(Dip_3=Dip_7) %>%
	relocate(Dip_3, .after = Dip_2) %>%
	select(!Dip_1) 

rlogdad <- rlogfull %>% filter(subgenome=="dad") %>% 
	select(contains(c("gene","Cs"))) %>% 
	setNames(gsub("Cs", "Dip", names(.)))

rlogdip <- bind_rows(rlogmom,rlogdad) %>% 
	right_join(.,rlog,by='gene') %>% 
	select(!subgenome) 


##### count difference #####
# for uniformity, always mom - dad or mom/dad
# compared to momDip - dadDip or momDip/dadDip
homoeoPairs <- pairs[,3:4] %>% filter(dad %in% rlog$gene & mom %in% rlog$gene)

## make the delta
Drlog <- rlogdip[1:2,]
Drlog <- Drlog[-(1:2),]

for (i in 1:nrow(homoeoPairs)) {
    Drlog[nrow(Drlog)+1,] <- rbind(
        cbind(gene=paste0(homoeoPairs[i,2],"_",homoeoPairs[i,1]),
        rlogdip[rlogdip$gene==homoeoPairs[i,2], -c(1,11,12,13)] - 
        rlogdip[rlogdip$gene==homoeoPairs[i,1], -c(1,11,12,13)],
        rlogdip[rlogdip$gene==homoeoPairs[i,1], c(11:13)]))
}

Ddata <- Drlog %>% 
	select(!contains("CyMira")) %>%
	rename(Targeting=category) %>% 
	unite("Plant2", c(Dip_2,Cq_1), sep="X") %>% 
	unite("Plant3", c(Dip_3,Cq_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,Cq_4), sep="X") %>% 
	unite("Plant5", c(Dip_5,Cq_5), sep="X") %>%
	select(!Cq_6) %>%	 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipDelta","tetDelta"), sep="X") %>% 
	mutate(across(ends_with("Delta"), as.double))

fwrite(Ddata, file="Cq.delta.tsv", quote=F, row.names=F, sep="\t")


## make the ratio 
Rrlog <- rlogdip[1:2,]
Rrlog <- Rrlog[-(1:2),]

for (i in 1:nrow(homoeoPairs)) {
    Rrlog[nrow(Rrlog)+1,] <- rbind(
        cbind(gene=paste0(homoeoPairs[i,2],"_",homoeoPairs[i,1]),
        rlogdip[rlogdip$gene==homoeoPairs[i,2], -c(1,11,12,13)] / 
        rlogdip[rlogdip$gene==homoeoPairs[i,1], -c(1,11,12,13)],
        rlogdip[rlogdip$gene==homoeoPairs[i,1], c(11:13)]))
}


Rdata <- Rrlog %>% 
	select(!contains("CyMira")) %>%
	rename(Targeting=category) %>% 
	unite("Plant2", c(Dip_2,Cq_1), sep="X") %>% 
	unite("Plant3", c(Dip_3,Cq_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,Cq_4), sep="X") %>% 
	unite("Plant5", c(Dip_5,Cq_5), sep="X") %>% 
	select(!Cq_6) %>%	 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipRatio","tetRatio"), sep="X") %>%
	mutate(across(ends_with("Ratio"), as.double))
	
	
fwrite(Rdata, file="Cq.ratio.tsv", quote=F, row.names=F, sep="\t")


##### run model 2, Asue #####

### run mixed model, ratio ###
Cqm2R <- lmer(tetRatio ~ Targeting + (1 | plantID), data=Rdata, REML=T)
summary(Cqm2R)

# rerun, removing specified random effects
# this did not strongly influence first model
Cqm2R <- lm(tetRatio ~ Targeting, data=Rdata)
summary(Cqm2R)


## save the model plot
jpeg(file="Cqm2.Ratio.modelPlot.jpeg")
plot(Cqm2R)
dev.off()

## make a table of fixed effects
Cqm2R.fixedEffects <- summary(Cqm2R)$coefficients
write.table(Cqm2R.fixedEffects, file="Cqm2.Ratio.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Cqm2R, type=2)
Cqm2R.anova <- anova_summary(Anova(Cqm2R, type=3))
write.table(Cqm2R.anova, file="Cqm2.Ratio.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## perform contrasts and generate table
Cqm2R.emm <- emmeans(Cqm2R, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Cqm2R.contrasts <- as.data.frame(Cqm2R.emm$contrasts)

Cqm2R.tbl <- Cqm2R.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") # %>%
#	filter( ((Parent1 == Parent2 & Category1 != Category2) | (Parent1 != Parent2 & Category1 == Category2)) & p.value<0.05)

write.table(Cqm2R.tbl, file="Cqm2.Ratio.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


### run mixed model, delta ###
Cqm2D <- lmer(tetDelta ~ Targeting + (1 | plantID), data=Ddata, REML=T)
summary(Cqm2D)

Cqm2D <- lm(tetDelta ~ Targeting, data=Ddata)
summary(Cqm2D)

## save the model plot
jpeg(file="Cqm2.Delta.modelPlot.jpeg")
plot(Cqm2D)
dev.off()

## make a table of fixed effects
Cqm2D.fixedEffects <- summary(Cqm2D)$coefficients
write.table(Cqm2D.fixedEffects, file="Cqm2.Delta.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Cqm2D, type=2)
Cqm2D.anova <- anova_summary(Anova(Cqm2D, type=3))
write.table(Cqm2D.anova, file="Cqm2.Delta.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## perform contrasts and generate table
Cqm2D.emm <- emmeans(Cqm2D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Cqm2D.contrasts <- as.data.frame(Cqm2D.emm$contrasts)

Cqm2D.tbl <- Cqm2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") # %>%
#	filter( ((Parent1 == Parent2 & Category1 != Category2) | (Parent1 != Parent2 & Category1 == Category2)) & p.value<0.05)

write.table(Cqm2D.tbl, file="Cqm2.Delta.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


##### run model 3, Cq #####

### run mixed model 3, ratio ###
Cqm3R <- lmer(tetRatio ~ Targeting*dipRatio + (1 | plantID), data=Rdata, REML=T)
summary(Cqm3R)

Cqm3R <- lm(tetRatio ~ Targeting*dipRatio, data=Rdata)
summary(Cqm3R)

Cqm3Rr <- lm(tetRatio ~ dipRatio*Targeting, data=Rdata)
summary(Cqm3Rr)


## save the model plot
jpeg(file="Cqm3.Ratio.modelPlot.jpeg")
plot(Cqm3R)
dev.off()

## make a table of fixed effects
Cqm3R.fixedEffects <- summary(Cqm3R)$coefficients
write.table(Cqm3R.fixedEffects, file="Cqm3.Ratio.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Cqm3R, type=2)
Cqm3R.anova <- anova_summary(Anova(Cqm3R, type=3))
write.table(Cqm3R.anova, file="Cqm3.Ratio.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## graph fixed effects with CI; does not apply bc no categorical interacting terms
# Cqm3R.emmp <- emmip(Cqm3R, dipRatio ~ Targeting, plotit=F, CIs=T)
# Cqm3R.fixedPlot <- emmip_ggplot(Cqm3R.emmp, style = "factor", dodge = 0.1, 
#   facetlab = "label_context", CIarg = list(lwd = 2, alpha = 0.5),
#   PIarg = list(lwd = 1.25, alpha = 0.33)) + theme(axis.text.x = element_text(angle = 45, hjust=1))
# ggsave("Cqm3.Ratio.fixedEffects.jpg", device="jpeg", units="in", height=8, width=12, plot=Cqm3R.fixedPlot)

## perform contrasts and generate table
Cqm3R.emm <- emmeans(Cqm3R, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Cqm3R.contrasts <- as.data.frame(Cqm3R.emm$contrasts)

Cqm3R.tbl <- Cqm3R.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") # %>%
#	filter( ((Parent1 == Parent2 & Category1 != Category2) | (Parent1 != Parent2 & Category1 == Category2)) & p.value<0.05)

write.table(Cqm3R.tbl, file="Cqm3.Ratio.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


### run mixed model 3, delta ###
Cqm3D <- lmer(tetDelta ~ Targeting*dipDelta + (1 | plantID), data=Ddata, REML=T)
summary(Cqm3D)

## save the model plot
jpeg(file="Cqm3.Delta.modelPlot.jpeg")
plot(Cqm3D)
dev.off()

## make a table of fixed effects
Cqm3D.fixedEffects <- summary(Cqm3D)$coefficients
write.table(Cqm3D.fixedEffects, file="Cqm3.Delta.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Cqm3D, type=2)
Cqm3D.anova <- anova_summary(Anova(Cqm3D, type=3))
write.table(Cqm3D.anova, file="Cqm3.Delta.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## graph fixed effects with CI; does not apply bc no categorical interacting terms
# Cqm3D.emmp <- emmip(Cqm3D, Subgenome ~ Targeting, plotit=F, CIs=T)
# Cqm3D.fixedPlot <- emmip_ggplot(Cqm3D.emmp, style = "factor", dodge = 0.1, 
#   facetlab = "label_context", CIarg = list(lwd = 2, alpha = 0.5),
#   PIarg = list(lwd = 1.25, alpha = 0.33)) + theme(axis.text.x = element_text(angle = 45, hjust=1))
# ggsave("Cqm3.Delta.fixedEffects.jpg", device="jpeg", units="in", height=8, width=12, plot=Cqm3D.fixedPlot)

## perform contrasts and generate table
Cqm3D.emm <- emmeans(Cqm3D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Cqm3D.contrasts <- as.data.frame(Cqm3D.emm$contrasts)

Cqm3D.tbl <- Cqm3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") # %>%
#	filter( ((Parent1 == Parent2 & Category1 != Category2) | (Parent1 != Parent2 & Category1 == Category2)) & p.value<0.05)

write.table(Cqm3D.tbl, file="Cqm3.Delta.SigContrasts.tsv", quote=F, row.names=F, sep="\t")




#####
# it weirds me out that delta and ratio give different results
# combining mom & dad doesn't speak to the subgenomes directly
# let's model tetMom, tetDad separately together for Model 4
#####


SSrlogmom <- rlogfull %>% filter(subgenome=="mom") %>% 
	select(contains(c("gene","Cp","Cq"))) %>% 
	setNames(gsub("Cp", "dipMom", names(.))) %>% 
	setNames(gsub("Cq", "tetMom", names(.))) %>% 
	rename(momGene = gene)

SSrlogdad <- rlogfull %>% filter(subgenome=="dad") %>% 
	select(contains(c("gene","Cs","Cq"))) %>% 
	setNames(gsub("Cs", "dipDad", names(.))) %>% 
	setNames(gsub("Cq", "tetDad", names(.))) %>% 
	rename(dadGene = gene)


SSrlog <- pairs %>% select(!genePairID) %>%
	rename(Targeting = category) %>%
	rename(dadGene = dad) %>%
	rename(momGene = mom)  %>% 
	left_join(.,SSrlogmom,by='momGene') %>% 
	left_join(.,SSrlogdad,by='dadGene') %>%
	filter(!(is.na(dipMom_1) | is.na(dipDad_2)))



SSdata <- SSrlog %>% 
	unite("PlantA", c(dipMom_1, dipDad_2, tetMom_1, tetDad_1), sep="X") %>% 
	unite("PlantB", c(dipMom_2, dipDad_3, tetMom_3, tetDad_3), sep="X") %>% 
	unite("PlantC", c(dipMom_4, dipDad_4, tetMom_4, tetDad_4), sep="X") %>% 
	unite("PlantD", c(dipMom_5, dipDad_5, tetMom_5, tetDad_5), sep="X") %>%
	select(!c(dipMom_7,tetMom_6,tetDad_6)) %>%	 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipMom","dipDad","tetMom","tetDad"), sep="X") %>% 
	mutate(across(ends_with("Mom"), as.double)) %>% 
	mutate(across(ends_with("Dad"), as.double)) %>%
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-")) # put NOT first alphabetically

fwrite(SSdata, file="Cq.SS.tsv", quote=F, row.names=F, sep="\t")

# check data
plot(SSdata$tetDad ~ SSdata$tetMom)

SSplot <- SSdata %>% 
	select(Targeting, tetMom, tetDad) %>%
	mutate(Targeting = str_replace(Targeting,"-targeted_Non-interacting","")) %>%
	mutate(Targeting = str_replace(Targeting,"-targeted_Interacting","")) %>%
	mutate(Targeting = str_replace(Targeting,"-organelle-targeted","")) %>% 
	arrange(Targeting) %>%
	ggplot(aes(x=tetMom, y=tetDad, color=Targeting)) + 
	   geom_point(size=0.9) + 
	   scale_color_manual(values=c('lightgrey','red3','blue4', 'green4')) +
	   geom_abline(intercept = 0, slope = 1) +
	   ggtitle("Chenopodium quinoa") + 
	   theme(plot.title = element_text(face = "italic")) + 
	   labs(x="Maternal Homoeolog Expression (rlog)", y="Paternal Homoeolog Expression (rlog)")

ggsave("CqMomvsDad.jpg",SSplot)



##### run model 4, Cq #####

Cqm4 <- lm(cbind(tetMom, tetDad) ~ Targeting*dipMom*dipDad, data = SSdata)
summary(Cqm4)

momPlot <- emmip_ggplot(emmip(Cqm4, Targeting ~ dipMom | subgenome, mult.name = "subgenome", at=list(dipMom=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
dadPlot <- emmip_ggplot(emmip(Cqm4, Targeting ~ dipDad | subgenome, mult.name = "subgenome", at=list(dipDad=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))

momdadPlot <- grid.arrange(momPlot,dadPlot,ncol=1)
ggsave("Cqm4.linearPrediction.jpg",momdadPlot)


# also run the models independently to access results easier
# will give same results
Cqm4Mom <- lm(tetMom ~ Targeting*dipMom*dipDad, data = SSdata)
Cqm4Dad <- lm(tetDad ~ Targeting*dipMom*dipDad, data = SSdata)


# check residual normality
SSdata$predictMom <- predict(Cqm4)[,1]
SSdata$predictDad <- predict(Cqm4)[,2]

SSdata$residMom <- SSdata$tetMom - SSdata$predictMom
SSdata$residDad <- SSdata$tetDad - SSdata$predictDad

SSdata$zresidMom <- (SSdata$residMom - mean(SSdata$residMom))/sd(SSdata$residMom)
SSdata$zresidDad <- (SSdata$residDad - mean(SSdata$residDad))/sd(SSdata$residDad)

histMom <- hist(SSdata$zresidMom, plot=F)
histDad <- hist(SSdata$zresidDad, plot=F)

cDad <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
cMom <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")

#jpeg(file="Cqm4.modelPlot.normality.jpeg")
plot(histMom, col=cMom)
plot(histDad, col=cDad, add=T)
#dev.off()


# heteroscedasticity
plot(SSdata$zresidMom ~ as.factor(SSdata$Targeting), col=cMom)
plot(SSdata$zresidDad ~ as.factor(SSdata$Targeting), col=cDad, add=T)

plot(SSdata$zresidMom ~ SSdata$tetDad, col=cMom)
plot(SSdata$zresidDad ~ SSdata$tetMom, col=cDad)



## make a table of fixed effects
Cqm4.fixedEffects <- as.data.frame(tidy(Cqm4) %>%
	mutate(term= str_replace(term,"_Non-interacting","NI")) %>% 
	mutate(term= str_replace(term,"Targeting","")) %>% 
	mutate(term= str_replace(term,"Dual-targeted","D")) %>% 
	mutate(term= str_replace(term,"Mitochondria-targeted","M")) %>% 
	mutate(term= str_replace(term,"Plastid-targeted","P")) %>% 
	mutate(term= str_replace(term,"_Interacting","I"))
)

write.table(Cqm4.fixedEffects, file="Cqm4.fixedEffects.tsv", quote=F, row.names=F, sep="\t") 

Cqm4.fixedEffects.sig <- Cqm4.fixedEffects %>%
	filter(p.value<0.05)



## run a type III anova and make table of summary
Cqm4.anova <- Anova(Cqm4, type=3)
Cqm4Mom.anova <- anova_summary(Anova(Cqm4Mom, type=3))
Cqm4Dad.anova <- anova_summary(Anova(Cqm4Dad, type=3))

write.table(Cqm4Mom.anova, file="Cqm4.anovaSummary.Mom.tsv", quote=F, row.names=F, sep="\t")
write.table(Cqm4Dad.anova, file="Cqm4.anovaSummary.Dad.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Cqm4Mom.emm <- emmeans(Cqm4Mom, pairwise ~ Targeting*dipMom*dipDad, lmer.df='satterthwaite',adjust='Holm')
Cqm4Mom.contrasts <- as.data.frame(Cqm4Mom.emm$contrasts)
Cqm4Mom.contrasts$model <- "Cqm4Mom"

Cqm4Dad.emm <- emmeans(Cqm4Dad, pairwise ~ Targeting*dipMom*dipDad, lmer.df='satterthwaite',adjust='Holm')
Cqm4Dad.contrasts <- as.data.frame(Cqm4Dad.emm$contrasts)
Cqm4Dad.contrasts$model <- "Cqm4Dad"

Cqm4.emm <- emmeans(Cqm4, pairwise ~ Targeting*dipMom*dipDad, lmer.df='satterthwaite',adjust='Holm')
Cqm4.contrasts <- as.data.frame(Cqm4.emm$contrasts)
Cqm4.contrasts$model <- "Cqm4full"

Cqm4.tbl <- rbind(Cqm4.contrasts,Cqm4Mom.contrasts,Cqm4Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,NA,"Category2",NA,NA),sep=" ") %>%
	filter(p.value<0.05)

write.table(Cqm4.tbl, file="Cqm4.SigContrasts.tsv", quote=F, row.names=F, sep="\t")









##### construct reciprocal df for model 5, Cq #####

deldata <- SSdata %>% 
	mutate(deltaMomDad = dipMom - dipDad) %>%
	mutate(deltaDadMom = dipDad - dipMom) %>%
	select(c(Targeting,dadGene,momGene,plantID,tetMom,tetDad,deltaMomDad,deltaDadMom))


##### run model 5, Cq #####

Cqm5 <- lm(cbind(tetMom, tetDad) ~ Targeting*deltaMomDad, data = deldata)
summary(Cqm5)

Plotm5 <- emmip_ggplot(emmip(Cqm5, Targeting ~ deltaMomDad | subgenome, mult.name = "subgenome", at=list(deltaMomDad=seq(-10,10,by=2)),plotit=F))
ggsave("Cqm5.linearPrediction.jpg",Plotm5)
goodPlot <- Plotm5 + scale_color_discrete(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"))


# also run the models independently to access results easier
# will give same results
Cqm5Mom <- lm(tetMom ~ Targeting*deltaMomDad, data = deldata)
Cqm5Dad <- lm(tetDad ~ Targeting*deltaMomDad, data = deldata)


# check residual normality
deldata$predictMom <- predict(Cqm5)[,1]
deldata$predictDad <- predict(Cqm5)[,2]

deldata$residMom <- deldata$tetMom - deldata$predictMom
deldata$residDad <- deldata$tetDad - deldata$predictDad

deldata$zresidMom <- (deldata$residMom - mean(deldata$residMom))/sd(deldata$residMom)
deldata$zresidDad <- (deldata$residDad - mean(deldata$residDad))/sd(deldata$residDad)

histMom5 <- hist(deldata$zresidMom, plot=F)
histDad5 <- hist(deldata$zresidDad, plot=F)

#jpeg(file="Cqm5.modelPlot.normality.jpeg")
plot(histMom5, col=cMom)
plot(histDad5, col=cDad, add=T)
#dev.off()


# heteroscedasticity
plot(deldata$zresidMom ~ as.factor(deldata$Targeting), col=cMom)
plot(deldata$zresidDad ~ as.factor(deldata$Targeting), col=cDad, add=T)

plot(deldata$zresidMom ~ deldata$tetDad, col=cMom)
plot(deldata$zresidDad ~ deldata$tetMom, col=cDad)



## make a table of fixed effects
Cqm5.fixedEffects <- as.data.frame(tidy(Cqm5) %>%
	mutate(term= str_replace(term,"_Non-interacting","NI")) %>% 
	mutate(term= str_replace(term,"Targeting","")) %>% 
	mutate(term= str_replace(term,"Dual-targeted","D")) %>% 
	mutate(term= str_replace(term,"Mitochondria-targeted","M")) %>% 
	mutate(term= str_replace(term,"Plastid-targeted","P")) %>% 
	mutate(term= str_replace(term,"_Interacting","I"))
)

write.table(Cqm5.fixedEffects, file="Cqm5.fixedEffects.tsv", quote=F, row.names=F, sep="\t") 

Cqm5.fixedEffects.sig <- Cqm5.fixedEffects %>%
	filter(p.value<0.05)



## run a type III anova and make table of summary
Cqm5.anova <- Anova(Cqm5, type=3)
Cqm5Mom.anova <- anova_summary(Anova(Cqm5Mom, type=3))
Cqm5Dad.anova <- anova_summary(Anova(Cqm5Dad, type=3))

write.table(Cqm5Mom.anova, file="Cqm5.anovaSummary.Mom.tsv", quote=F, row.names=F, sep="\t")
write.table(Cqm5Dad.anova, file="Cqm5.anovaSummary.Dad.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Cqm5Mom.emm <- emmeans(Cqm5Mom, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Cqm5Mom.contrasts <- as.data.frame(Cqm5Mom.emm$contrasts)
Cqm5Mom.contrasts$model <- "Cqm5Mom"

Cqm5Dad.emm <- emmeans(Cqm5Dad, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Cqm5Dad.contrasts <- as.data.frame(Cqm5Dad.emm$contrasts)
Cqm5Dad.contrasts$model <- "Cqm5Dad"

Cqm5.emm <- emmeans(Cqm5, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Cqm5.contrasts <- as.data.frame(Cqm5.emm$contrasts)
Cqm5.contrasts$model <- "Cqm5full"

Cqm5.tbl <- rbind(Cqm5.contrasts,Cqm5Mom.contrasts,Cqm5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") %>%
	filter(p.value<0.05)

write.table(Cqm5.tbl, file="Cqm5.SigContrasts.tsv", quote=F, row.names=F, sep="\t")




