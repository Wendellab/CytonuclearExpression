setwd("w:/corrinne/Cytonuclear/cotton/DEanalysis/linearModels")

library(tidyverse)
library(lme4)
library(car)
library(sjPlot)
library(rstatix)
library(emmeans)
library(data.table)
library(gridExtra)
library(lmtest)

##### global options #####
emm_options(lmerTest.limit = 100000)
setDTthreads(10)
options(datatable.fread.datatable=FALSE)

# read and filter the table
rlogfull <-fread("cotton.lm.tbl", header=T, sep="\t") %>% 
	filter(str_detect(gene,"Gohir")) %>% 
	filter(!is.na(subgenome)) %>%
	filter(!is.na(AD1_1)) %>% 
	rename(Subgenome=subgenome) %>%
	rename(Targeting=category)%>%
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-"))

rlog <- rlogfull %>% 
	select(!contains(c("A2","D5")))

# read in and number the pairs then extract into two 2-column df
pairs <- fread("cotton.gene.pairs", col.names=c("category","dad","mom"), sep="\t") %>% rowid_to_column("genePairID")


## make a df that has diploids combined
rlogmom <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","A2"))) %>% 
	setNames(gsub("A2", "Dip", names(.)))

rlogdad <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","D5"))) %>% 
	setNames(gsub("D5", "Dip", names(.)))

rlogdip <- bind_rows(rlogmom,rlogdad) %>% 
	right_join(.,rlog,by='gene') %>%
	select(!Subgenome) 




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
        rlogdip[rlogdip$gene==homoeoPairs[i,2], -c(1,17)] - 
        rlogdip[rlogdip$gene==homoeoPairs[i,1], -c(1,17)],
        rlogdip[rlogdip$gene==homoeoPairs[i,1], c(17)]))
}

DdataAD1 <- Drlog %>% 
	select(!contains("AD2")) %>%
	unite("Plant1", c(Dip_1,AD1_1), sep="X") %>% 
	unite("Plant2", c(Dip_2,AD1_2), sep="X") %>% 
	unite("Plant3", c(Dip_3,AD1_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,AD1_4), sep="X") %>% 
	unite("Plant5", c(Dip_5,AD1_5), sep="X") %>% 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipDelta","tetDelta"), sep="X") %>%
	mutate(across(ends_with("Delta"), as.double))



DdataAD2 <- Drlog %>% 
	select(!contains("AD1")) %>%
	unite("Plant1", c(Dip_1,AD2_1), sep="X") %>% 
	unite("Plant2", c(Dip_2,AD2_2), sep="X") %>% 
	unite("Plant3", c(Dip_3,AD2_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,AD2_4), sep="X") %>% 
	unite("Plant5", c(Dip_5,AD2_5), sep="X") %>% 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipDelta","tetDelta"), sep="X") %>%
	mutate(across(ends_with("Delta"), as.double))

fwrite(DdataAD1, file="AD1.delta.tbl", quote=F, row.names=F, sep="\t")
fwrite(DdataAD2, file="AD2.delta.tbl", quote=F, row.names=F, sep="\t")




##### run model 2, AD1 #####

### run mixed model, delta ###
AD1m2D <- lmer(tetDelta ~ Targeting + (1 | plantID), data=DdataAD1, REML=T)
summary(AD1m2D)

AD1m2D <- lm(tetDelta ~ Targeting, data=DdataAD1)
summary(AD1m2D)

## save the model plot
jpeg(file="AD1m2.Delta.modelPlot.jpeg")
plot(AD1m2D)
dev.off()

## make a table of fixed effects
AD1m2D.fixedEffects <- summary(AD1m2D)$coefficients
write.table(AD1m2D.fixedEffects, file="AD1m2.Delta.fixedEffects.tbl", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(AD1m2D, type=2)
AD1m2D.anova <- anova_summary(Anova(AD1m2D, type=3))
write.table(AD1m2D.anova, file="AD1m2.Delta.anovaSummary.tbl", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
AD1m2D.emm <- emmeans(AD1m2D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
AD1m2D.contrasts <- as.data.frame(AD1m2D.emm$contrasts)

AD1m2D.tbl <- AD1m2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ")  %>%
	filter(p.value<0.05)

write.table(AD1m2D.tbl, file="AD1m2.Delta.SigContrasts.tbl", quote=F, row.names=F, sep="\t")


AD1m2D.all.tbl <- AD1m2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ")  
write.table(AD1m2D.all.tbl, file="AD1m2.Delta.AllContrasts.tbl", quote=F, row.names=F, sep="\t")



##### run model 2, AD2 #####
### run mixed model, delta ###
AD2m2D <- lmer(tetDelta ~ Targeting + (1 | plantID), data=DdataAD2, REML=T)
summary(AD2m2D)

AD2m2D <- lm(tetDelta ~ Targeting, data=DdataAD2)
summary(AD2m2D)

## save the model plot
jpeg(file="AD2m2.Delta.modelPlot.jpeg")
plot(AD2m2D)
dev.off()

## make a table of fixed effects
AD2m2D.fixedEffects <- summary(AD2m2D)$coefficients
write.table(AD2m2D.fixedEffects, file="AD2m2.Delta.fixedEffects.tbl", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(AD2m2D, type=2)
AD2m2D.anova <- anova_summary(Anova(AD2m2D, type=3))
write.table(AD2m2D.anova, file="AD2m2.Delta.anovaSummary.tbl", quote=F, row.names=F, sep="\t")

## perform contrasts and generate table
AD2m2D.emm <- emmeans(AD2m2D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
AD2m2D.contrasts <- as.data.frame(AD2m2D.emm$contrasts)

AD2m2D.tbl <- AD2m2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") %>%
	filter(p.value<0.05)

write.table(AD2m2D.tbl, file="AD2m2.Delta.SigContrasts.tbl", quote=F, row.names=F, sep="\t")


AD2m2D.all.tbl <- AD2m2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") 
write.table(AD2m2D.all.tbl, file="AD2m2.Delta.AllContrasts.tbl", quote=F, row.names=F, sep="\t")

##### run model 3, AD1 #####

### run mixed model 3, delta ###
AD1m3D <- lmer(tetDelta ~ Targeting*dipDelta + (1 | plantID), data=DdataAD1, REML=T)
summary(AD1m3D)

AD1m3D <- lm(tetDelta ~ Targeting*dipDelta, data=DdataAD1)
summary(AD1m3D)

## save the model plot
jpeg(file="AD1m3.Delta.modelPlot.jpeg")
plot(AD1m3D)
dev.off()

## make a table of fixed effects
AD1m3D.fixedEffects <- summary(AD1m3D)$coefficients
write.table(AD1m3D.fixedEffects, file="AD1m3.Delta.fixedEffects.tbl", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(AD1m3D, type=2)
AD1m3D.anova <- anova_summary(Anova(AD1m3D, type=3))
write.table(AD1m3D.anova, file="AD1m3.Delta.anovaSummary.tbl", quote=F, row.names=F, sep="\t")

## graph fixed effects 
emmip(AD1m3D, Subgenome ~ Targeting, cov.reduce=F)


## perform contrasts and generate table
AD1m3D.emm <- emmeans(AD1m3D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
AD1m3D.contrasts <- as.data.frame(AD1m3D.emm$contrasts)

AD1m3D.tbl <- AD1m3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") %>%
	filter( p.value<0.05)

write.table(AD1m3D.tbl, file="AD1m3.Delta.SigContrasts.tbl", quote=F, row.names=F, sep="\t")



AD1m3D.all.tbl <- AD1m3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") 
write.table(AD1m3D.all.tbl, file="AD1m3.Delta.AllContrasts.tbl", quote=F, row.names=F, sep="\t")


##### run model 3, AD2 #####

### run mixed model 3, delta ###
AD2m3D <- lmer(tetDelta ~ Targeting*dipDelta + (1 | plantID), data=DdataAD2, REML=T)
summary(AD2m3D)

AD2m3D <- lm(tetDelta ~ Targeting*dipDelta, data=DdataAD2)
summary(AD2m3D)


## save the model plot
jpeg(file="AD2m3.Delta.modelPlot.jpeg")
plot(AD2m3D)
dev.off()

## make a table of fixed effects
AD2m3D.fixedEffects <- summary(AD2m3D)$coefficients
write.table(AD2m3D.fixedEffects, file="AD2m3.Delta.fixedEffects.tbl", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(AD2m3D, type=2)
AD2m3D.anova <- anova_summary(Anova(AD2m3D, type=3))
write.table(AD2m3D.anova, file="AD2m3.Delta.anovaSummary.tbl", quote=F, row.names=F, sep="\t")

## graph fixed effects 
emmip(AD2m3D, Subgenome ~ Targeting, cov.reduce=F)


## perform contrasts and generate table
AD2m3D.emm <- emmeans(AD2m3D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
AD2m3D.contrasts <- as.data.frame(AD2m3D.emm$contrasts)

AD2m3D.tbl <- AD2m3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ")  %>%
	filter( p.value<0.05)

write.table(AD2m3D.tbl, file="AD2m3.Delta.SigContrasts.tbl", quote=F, row.names=F, sep="\t")



AD2m3D.all.tbl <- AD2m3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ")  
write.table(AD2m3D.all.tbl, file="AD2m3.Delta.AllContrasts.tbl", quote=F, row.names=F, sep="\t")

lrtest(AD1m3D,AD1m2D)
lrtest(AD2m3D,AD2m2D)

#####
# let's model tetMom, tetDad separately together for Model 4
#####


SSrlogmom <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","A2","AD"))) %>% 
	setNames(gsub("A2", "dipMom", names(.))) %>% 
	setNames(gsub("AD", "tetMom", names(.))) %>% 
	rename(momGene = gene)

SSrlogdad <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","D5","AD"))) %>% 
	setNames(gsub("D5", "dipDad", names(.))) %>% 
	setNames(gsub("AD", "tetDad", names(.))) %>% 
	rename(dadGene = gene)


SSrlog <- pairs %>% select(!genePairID) %>%
	rename(Targeting = category) %>%
	rename(dadGene = dad) %>%
	rename(momGene = mom)  %>% 
	left_join(.,SSrlogmom,by='momGene') %>% 
	left_join(.,SSrlogdad,by='dadGene') %>%
	filter(!(is.na(dipMom_1) | is.na(dipDad_1)))%>%
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-"))




SSdataAD1 <- SSrlog %>% 
	select(!contains("2_")) %>% 
	setNames(gsub("1_", "_", names(.))) %>% 
	unite("PlantA", c(dipMom_1, dipDad_1, tetMom_1, tetDad_1), sep="X") %>% 
	unite("PlantB", c(dipMom_2, dipDad_2, tetMom_2, tetDad_2), sep="X") %>% 
	unite("PlantC", c(dipMom_3, dipDad_3, tetMom_3, tetDad_3), sep="X") %>% 
	unite("PlantD", c(dipMom_4, dipDad_4, tetMom_4, tetDad_4), sep="X") %>% 
	unite("PlantE", c(dipMom_5, dipDad_5, tetMom_5, tetDad_5), sep="X") %>%
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipMom","dipDad","tetMom","tetDad"), sep="X") %>% 
	mutate(across(ends_with("Mom"), as.double)) %>% 
	mutate(across(ends_with("Dad"), as.double))


SSdataAD2 <- SSrlog %>% 
	select(!contains("1_")) %>% 
	setNames(gsub("2_", "_", names(.))) %>% 
	unite("PlantA", c(dipMom_1, dipDad_1, tetMom_1, tetDad_1), sep="X") %>% 
	unite("PlantB", c(dipMom_2, dipDad_2, tetMom_2, tetDad_2), sep="X") %>% 
	unite("PlantC", c(dipMom_3, dipDad_3, tetMom_3, tetDad_3), sep="X") %>% 
	unite("PlantD", c(dipMom_4, dipDad_4, tetMom_4, tetDad_4), sep="X") %>% 
	unite("PlantE", c(dipMom_5, dipDad_5, tetMom_5, tetDad_5), sep="X") %>%
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipMom","dipDad","tetMom","tetDad"), sep="X") %>% 
	mutate(across(ends_with("Mom"), as.double)) %>% 
	mutate(across(ends_with("Dad"), as.double))


fwrite(SSdataAD1, file="AD1.SS.tsv", quote=F, row.names=F, sep="\t")
fwrite(SSdataAD2, file="AD2.SS.tsv", quote=F, row.names=F, sep="\t")

# check data
plot(SSdataAD1$tetDad ~ SSdataAD1$tetMom)
plot(SSdataAD2$tetDad ~ SSdataAD2$tetMom)

SSplotAD1 <- SSdataAD1 %>% 
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

ggsave("GhMomvsDad.jpg",SSplotAD1)


SSplotAD2 <- SSdataAD1 %>% 
	select(Targeting, tetMom, tetDad) %>%
	mutate(Targeting = str_replace(Targeting,"-targeted_Non-interacting","")) %>%
	mutate(Targeting = str_replace(Targeting,"-targeted_Interacting","")) %>%
	mutate(Targeting = str_replace(Targeting,"-organelle-targeted","")) %>% 
	arrange(Targeting) %>%
	ggplot(aes(x=tetMom, y=tetDad, color=Targeting)) + 
	   geom_point(size=0.9) + 
	   scale_color_manual(values=c('lightgrey','red3','blue4', 'green4')) +
	   geom_abline(intercept = 0, slope = 1) +
	   ggtitle("Gossypium barbadense") + 
	   theme(plot.title = element_text(face = "italic")) + 
	   labs(x="Maternal Homoeolog Expression (rlog)", y="Paternal Homoeolog Expression (rlog)")

ggsave("GbMomvsDad.jpg",SSplotAD2)





##### construct reciprocal df for model 5, AD1 #####

deldataAD1 <- SSdataAD1 %>% 
	mutate(deltaMomDad = dipMom - dipDad) %>%
	mutate(deltaDadMom = dipDad - dipMom) %>%
	select(c(Targeting,dadGene,momGene,plantID,tetMom,tetDad,deltaMomDad,deltaDadMom))


##### run model 5, AD1 #####

AD1m5 <- lm(cbind(tetMom, tetDad) ~ Targeting*deltaMomDad, data = deldataAD1)
summary(AD1m5)

Plotm5AD1 <- emmip_ggplot(emmip(AD1m5, Targeting ~ deltaMomDad | subgenome, mult.name = "subgenome", at=list(deltaMomDad=seq(-10,10,by=2)),plotit=F))
goodPlotAD1 <- Plotm5AD1 + scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) 
ggsave("AD1m5.linearPrediction.jpg",goodPlotAD1,units="in", width=10, height=5, dpi=400)



# also run the models independently to access results easier
# will give same results
AD1m5Mom <- lm(tetMom ~ Targeting*deltaMomDad, data = deldataAD1)
AD1m5Dad <- lm(tetDad ~ Targeting*deltaMomDad, data = deldataAD1)


# check residual normality
deldataAD1$predictMom <- predict(AD1m5)[,1]
deldataAD1$predictDad <- predict(AD1m5)[,2]

deldataAD1$residMom <- deldataAD1$tetMom - deldataAD1$predictMom
deldataAD1$residDad <- deldataAD1$tetDad - deldataAD1$predictDad

deldataAD1$zresidMom <- (deldataAD1$residMom - mean(deldataAD1$residMom))/sd(deldataAD1$residMom)
deldataAD1$zresidDad <- (deldataAD1$residDad - mean(deldataAD1$residDad))/sd(deldataAD1$residDad)

histMom5AD1 <- hist(deldataAD1$zresidMom, plot=F)
histDad5AD1 <- hist(deldataAD1$zresidDad, plot=F)

## make a table of fixed effects
AD1m5.fixedEffects <- as.data.frame(tidy(AD1m5) %>%
	mutate(term= str_replace(term,"_Non-interacting","NI")) %>% 
	mutate(term= str_replace(term,"Targeting","")) %>% 
	mutate(term= str_replace(term,"Dual-targeted","D")) %>% 
	mutate(term= str_replace(term,"Mitochondria-targeted","M")) %>% 
	mutate(term= str_replace(term,"Plastid-targeted","P")) %>% 
	mutate(term= str_replace(term,"_Interacting","I"))
)

write.table(AD1m5.fixedEffects, file="AD1m5.fixedEffects.tsv", quote=F, row.names=F, sep="\t") 

AD1m5.fixedEffects.sig <- AD1m5.fixedEffects %>%
	filter(p.value<0.05)



## run a type III anova and make table of summary
AD1m5.anova <- Anova(AD1m5, type=3)
AD1m5Mom.anova <- anova_summary(Anova(AD1m5Mom, type=3))
AD1m5Dad.anova <- anova_summary(Anova(AD1m5Dad, type=3))

write.table(AD1m5Mom.anova, file="AD1m5.anovaSummary.Mom.tsv", quote=F, row.names=F, sep="\t")
write.table(AD1m5Dad.anova, file="AD1m5.anovaSummary.Dad.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
AD1m5Mom.emm <- emmeans(AD1m5Mom, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
AD1m5Mom.contrasts <- as.data.frame(AD1m5Mom.emm$contrasts)
AD1m5Mom.contrasts$model <- "AD1m5Mom"

AD1m5Dad.emm <- emmeans(AD1m5Dad, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
AD1m5Dad.contrasts <- as.data.frame(AD1m5Dad.emm$contrasts)
AD1m5Dad.contrasts$model <- "AD1m5Dad"

AD1m5.emm <- emmeans(AD1m5, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
AD1m5.contrasts <- as.data.frame(AD1m5.emm$contrasts)
AD1m5.contrasts$model <- "AD1m5full"

AD1m5.tbl <- rbind(AD1m5.contrasts,AD1m5Mom.contrasts,AD1m5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") %>%
	filter(p.value<0.05)

write.table(AD1m5.tbl, file="AD1m5.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


AD1m5.all.tbl <- rbind(AD1m5.contrasts,AD1m5Mom.contrasts,AD1m5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") 

write.table(AD1m5.all.tbl, file="AD1m5.AllContrasts.tsv", quote=F, row.names=F, sep="\t")




##### construct reciprocal df for model 5, AD2 #####

deldataAD2 <- SSdataAD2 %>% 
	mutate(deltaMomDad = dipMom - dipDad) %>%
	mutate(deltaDadMom = dipDad - dipMom) %>%
	select(c(Targeting,dadGene,momGene,plantID,tetMom,tetDad,deltaMomDad,deltaDadMom))


##### run model 5, AD2 #####

AD2m5 <- lm(cbind(tetMom, tetDad) ~ Targeting*deltaMomDad, data = deldataAD2)
summary(AD2m5)

Plotm5AD2 <- emmip_ggplot(emmip(AD2m5, Targeting ~ deltaMomDad | subgenome, mult.name = "subgenome", at=list(deltaMomDad=seq(-10,10,by=2)),plotit=F))
goodPlotAD2 <- Plotm5AD2 + scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) 
ggsave("AD2m5.linearPrediction.jpg",goodPlotAD2,units="in", width=10, height=5, dpi=400)


# also run the models independently to access results easier
# will give same results
AD2m5Mom <- lm(tetMom ~ Targeting*deltaMomDad, data = deldataAD2)
AD2m5Dad <- lm(tetDad ~ Targeting*deltaMomDad, data = deldataAD2)


# check residual normality
deldataAD2$predictMom <- predict(AD2m5)[,1]
deldataAD2$predictDad <- predict(AD2m5)[,2]

deldataAD2$residMom <- deldataAD2$tetMom - deldataAD2$predictMom
deldataAD2$residDad <- deldataAD2$tetDad - deldataAD2$predictDad

deldataAD2$zresidMom <- (deldataAD2$residMom - mean(deldataAD2$residMom))/sd(deldataAD2$residMom)
deldataAD2$zresidDad <- (deldataAD2$residDad - mean(deldataAD2$residDad))/sd(deldataAD2$residDad)

histMom5AD2 <- hist(deldataAD2$zresidMom, plot=F)
histDad5AD2 <- hist(deldataAD2$zresidDad, plot=F)

## make a table of fixed effects
AD2m5.fixedEffects <- as.data.frame(tidy(AD2m5) %>%
	mutate(term= str_replace(term,"_Non-interacting","NI")) %>% 
	mutate(term= str_replace(term,"Targeting","")) %>% 
	mutate(term= str_replace(term,"Dual-targeted","D")) %>% 
	mutate(term= str_replace(term,"Mitochondria-targeted","M")) %>% 
	mutate(term= str_replace(term,"Plastid-targeted","P")) %>% 
	mutate(term= str_replace(term,"_Interacting","I"))
)

write.table(AD2m5.fixedEffects, file="AD2m5.fixedEffects.tsv", quote=F, row.names=F, sep="\t") 

AD2m5.fixedEffects.sig <- AD2m5.fixedEffects %>%
	filter(p.value<0.05)



## run a type III anova and make table of summary
AD2m5.anova <- Anova(AD2m5, type=3)
AD2m5Mom.anova <- anova_summary(Anova(AD2m5Mom, type=3))
AD2m5Dad.anova <- anova_summary(Anova(AD2m5Dad, type=3))

write.table(AD2m5Mom.anova, file="AD2m5.anovaSummary.Mom.tsv", quote=F, row.names=F, sep="\t")
write.table(AD2m5Dad.anova, file="AD2m5.anovaSummary.Dad.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
AD2m5Mom.emm <- emmeans(AD2m5Mom, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
AD2m5Mom.contrasts <- as.data.frame(AD2m5Mom.emm$contrasts)
AD2m5Mom.contrasts$model <- "AD2m5Mom"

AD2m5Dad.emm <- emmeans(AD2m5Dad, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
AD2m5Dad.contrasts <- as.data.frame(AD2m5Dad.emm$contrasts)
AD2m5Dad.contrasts$model <- "AD2m5Dad"

AD2m5.emm <- emmeans(AD2m5, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
AD2m5.contrasts <- as.data.frame(AD2m5.emm$contrasts)
AD2m5.contrasts$model <- "AD2m5full"

AD2m5.tbl <- rbind(AD2m5.contrasts,AD2m5Mom.contrasts,AD2m5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") %>%
	filter(p.value<0.05)

write.table(AD2m5.tbl, file="AD2m5.SigContrasts.tsv", quote=F, row.names=F, sep="\t")



AD2m5.all.tbl <- rbind(AD2m5.contrasts,AD2m5Mom.contrasts,AD2m5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") 

write.table(AD2m5.all.tbl, file="AD2m5.AllContrasts.tsv", quote=F, row.names=F, sep="\t")


#### lm with subgenome vs parent only ####

AD1m4Mom <- lm(tetMom ~ Targeting*dipMom, data = SSdataAD1)
AD1m4Dad <- lm(tetDad ~ Targeting*dipDad, data = SSdataAD1)

momPlot <- emmip_ggplot(emmip(AD1m4Mom, Targeting ~ dipMom, at=list(dipMom=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
momPlot <- momPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) + 
	theme(legend.position = "none") +
	ylim(-2,20)

dadPlot <- emmip_ggplot(emmip(AD1m4Dad, Targeting ~ dipDad, at=list(dipDad=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
dadPlot <- dadPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) + 
	theme(legend.position = "none") +
	ylim(-2,20) + 
	ylab("") +
	theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
	theme(plot.margin = unit(c(5,5,5,0),"pt"))


momdadPlot <- grid.arrange(momPlot,dadPlot,ncol=2)
ggsave("AD1m4sep.linearPrediction.jpg",momdadPlot,units="in", width=10, height=5, dpi=400)

#### lm with subgenome vs parent only ####

AD2m4Mom <- lm(tetMom ~ Targeting*dipMom, data = SSdataAD2)
AD2m4Dad <- lm(tetDad ~ Targeting*dipDad, data = SSdataAD2)

momPlot <- emmip_ggplot(emmip(AD2m4Mom, Targeting ~ dipMom, at=list(dipMom=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
momPlot <- momPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) +
	theme(legend.position = "none") +
	ylim(-2,20)
dadPlot <- emmip_ggplot(emmip(AD2m4Dad, Targeting ~ dipDad, at=list(dipDad=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
dadPlot <- dadPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) + 
	theme(legend.position = "none") +
	ylim(-2,20) + 
	ylab("") +
	theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
	theme(plot.margin = unit(c(5,5,5,0),"pt"))


momdadPlot <- grid.arrange(momPlot,dadPlot,ncol=2)
ggsave("AD2m4sep.linearPrediction.jpg",momdadPlot,units="in", width=10, height=5, dpi=400)



