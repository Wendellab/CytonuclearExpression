setwd("w:/corrinne/Cytonuclear/arachis/DEanalysis/linearModels")

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
his
fwrite(DdataAid, file="Aid.delta.tsv", quote=F, row.names=F, sep="\t")




##### run model 2, Ah #####


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


lrtest(Ahm3D,Ahm2D)
lrtest(Aidm3D,Aidm2D)


#####
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





##### construct reciprocal df for model 5, Ah #####

deldataAh <- SSdataAh %>% 
	mutate(deltaMomDad = dipMom - dipDad) %>%
	select(c(Targeting,dadGene,momGene,plantID,tetMom,tetDad,deltaMomDad))


##### run model 5, Ah #####

Ahm5 <- lm(cbind(tetMom, tetDad) ~ Targeting*deltaMomDad, data = deldataAh)
summary(Ahm5)

Plotm5Ah <- emmip_ggplot(emmip(Ahm5, Targeting ~ deltaMomDad | subgenome, mult.name = "subgenome", at=list(deltaMomDad=seq(-10,10,by=2)),plotit=F))
goodPlotAh <- Plotm5Ah + scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) 
ggsave("Ahm5.linearPrediction.jpg",goodPlotAh,units="in", width=10, height=5, dpi=400)


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
goodPlotAid <- Plotm5Aid + scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) 
ggsave("Aidm5.linearPrediction.jpg",goodPlotAid,units="in", width=10, height=5, dpi=400)


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
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ")

write.table(Aidm5.all.tbl, file="Aidm5.AllContrasts.tsv", quote=F, row.names=F, sep="\t")





#### lm with subgenome vs parent only ####

Ahm4Mom <- lm(tetMom ~ Targeting*dipMom, data = SSdataAh)
Ahm4Dad <- lm(tetDad ~ Targeting*dipDad, data = SSdataAh)

momPlot <- emmip_ggplot(emmip(Ahm4Mom, Targeting ~ dipMom, at=list(dipMom=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
momPlot <- momPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) +
	theme(legend.position = "none") +
	ylim(-3,22)
dadPlot <- emmip_ggplot(emmip(Ahm4Dad, Targeting ~ dipDad, at=list(dipDad=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
dadPlot <- dadPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) + 
	theme(legend.position = "none") +
	ylim(-3,22) + 
	ylab("") +
	theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
	theme(plot.margin = unit(c(5,5,5,0),"pt"))

momdadPlot <- grid.arrange(momPlot,dadPlot,ncol=2)
ggsave("Ahm4sep.linearPrediction.jpg",momdadPlot,units="in", width=10, height=5, dpi=400)



#### lm with subgenome vs parent only ####

Aidm4Mom <- lm(tetMom ~ Targeting*dipMom, data = SSdataAid)
Aidm4Dad <- lm(tetDad ~ Targeting*dipDad, data = SSdataAid)

momPlot <- emmip_ggplot(emmip(Aidm4Mom, Targeting ~ dipMom, at=list(dipMom=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
momPlot <- momPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) +
	theme(legend.position = "none") +
	ylim(-3,22)
dadPlot <- emmip_ggplot(emmip(Aidm4Dad, Targeting ~ dipDad, at=list(dipDad=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
dadPlot <- dadPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) + 
	theme(legend.position = "none") +
	ylim(-3,22) + 
	ylab("") +
	theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
	theme(plot.margin = unit(c(5,5,5,0),"pt"))

momdadPlot <- grid.arrange(momPlot,dadPlot,ncol=2)
ggsave("Aidm4sep.linearPrediction.jpg",momdadPlot,units="in", width=10, height=5, dpi=400)





