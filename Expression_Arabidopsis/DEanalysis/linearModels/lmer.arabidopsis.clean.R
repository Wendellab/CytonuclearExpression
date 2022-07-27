setwd("w:/corrinne/Cytonuclear/arabidopsis/DEanalysis12kV2/linearModels")

library(tidyverse)
library(lme4)
library(car)
library(sjPlot)
library(rstatix)
library(emmeans)
library(data.table)
library(gridExtra)

##### global options #####
emm_options(lmerTest.limit = 88020)
setDTthreads(10)
options(datatable.fread.datatable=FALSE)

# read and filter the table
rlogfull <- fread("arabidopsis.lm.tbl", header=T, sep="\t") %>% 
	filter(str_detect(gene,"As")) %>%
	filter(!is.na(subgenome)) %>%
	filter(!is.na(Asue_1)) %>% 
	dplyr::rename(Subgenome=subgenome) %>%
	dplyr::rename(Targeting=category)%>%
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-"))%>%
	select(!contains("CyMira")) # put NOT first alphabetically


rlog <- rlogfull %>% 
	select(!contains(c("Aare","Atha"))) 

# read in and number the pairs then extract into two 2-column df
pairs <- fread("arabidopsis.gene.pairs", col.names=c("category","dad","mom"), sep="\t") %>% rowid_to_column("genePairID")

## make a df that has diploids combined
rlogmom <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","Atha"))) %>% 
	setNames(gsub("Atha", "Dip", names(.)))

rlogdad <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","Aare"))) %>% 
	setNames(gsub("Aare", "Dip", names(.)))

rlogdip <- bind_rows(rlogmom,rlogdad) %>% 
	right_join(.,rlog,by='gene') %>%
	select(!c("Dip_6","Dip_7","Subgenome")) 


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
        rlogdip[rlogdip$gene==homoeoPairs[i,2], -c(1,12)] - 
        rlogdip[rlogdip$gene==homoeoPairs[i,1], -c(1,12)],
        rlogdip[rlogdip$gene==homoeoPairs[i,1], c(12)]))
}

Ddata <- Drlog %>% 
	unite("Plant1", c(Dip_1,Asue_1), sep="X") %>% 
	unite("Plant2", c(Dip_2,Asue_2), sep="X") %>% 
	unite("Plant3", c(Dip_3,Asue_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,Asue_4), sep="X") %>% 
	unite("Plant5", c(Dip_5,Asue_5), sep="X") %>% 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipDelta","tetDelta"), sep="X") %>%
	mutate(across(ends_with("Delta"), as.double))

fwrite(Ddata, file="As.delta.tbl", quote=F, row.names=F, sep="\t")


##### run model 2, Asue #####

### run mixed model, delta ###
Asm2D <- lmer(tetDelta ~ Targeting + (1 | plantID), data=Ddata, REML=T)
summary(Asm2D)

Asm2D <- lm(tetDelta ~ Targeting, data=Ddata)
summary(Asm2D)

## save the model plot
jpeg(file="Asm2.Delta.modelPlot.jpeg")
plot(Asm2D)
dev.off()

## make a table of fixed effects
Asm2D.fixedEffects <- summary(Asm2D)$coefficients
write.table(Asm2D.fixedEffects, file="Asm2.Delta.fixedEffects.tbl", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Asm2D, type=2)
Asm2D.anova <- anova_summary(Anova(Asm2D, type=3))
write.table(Asm2D.anova, file="Asm2.Delta.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## perform contrasts and generate table
Asm2D.emm <- emmeans(Asm2D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Asm2D.contrasts <- as.data.frame(Asm2D.emm$contrasts)

Asm2D.tbl <- Asm2D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") # %>%
#	filter( ((Parent1 == Parent2 & Category1 != Category2) | (Parent1 != Parent2 & Category1 == Category2)) & p.value<0.05)

write.table(Asm2D.tbl, file="Asm2.Delta.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


##### run model 3, Asue #####

### run mixed model 3, delta ###
Asm3D <- lmer(tetDelta ~ Targeting*dipDelta + (1 | plantID), data=Ddata, REML=T)
summary(Asm3D)

Asm3D <- lm(tetDelta ~ Targeting*dipDelta, data=Ddata)
summary(Asm3D)


## save the model plot
jpeg(file="Asm3.Delta.modelPlot.jpeg")
plot(Asm3D)
dev.off()

## make a table of fixed effects
Asm3D.fixedEffects <- summary(Asm3D)$coefficients
write.table(Asm3D.fixedEffects, file="Asm3.Delta.fixedEffects.tsv", quote=F, row.names=T, sep="\t") 

## run a type III anova and make table of summary
Anova(Asm3D, type=2)
Asm3D.anova <- anova_summary(Anova(Asm3D, type=3))
write.table(Asm3D.anova, file="Asm3.Delta.anovaSummary.tsv", quote=F, row.names=F, sep="\t")

## graph fixed effects
emmip(Asm3D, dipDelta ~ Targeting, cov.reduce=F)

## perform contrasts and generate table
Asm3D.emm <- emmeans(Asm3D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Asm3D.contrasts <- as.data.frame(Asm3D.emm$contrasts)

Asm3D.tbl <- Asm3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") # %>%
	#filter(p.value<0.05)

write.table(Asm3D.tbl, file="Asm3.Delta.SigContrasts.tsv", quote=F, row.names=F, sep="\t")


#####
# let's model tetMom, tetDad separately together for Model 4
#####


SSrlogmom <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","Atha","Asue"))) %>% 
	setNames(gsub("Atha", "dipMom", names(.))) %>% 
	setNames(gsub("Asue", "tetMom", names(.))) %>% 
	rename(momGene = gene)

SSrlogdad <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","Aare","Asue"))) %>% 
	setNames(gsub("Aare", "dipDad", names(.))) %>% 
	setNames(gsub("Asue", "tetDad", names(.))) %>% 
	rename(dadGene = gene)


SSrlog <- pairs %>% select(!genePairID) %>%
	rename(Targeting = category) %>%
	rename(dadGene = dad) %>%
	rename(momGene = mom)  %>% 
	left_join(.,SSrlogmom,by='momGene') %>% 
	left_join(.,SSrlogdad,by='dadGene') %>%
	filter(!(is.na(dipMom_1) | is.na(dipDad_1))) %>%
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-"))



SSdata <- SSrlog %>% 
	unite("PlantA", c(dipMom_1, dipDad_1, tetMom_1, tetDad_1), sep="X") %>% 
	unite("PlantB", c(dipMom_2, dipDad_2, tetMom_2, tetDad_2), sep="X") %>% 
	unite("PlantC", c(dipMom_3, dipDad_3, tetMom_3, tetDad_3), sep="X") %>% 
	unite("PlantD", c(dipMom_4, dipDad_4, tetMom_4, tetDad_4), sep="X") %>% 
	unite("PlantE", c(dipMom_5, dipDad_5, tetMom_5, tetDad_5), sep="X") %>%
	select(!c(dipDad_7,dipDad_6)) %>%	 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipMom","dipDad","tetMom","tetDad"), sep="X") %>% 
	mutate(across(ends_with("Mom"), as.double)) %>% 
	mutate(across(ends_with("Dad"), as.double))

fwrite(SSdata, file="As.SS.tsv", quote=F, row.names=F, sep="\t")

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
	   ggtitle("Arabidopsis suecica") + 
	   theme(plot.title = element_text(face = "italic")) + 
	   labs(x="Maternal Homoeolog Expression (rlog)", y="Paternal Homoeolog Expression (rlog)")

ggsave("AsMomvsDad.jpg",SSplot)





##### construct df for model 5, As #####

deldata <- SSdata %>% 
	mutate(deltaMomDad = dipMom - dipDad) %>%
	mutate(deltaDadMom = dipDad - dipMom) %>%
	select(c(Targeting,dadGene,momGene,plantID,tetMom,tetDad,deltaMomDad,deltaDadMom))


##### run model 5, As #####

Asm5 <- lm(cbind(tetMom, tetDad) ~ Targeting*deltaMomDad, data = deldata)
summary(Asm5)

Plotm5 <- emmip_ggplot(emmip(Asm5, Targeting ~ deltaMomDad | subgenome, mult.name = "subgenome", at=list(deltaMomDad=seq(-10,10,by=2)),plotit=F))
goodPlot <- Plotm5 + scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) 
ggsave("Asm5.linearPrediction.jpg",goodPlot,units="in", width=10, height=5, dpi=400)

emmip_ggplot(emmip(Asm5, Targeting ~ deltaMomDad | subgenome, mult.name = "subgenome", at=list(deltaMomDad=seq(-10,10,by=2))))



# also run the models independently to access results easier
# will give same results
Asm5Mom <- lm(tetMom ~ Targeting*deltaMomDad, data = deldata)
Asm5Dad <- lm(tetDad ~ Targeting*deltaMomDad, data = deldata)


# check residual normality
deldata$predictMom <- predict(Asm5)[,1]
deldata$predictDad <- predict(Asm5)[,2]

deldata$residMom <- deldata$tetMom - deldata$predictMom
deldata$residDad <- deldata$tetDad - deldata$predictDad

deldata$zresidMom <- (deldata$residMom - mean(deldata$residMom))/sd(deldata$residMom)
deldata$zresidDad <- (deldata$residDad - mean(deldata$residDad))/sd(deldata$residDad)

histMom5 <- hist(deldata$zresidMom, plot=F)
histDad5 <- hist(deldata$zresidDad, plot=F)


## make a table of fixed effects
Asm5.fixedEffects <- as.data.frame(tidy(Asm5) %>%
	mutate(term= str_replace(term,"_Non-interacting","NI")) %>% 
	mutate(term= str_replace(term,"Targeting","")) %>% 
	mutate(term= str_replace(term,"Dual-targeted","D")) %>% 
	mutate(term= str_replace(term,"Mitochondria-targeted","M")) %>% 
	mutate(term= str_replace(term,"Plastid-targeted","P")) %>% 
	mutate(term= str_replace(term,"_Interacting","I"))
)

write.table(Asm5.fixedEffects, file="Asm5.fixedEffects.tsv", quote=F, row.names=F, sep="\t") 

Asm5.fixedEffects.sig <- Asm5.fixedEffects %>%
	filter(p.value<0.05)



## run a type III anova and make table of summary
Asm5.anova <- Anova(Asm5, type=3)
Asm5Mom.anova <- anova_summary(Anova(Asm5Mom, type=3))
Asm5Dad.anova <- anova_summary(Anova(Asm5Dad, type=3))

write.table(Asm5Mom.anova, file="Asm5.anovaSummary.Mom.tsv", quote=F, row.names=F, sep="\t")
write.table(Asm5Dad.anova, file="Asm5.anovaSummary.Dad.tsv", quote=F, row.names=F, sep="\t")


## perform contrasts and generate table
Asm5Mom.emm <- emmeans(Asm5Mom, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Asm5Mom.contrasts <- as.data.frame(Asm5Mom.emm$contrasts)
Asm5Mom.contrasts$model <- "Asm5Mom"

Asm5Dad.emm <- emmeans(Asm5Dad, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Asm5Dad.contrasts <- as.data.frame(Asm5Dad.emm$contrasts)
Asm5Dad.contrasts$model <- "Asm5Dad"

Asm5.emm <- emmeans(Asm5, pairwise ~ Targeting*deltaMomDad, lmer.df='satterthwaite',adjust='Holm')
Asm5.contrasts <- as.data.frame(Asm5.emm$contrasts)
Asm5.contrasts$model <- "Asm5full"

Asm5.tbl <- rbind(Asm5.contrasts,Asm5Mom.contrasts,Asm5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") #%>%
	#filter(p.value<0.05)

write.table(Asm5.tbl, file="Asm5.SigContrasts.tsv", quote=F, row.names=F, sep="\t")





#### lm with subgenome vs parent only ####

Asm4Mom <- lm(tetMom ~ Targeting*dipMom, data = SSdata)
Asm4Dad <- lm(tetDad ~ Targeting*dipDad, data = SSdata)

momPlot <- emmip_ggplot(emmip(Asm4Mom, Targeting ~ dipMom, at=list(dipMom=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
momPlot <- momPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) +
	theme(legend.position = "none") +
	ylim(-3,20)
dadPlot <- emmip_ggplot(emmip(Asm4Dad, Targeting ~ dipDad, at=list(dipDad=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
dadPlot <- dadPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) + 
	theme(legend.position = "none") +
	ylim(-3,20) + 
	ylab("") +
	theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
	theme(plot.margin = unit(c(5,5,5,0),"pt"))

momdadPlot <- grid.arrange(momPlot,dadPlot,ncol=2)
ggsave("Asm4sep.linearPrediction.jpg",momdadPlot,units="in", width=10, height=5, dpi=400)




