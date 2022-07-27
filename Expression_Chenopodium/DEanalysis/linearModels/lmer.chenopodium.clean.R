setwd("W:/corrinne/Cytonuclear/chenopodium/DEanalysis/linearModels")

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
rlogfull <- fread("chenopodium.lm.tbl", header=T, sep="\t") %>% 
	filter(str_detect(gene,"AUR")) %>%
	filter(!is.na(subgenome)) %>%
	filter(!is.na(Cq_1)) %>%
	rename(Subgenome=subgenome) %>%
	rename(Targeting=category) %>% 
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-"))

rlog <- rlogfull %>% 
	select(!contains(c("Cp","Cs")))
		
# read in and number the pairs then extract into two 2-column df
pairs <- fread("chenopodium.gene.pairs", col.names=c("category","dad","mom"), sep="\t") 

## make a df that has diploids combined
rlogmom <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","Cp"))) %>% 
	setNames(gsub("Cp", "Dip", names(.))) %>%
	rename(Dip_3=Dip_7) %>%
	relocate(Dip_3, .after = Dip_2) %>%
	select(!Dip_1) 

rlogdad <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","Cs"))) %>% 
	setNames(gsub("Cs", "Dip", names(.)))

rlogdip <- bind_rows(rlogmom,rlogdad) %>% 
	right_join(.,rlog,by='gene') %>% 
	select(!Subgenome) 


##### count difference #####
# for uniformity, always mom - dad or mom/dad
# compared to momDip - dadDip or momDip/dadDip
homoeoPairs <- pairs[,2:3] %>% filter(dad %in% rlog$gene & mom %in% rlog$gene)

## make the delta
Drlog <- rlogdip[1:2,]
Drlog <- Drlog[-(1:2),]

for (i in 1:nrow(homoeoPairs)) {
    Drlog[nrow(Drlog)+1,] <- rbind(
        cbind(gene=paste0(homoeoPairs[i,2],"_",homoeoPairs[i,1]),
        rlogdip[rlogdip$gene==homoeoPairs[i,2], -c(1,11)] - 
        rlogdip[rlogdip$gene==homoeoPairs[i,1], -c(1,11)],
        rlogdip[rlogdip$gene==homoeoPairs[i,1], c(11)]))
}

Ddata <- Drlog %>% 
	select(!contains("CyMira")) %>%
	unite("Plant2", c(Dip_2,Cq_1), sep="X") %>% 
	unite("Plant3", c(Dip_3,Cq_3), sep="X") %>% 
	unite("Plant4", c(Dip_4,Cq_4), sep="X") %>% 
	unite("Plant5", c(Dip_5,Cq_5), sep="X") %>%
	select(!Cq_6) %>%	 
	pivot_longer(cols=starts_with("P"), names_to="plantID", values_to="DipXTetlog") %>% 
	separate(DipXTetlog, c("dipDelta","tetDelta"), sep="X") %>% 
	mutate(across(ends_with("Delta"), as.double))

fwrite(Ddata, file="Cq.delta.tsv", quote=F, row.names=F, sep="\t")




##### run model 2, Asue #####
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
	separate(contrast,c("Category1","Category2"),sep=" ") 

write.table(Cqm2D.tbl, file="Cqm2.Delta.AllContrasts.tsv", quote=F, row.names=F, sep="\t")


##### run model 3, Cq #####
Cqm3D <- lmer(tetDelta ~ Targeting*dipDelta + (1 | plantID), data=Ddata, REML=T)
summary(Cqm3D)

Cqm3D <- lm(tetDelta ~ Targeting*dipDelta, data=Ddata)
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

## perform contrasts and generate table
Cqm3D.emm <- emmeans(Cqm3D, pairwise ~ Targeting, lmer.df='satterthwaite',adjust='Holm')
Cqm3D.contrasts <- as.data.frame(Cqm3D.emm$contrasts)

Cqm3D.tbl <- Cqm3D.contrasts %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1","Category2"),sep=" ") # %>%

write.table(Cqm3D.tbl, file="Cqm3.Delta.AllContrasts.tsv", quote=F, row.names=F, sep="\t")


lrtest(Cqm3D, Cqm2D)

#####
# let's model tetMom, tetDad separately 
#####


SSrlogmom <- rlogfull %>% filter(Subgenome=="mom") %>% 
	select(contains(c("gene","Cp","Cq"))) %>% 
	setNames(gsub("Cp", "dipMom", names(.))) %>% 
	setNames(gsub("Cq", "tetMom", names(.))) %>% 
	rename(momGene = gene)

SSrlogdad <- rlogfull %>% filter(Subgenome=="dad") %>% 
	select(contains(c("gene","Cs","Cq"))) %>% 
	setNames(gsub("Cs", "dipDad", names(.))) %>% 
	setNames(gsub("Cq", "tetDad", names(.))) %>% 
	rename(dadGene = gene)


SSrlog <- pairs %>% 
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
	mutate(Targeting = str_replace(Targeting,"Not-","AreNot-"))

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




##### construct reciprocal df for model 5, Cq #####

deldata <- SSdata %>% 
	mutate(deltaMomDad = dipMom - dipDad) %>%
	mutate(deltaDadMom = dipDad - dipMom) %>%
	select(c(Targeting,dadGene,momGene,plantID,tetMom,tetDad,deltaMomDad,deltaDadMom))


##### run model 5, Cq #####

Cqm5 <- lm(cbind(tetMom, tetDad) ~ Targeting*deltaMomDad, data = deldata)
summary(Cqm5)

Plotm5 <- emmip_ggplot(emmip(Cqm5, Targeting ~ deltaMomDad | subgenome, mult.name = "subgenome", at=list(deltaMomDad=seq(-10,10,by=2)),plotit=F))
goodPlot <- Plotm5 + scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) 
ggsave("Cqm5.linearPrediction.jpg",goodPlot,units="in", width=10, height=5, dpi=400)



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


Cqm5.all.tbl <- rbind(Cqm5.contrasts,Cqm5Mom.contrasts,Cqm5Dad.contrasts) %>% 
	mutate(contrast = str_replace(contrast,"\\) - \\("," ")) %>%
	mutate(contrast = str_replace_all(contrast,"[()]","")) %>% 
	separate(contrast,c("Category1",NA,"Category2",NA),sep=" ") 

write.table(Cqm5.all.tbl, file="Cqm5.AllContrasts.tsv", quote=F, row.names=F, sep="\t")



#### lm with subgenome vs parent only ####


Cqm4Mom <- lm(tetMom ~ Targeting*dipMom, data = SSdata)
Cqm4Dad <- lm(tetDad ~ Targeting*dipDad, data = SSdata)

momPlot <- emmip_ggplot(emmip(Cqm4Mom, Targeting ~ dipMom, at=list(dipMom=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
momPlot <- momPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) +
	theme(legend.position = "none") +
	ylim(-2,20)
dadPlot <- emmip_ggplot(emmip(Cqm4Dad, Targeting ~ dipDad, at=list(dipDad=c(0,2,4,6,8,10,12,14,16,18,20)),plotit=F))
dadPlot <- dadPlot + 
	scale_color_manual(labels = c("NOT", "DI","DNI","MI","MNI","PI","PNI"),values=c("black","mediumblue","skyblue1","darkred","lightcoral","darkgreen","palegreen3")) + 
	theme(legend.position = "none") +
	ylim(-2,20) + 
	ylab("") +
	theme(axis.text.y = element_blank(),axis.ticks.y = element_blank()) +
	theme(plot.margin = unit(c(5,5,5,0),"pt"))

momdadPlot <- grid.arrange(momPlot,dadPlot,ncol=2)
ggsave("Cqm4sep.linearPrediction.jpg",momdadPlot,units="in", width=10, height=5, dpi=400)





