library(data.table)
library(tidyverse)
setDTthreads(12)

cotton <- fread("Gossypium.mergedQuartets.one-file-to-rule-them-all.txt")

cotton <- cotton %>% 
	select(Targeting,paternalDiploid,paternalTetraploid,maternalTetraploid,maternalDiploid) %>%
	filter(!str_detect(paternalDiploid, ',')) %>%
	filter(!str_detect(paternalTetraploid, ',')) %>%
	filter(!str_detect(maternalTetraploid, ',')) %>%
	filter(!str_detect(maternalDiploid, ',')) %>%
	filter(str_detect(paternalDiploid, 'D5')) %>%
	filter(str_detect(maternalDiploid, 'G')) %>% 
	filter(startsWith(paternalTetraploid,"G") | startsWith(maternalTetraploid,"G")) %>%
	mutate(presentIn = case_when(startsWith(paternalTetraploid,"G") & startsWith(maternalTetraploid,"G") ~ "both", 
	startsWith(paternalTetraploid,"G") & startsWith(maternalTetraploid,"") ~ "dad",
	startsWith(paternalTetraploid,"") & startsWith(maternalTetraploid,"G") ~ "mom"))

cottonTable <- t(table(cotton$presentIn,cotton$Targeting))

for (i in c(1,2,3,4,6,7)) {
	ft <- fisher.test(cottonTable[c(i,5),c(2,3)])
	print(paste0(i," ", ft$p.value)) }



chenopodium <- fread("Chenopodium.mergedQuartets.one-file-to-rule-them-all.txt")

chenopodium <- chenopodium %>% 
	select(Targeting,paternalDiploid,PHY_paternal_homeolog_in_tet,PHY_maternal_homeolog_in_tet,maternalDiploid) %>%
	rename(paternalTetraploid=PHY_paternal_homeolog_in_tet) %>%
	rename(maternalTetraploid=PHY_maternal_homeolog_in_tet) %>%
	mutate(across(everything(), ~replace_na(.x, ""))) %>%
	filter(!str_detect(paternalDiploid, ',')) %>%
	filter(!str_detect(paternalTetraploid, ',')) %>%
	filter(!str_detect(maternalTetraploid, ',')) %>%
	filter(!str_detect(maternalDiploid, ',')) %>%
	filter(str_detect(paternalDiploid, 'Cs')) %>%
	filter(str_detect(maternalDiploid, 'C')) %>% 
	filter(startsWith(paternalTetraploid,"C") | startsWith(maternalTetraploid,"C")) %>%
	mutate(presentIn = case_when(startsWith(paternalTetraploid,"C") & startsWith(maternalTetraploid,"C") ~ "both", 
	startsWith(paternalTetraploid,"C") & startsWith(maternalTetraploid,"") ~ "dad",
	startsWith(paternalTetraploid,"") & startsWith(maternalTetraploid,"C") ~ "mom"))

chenopodiumTable <- t(table(chenopodium$presentIn,chenopodium$Targeting))


arachis <- fread("Arachis.mergedQuartets.one-file-to-rule-them-all.txt")

arachis <- arachis %>% 
	select(Targeting,paternalDiploid,PHY_paternal_homeolog_in_tet,PHY_maternal_homeolog_in_tet,maternalDiploid) %>%
	rename(paternalTetraploid=PHY_paternal_homeolog_in_tet) %>%
	rename(maternalTetraploid=PHY_maternal_homeolog_in_tet) %>%
	mutate(across(everything(), ~replace_na(.x, ""))) %>%
	filter(!str_detect(paternalDiploid, ',')) %>%
	filter(!str_detect(paternalTetraploid, ',')) %>%
	filter(!str_detect(maternalTetraploid, ',')) %>%
	filter(!str_detect(maternalDiploid, ',')) %>%
	filter(str_detect(paternalDiploid, 'a')) %>%
	filter(str_detect(maternalDiploid, 'a')) %>% 
	filter(startsWith(paternalTetraploid,"a") | startsWith(maternalTetraploid,"a")) %>%
	mutate(presentIn = case_when(startsWith(paternalTetraploid,"a") & startsWith(maternalTetraploid,"a") ~ "both", 
	startsWith(paternalTetraploid,"a") & startsWith(maternalTetraploid,"") ~ "dad",
	startsWith(paternalTetraploid,"") & startsWith(maternalTetraploid,"a") ~ "mom"))

arachisTable <- t(table(arachis$presentIn,arachis$Targeting))


arabidopsis <- fread("Arabidopsis.mergedQuartets.one-file-to-rule-them-all.txt")

arabidopsis <- arabidopsis %>% 
	select(Targeting,paternalDiploid,paternalTetraploid,maternalTetraploid,maternalDiploid) %>%
	filter(!str_detect(paternalDiploid, ',')) %>%
	filter(!str_detect(paternalTetraploid, ',')) %>%
	filter(!str_detect(maternalTetraploid, ',')) %>%
	filter(!str_detect(maternalDiploid, ',')) %>%
	filter(str_detect(paternalDiploid, 'A')) %>%
	filter(str_detect(maternalDiploid, 'A')) %>% 
	filter(startsWith(paternalTetraploid,"A") | startsWith(maternalTetraploid,"A")) %>%
	mutate(presentIn = case_when(startsWith(paternalTetraploid,"A") & startsWith(maternalTetraploid,"A") ~ "both", 
	startsWith(paternalTetraploid,"A") & startsWith(maternalTetraploid,"") ~ "dad",
	startsWith(paternalTetraploid,"") & startsWith(maternalTetraploid,"A") ~ "mom"))

arabidopsisTable <- t(table(arabidopsis$presentIn,arabidopsis$Targeting))

for (i in c(1,2,3,4,6,7)) {
	ft <- fisher.test(arabidopsisTable[c(i,5),c(2,3)])
	print(paste0(i," ", ft$p.value)) }

	
