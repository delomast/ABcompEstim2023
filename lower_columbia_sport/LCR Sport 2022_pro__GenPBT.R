# devtools::install_github("delomast/fishCompTools")

library(fishCompTools)
library(dplyr)
library(sqldf)

## Steelhead harvest stock composition for LCR Sport 2022. Run 5/3/23.
## set catch file = 1  to generate proportions

## only AD clipped fish could be retained.
## R script modified for harvest analysis.  Using GenPBT_Rgroup assignments. 
## Heir varibles = pbt_Rgroup and Stock. script codes GenPBT_Rgroup == pbt_Rgroup
## 10000 iterations for 90% CI

rm(list=ls())      
start.time<-Sys.time()


setwd("C:/Users/abyrne/Documents/TAC/Steelhead PBT-GSI sampling/2022 Fishery/LCR Sport/Output proportion/")

##set CI intervals for all estimates and marginalized/summary estimates (this is alph in Scobi script). Default is 90% CIs.
 ci = .9  # it gets converted to alph value as alph = (1 - ci) for script
  a_1 = 1 - ci # and set alph = a_1 in script

# input sample biodata and harvest estimate 
fishData <- read.csv("LCR Sport 2022 data.csv", stringsAsFactors = F)  # sample biodata
fishData <- fishData[!is.na(fishData$DataKey),]
head(fishData)
tail(fishData)

#clipharvest <- read.csv("LCR Sport 20xx Clip harvest.csv", stringsAsFactors = F) # use when there is clip and unclip harvest (samples)
  #unclipharvest <- read.csv("LCR Sport catch xxx.csv", stringsAsFactors = F)  # use when there is clip and unclip harvest (samples)
  allharvest<- read.csv("LCR Sport 2022 catch for pro.csv", stringsAsFactors = F)  ## use when there is only clipped harvest (samples)
  
# specify tagrate file
tagrate_info <- read.csv("C:/Users/abyrne/Documents/TAC/Steelhead PBT-GSI sampling/Tag_rate_8-11-22.csv", stringsAsFactors = F)

### look for physical tags in Comments field and a stubby dorsal in the StubbyDorsal field (for BONAFF data only) 

fishData$PhysTag <- "notag"
fishData$PhysTag[grepl("RV|LV|MR|ML|LP|RP|RM|LM|CW|Q1|Wire", fishData$Comments) | grepl("YES|Yes|yes", fishData$StubbyDorsal)] <- "tag"

#recode PBT variable appropriately 

fishData$pbt_Rgroup <- fishData$GenPBT_Rgroup  

fishData$pbt_Rgroup[fishData$pbt_Rgroup %in% c("Unassigned-H", "Unassigned-W","Unassigned","Unassign")] <- "Unassigned"

fishData$pbt_Rgroup[fishData$pbt_Rgroup %in% c("", "NG","NA")] <- NA

fishData$GenPBT_Rgroup[fishData$GenPBT_Rgroup %in% c("Unassigned-H", "Unassigned-W","Unassigned","Unassign")] <- "Unassigned"
  fishData$GenPBT_Rgroup[fishData$GenPBT_Rgroup %in% c("", "NG","NA")] <- NA

 unique(fishData$pbt_Rgroup)
  sum(is.na(fishData$pbt_Rgroup))

##classify fish as Large (B-Index) or Small (A-Index). 
 fishData$Size<-ifelse(fishData$Length<780,"Small", "Large")

 #fishData_clip<-fishData[!(fishData$Adipose=="AI"),]  ## remove AI samples--in tribal harvest only.
  #fishData_unclip<-fishData[!(fishData$Adipose=="AD"),]  ## remove AD samples--for tribal unclip analysis only

 #Generate a list of all distinct pbt_Rgroups with Stock, BY, Basin
 list<-distinct(data.frame(fishData$pbt_Rgroup,fishData$Stock,fishData$BY,fishData$Basin))
   colnames(list) <- c("pbt_Rgroup", "Stock", "BY","Basin")
 
 set.seed(9)   # so if you need to re-run this script you'll get same CIs from the bootstrap
    
## do for clipped fish (H) by pbt_Rgroup, Stock
 
SCOBI_deux_fast(adultData = fishData, windowData = allharvest,
Run = "LCR Sport 2022_pro_H", RTYPE = "clipped", Hierarch_variables = c("pbt_Rgroup","Stock"), 
SizeCut = NULL, alph = a_1, B = 10000, writeBoot = T, pbtRates = tagrate_info, 
adClipVariable = "Adipose", physTagsVariable = "PhysTag", pbtGroupVariable = "pbt_Rgroup", dataGroupVariable = "week", screenOutput = "LCR Sport 2022_pro_H.txt",spibetr=TRUE)

## Set output file names for summaries--add other name for HNC and Wild if unclipped samples present)
    Run = "LCR Sport 2022_pro_H.csv"  # for naming csv H clipped file with totals and CIs from AB sum Loops
   
 
##### Loop 1 - get total and CI by Assign type (all PBT and all GSI)

Results <- read.delim("LCR Sport 2022_pro_H_Estim_Grand_Totals_Hier_Stock.txt", header=TRUE)
 Results$assigntype<-ifelse(Results$pbt_Rgroup=="Unassigned","GSI assigned","PBT assigned") 
 
 Results_by_assigntype <- Results %>%
  group_by(assigntype) %>%
  summarise(Total = sum(Estimated_total_for_run))

# Now get CIs for all PBT assigned and all GSI assigned 

Boots <- read.delim("LCR Sport 2022_pro_H_Boot_Hier_Stock.txt", header=TRUE)  
 Boots$assigntype<-Results$assigntype
  Boots.assign.type = Boots %>% group_by(assigntype) %>% summarise_all(funs(sum))
   Boots.for.CI <- subset(Boots.assign.type, select = -assigntype)

# Write totals and CI to dataframe for all PBT assigned and all GSI assigned (PBT Unassigned)
assigntype.lci.uci<- data.frame(Assigntype=Boots.assign.type$assigntype,Estimate=Results_by_assigntype$Total,lci=apply(Boots.for.CI,1,quantile,(1-ci)/2),uci=apply(Boots.for.CI,1,quantile,((1-ci)/2)+ci), mean=apply(Boots.for.CI,1,mean),st.dev=apply(Boots.for.CI,1,sd))

write.csv(assigntype.lci.uci, file = paste("Results by Assign Type_",Run,sep="_"))

##### Loop 2 -- totals and CIs by Basin of PBT assigned stocks only
  
# Remove stock duplicate names in list
list1<-list[!duplicated(list$Stock), ]

#add Basin to Grand total dataframe
Basinpbt = sqldf("Select Results.Stock,list1.Basin,Results.assigntype,Results.Estimated_total_for_run From
                 Results INNER JOIN list1 ON(Results.Stock = list1.Stock)")

# add Basin to bootstrap iteration datframe
  Boots$Basin <-Basinpbt$Basin

#Remove the GSI assigned and only retain the PBT assigned
Basinpbt1<-Basinpbt[Basinpbt$assigntype != "GSI assigned", ] 
 Ci.basin.pbt<-Boots[Boots$assign != "GSI assigned", ] 


#Sum totals by Basin 
Totals_by_basin.pbtonly <- Basinpbt1 %>%
  group_by(Basin) %>%
  summarise(Total = sum(Estimated_total_for_run))

# CIs by Basin
Ci.basin.pbt<- subset(Ci.basin.pbt, select = -assigntype)    
 Ci.basin.pbt2 <- Ci.basin.pbt %>% group_by(Basin) %>% summarise_each(funs(sum))
  Ci.basin.pbt2<- subset(Ci.basin.pbt2, select = -Basin)

#Write Basin totals and CI to dataframe for PBT assigned only 

basin.pbt.lci.uci<- data.frame(Basin=Totals_by_basin.pbtonly$Basin,Estimate=Totals_by_basin.pbtonly$Total,lci=apply(Ci.basin.pbt2,1,quantile,(1-ci)/2),uci=apply(Ci.basin.pbt2,1,quantile,((1-ci)/2)+ci), mean=apply(Ci.basin.pbt2,1,mean),st.dev=apply(Ci.basin.pbt2,1,sd))

write.csv(basin.pbt.lci.uci, file = paste("Results by Basin_H_PBT only_",Run,sep="_"))
  
##### Loop 3 -- totals by Stock and BY for PBT assigned pbt_Rgroups 
  
  stocktotal <- read.delim("LCR Sport 2022_pro_H_Estim_Grand_Totals_Hier_pbt_Rgroup.txt", header=TRUE) ## use totals
   stock.boots <- read.delim("LCR Sport 2022_pro_H_Boot_Hier_pbt_Rgroup.txt", header=TRUE) ## iterations for CI
    stock.boots$pbt_Rgroup = stocktotal$pbt_Rgroup
	
 stocktotal<-stocktotal[stocktotal$pbt_Rgroup != "Unassigned", ] # remove GSI assigned, keep only PBT assigned
   stock.boots<-stock.boots[stock.boots$pbt_Rgroup != "Unassigned", ] # remove GSI assigned, keep only PBT assigned
   
   list_pbt<-list[list$pbt_Rgroup != "Unassigned", ] # retain only the PBT assigned groups
   
# add Stock and BY to dataframes
   
   stocktotal_pbt = sqldf("Select stocktotal.pbt_Rgroup,list_pbt.Stock,list_pbt.BY,stocktotal.Estimated_total_for_run From
   stocktotal INNER JOIN list_pbt ON(stocktotal.pbt_Rgroup = list_pbt.pbt_Rgroup)")
  
  stock.boots$Stock <-stocktotal_pbt$Stock
   stock.boots$BY <-stocktotal_pbt$BY
 
 # total by Stock and BY 
  Totals_by_stock.by <- stocktotal_pbt %>%
    group_by(Stock,BY) %>%
    summarise(Total = sum(Estimated_total_for_run))
  
# get CIs for Stock and BY
 
  stock.boots<-subset(stock.boots,select =  - pbt_Rgroup)
    stock.boots.ci = stock.boots %>% group_by(Stock,BY) %>% summarise_all(funs(sum))
     stock.boots.ci2<- subset(stock.boots.ci, select = c(-Stock, -BY)) 
  
  #Write totals and CI to dataframe for all PBT assigned by Stock and BY
  stock.by.pbt.lci.uci<- data.frame(Stock=Totals_by_stock.by$Stock,BY = Totals_by_stock.by$BY,Estimate=Totals_by_stock.by$Total,lci=apply(stock.boots.ci2,1,quantile,(1-ci)/2),uci=apply(stock.boots.ci2,1,quantile,((1-ci)/2)+ci), mean=apply(stock.boots.ci2,1,mean),st.dev=apply(stock.boots.ci2,1,sd))
  
  write.csv(stock.by.pbt.lci.uci, file = paste("Results by Hatchery stock and BY_",Run,sep="_"))
  
 ##### Loop 3.a  -- get totals and Cis by Stock (add all BYs)
  
  marginalize_SD(estimatesFile="LCR Sport 2022_pro_H_Estim_Grand_Totals_Hier_Stock.txt",bootHier="LCR Sport 2022_pro_H_Boot_Hier_Stock.txt",marginalize = c("pbt_Rgroup"),alph=a_1)
 
##### Loop 4  -- Totals and CIs by Basin for all samples 
  
stockbasin<- read.delim("LCR Sport 2022_pro_H_Estim_Grand_Totals_Hier_Stock.txt", header=TRUE) ## use totals
 stockbasinci<- read.delim("LCR Sport 2022_pro_H_Boot_Hier_Stock.txt", header=TRUE) ## for Cis
  stockbasin1<-subset(stockbasin,select = -pbt_Rgroup)
 
#add Basin to Grand total dataframe
  stockbasin2 = sqldf("Select stockbasin1.Stock,list1.Basin,stockbasin1.Estimated_total_for_run From
   stockbasin1 INNER JOIN list1 ON(stockbasin1.Stock= list1.Stock)")
  
#add Basin to Bootstrap iterations dataframe
   stockbasinci$Basin<-stockbasin2$Basin
  
# Get totals by Basin
  allBasin.total <- stockbasin2 %>%
    group_by(Basin) %>%
    summarise(Total = sum(Estimated_total_for_run))
  
##Now get CIs by Basin

 basin.ci1 <- stockbasinci %>% 
    group_by(Basin) %>% 
    summarise_each(funs(sum))
  
  for.basin.ci <- subset(basin.ci1, select = -Basin)

#Write totals and CI to dataframe for all PBT assigned by Stock and BY

  basin.all.lci.uci<- data.frame(Basin=allBasin.total$Basin,Estimate=allBasin.total$Total,lci=apply(for.basin.ci,1,quantile,(1-ci)/2),uci=apply(for.basin.ci,1,quantile,((1-ci)/2)+ci), mean=apply(for.basin.ci,1,mean),st.dev=apply(for.basin.ci,1,sd))
 
   write.csv(basin.all.lci.uci, file = paste("Results by Basin all samples_",Run,sep="_"))

## HNC and Wild routines not needed since only have AD-clipped fish 
  
##do for unclipped HNC fish.  
#SCOBI_deux(adultData = fishData, windowData = allharvest,
#Run = "Pound Net 2020_HNC", RTYPE = "noclip_H", Hierarch_variables = c("pbt_Rgroup","Stock"), SizeCut = NULL, alph = a_1, B = 10, writeBoot = T, pbtRates = tagrate_info, adClipVariable = "Adipose", physTagsVariable = "PhysTag", pbtGroupVariable = "pbt_Rgroup", dataGroupVariable = "week", screenOutput = "Pound Net 2020_HNC.txt",spibetr=TRUE)

##### add summary Loops from H clipped for HNC (to get same info for PBT assigned unclipped fish as for PBT assigned clipped fish)
##### when there are unclipped samples being analyzed
   
##do for unclipped WILD fish.  
#SCOBI_deux(adultData = fishData, windowData = allharvest,
#Run = "Pound Net 2020_Wild", RTYPE = "wild", Hierarch_variables = c("GenStock","Basin"), SizeCut = NULL, alph = a_1, B = 10000, writeBoot = T, pbtRates = tagrate_info, adClipVariable = "Adipose", physTagsVariable = "PhysTag", pbtGroupVariable = "pbt_Rgroup", dataGroupVariable = "week", screenOutput = "Pound Net 2020_Wild.txt",spibetr=TRUE)

# total and CIs for each wild Genstock from output Grand Estim and CI files for Wild


end.time<-Sys.time()
 runtime<-end.time-start.time
   runtime

   
##Valid Rtype = "clipped", "noclip_H", "wild"

## Note: if you do not have samples with missing data for the Hierarch_variables, you can use
##   the function "SCOBI_deux_fast" instead of "SCOBI_dux" which may speed things up

