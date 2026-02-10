#!/usr/bin/env Rscript
# SPDX-FileCopyrightText: 2025 Ovidio Garcia-Oliva
# SPDX-License-Identifier: CC-BY-4.0
# SPDX-FileContributor Ovidio Garcia-Oliva <ovidio.garcia@hereon.de
#
# This script consolidate and fix data from original files
# DO NOT EXECUTE: THIS IS INCLUDED JUST AS REFERENCE

library(lubridate)

CWG.nox = read.csv("/home/data/CWG_NOX.csv")
CWG.nox$date = dmy(CWG.nox$date)
CWG.nox$nitrate_plus_nitrite = as.numeric(CWG.nox$NOX)
CWG.nox = subset(CWG.nox,select=c(date,depth,nitrate_plus_nitrite))
CWG.nox$nitrate_plus_nitrite[CWG.nox$nitrate_plus_nitrite<2] = 1

CWG.nh4 = read.csv("/home/data/CWG_NH4.csv")
CWG.nh4$date = dmy(CWG.nh4$date)
CWG.nh4$ammonia = as.numeric(CWG.nh4$NH4)
CWG.nh4 = subset(CWG.nh4,select=c(date,depth,ammonia))
CWG.nh4$ammonia[CWG.nh4$ammonia<3] = 1.5 

CWG.po4 = read.csv("/home/data/CWG_PO4.csv")
CWG.po4$date = dmy(CWG.po4$date)
CWG.po4$`ortho-phosphate` = as.numeric(CWG.po4$PO4)
CWG.po4 = subset(CWG.po4,select=c(date,depth,`ortho-phosphate`))
CWG.po4$`ortho-phosphate`[CWG.po4$`ortho-phosphate`<1] = 0.5 

CWG.tp = read.csv("/home/data/CWG_TP.csv")
CWG.tp$date = mdy(CWG.tp$date)
CWG.tp$total_phosphorus = as.numeric(CWG.tp$TP)
CWG.tp = subset(CWG.tp,select=c(date,depth,total_phosphorus))
CWG.tp = na.omit(CWG.tp)
CWG.tp$total_phosphorus[CWG.tp$total_phosphorus<1] = 0.5 


CWG.do = read.csv("/home/data/Full_WG_DO.csv")
CWG.do$date = dmy(CWG.do$Fecha)
CWG.do$dissolved_oxygen = as.numeric(CWG.do$DO)
CWG.do = subset(CWG.do,select=c(date,depth,dissolved_oxygen))
CWG.do = na.omit(CWG.do)
CWG.do$dissolved_oxygen[CWG.do$dissolved_oxygen<0] = 0. 

CWG.temp = read.csv("/home/data/Full_WG_Temp.csv")
CWG.temp$date = dmy(CWG.temp$Fecha)
CWG.temp$temperature = as.numeric(CWG.temp$Temp)
CWG.temp = subset(CWG.temp,select=c(date,depth,temperature))

files = list(CWG.temp,
             CWG.do,
             CWG.nh4,
             CWG.nox,
             CWG.po4,
             CWG.tp
              )

merge.recursively = function(db,...){
  is.first = 1
  for(x in db){
    if(is.first){
      db = x 
      is.first = 0
    }
    else{
      db = merge(db,x,all=T)  
    } 
  }
  return(db)
}

for(x in files) print(colnames(x))

db = merge.recursively(files)

db$total_phosphorus[db$`ortho-phosphate`>db$total_phosphorus] = NA

colnames(db)

plot(db)
db[,2:8] = round(db[,2:8],digits=3)
write.csv(db,'Lake_Atitlan_water_quality_CWG_2010-2024.csv',row.names = F)


