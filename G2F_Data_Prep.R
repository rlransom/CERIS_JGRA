###################
#2017 G2f Data Prep

#Clear Environment
remove(list = ls())

##Load Data
#Hybrid PCs
Hybrid_PCs <- read.csv('./G2F/Data_Prep/PCs_from_hybrid_D.csv')
Hybrids <- Hybrid_PCs[,c(1,2)]

#Trait Records (Yield)
Phenotype <- read.csv('./G2F/Data_Prep/Hybrid_Phenotype_2017.csv')
Phenotype <- Phenotype[, c(1, 2, 5, 20, 21, 34)]

#Join
Trait_records <- inner_join(Hybrids, Phenotype, by = 'Pedigree')
Trait_records <- Trait_records[,c(1,4,7)]

#Field Metadata
Meta <- read.csv('./G2F/Data_Prep/g2f_2017_field_metadata.csv')

#Join
Env_meta_table <- left_join(Meta, Phenotype, by = 'Location')
Env_meta_table <- Env_meta_table[,c(1,2,3,6,7)]
Env_meta_table_w_ped <- unique(Env_meta_table)
Env_meta_table <- Env_meta_table[,c(1:3,5)]
Env_meta_table <- unique(Env_meta_table)

Env_meta_table$Date <- strptime(as.character(Env_meta_table$Date_Planted), "%m/%d/%Y")

Env_meta_table <- Env_meta_table[,-4]


#Environmental Parameters
Environ_Params <- read.csv('./G2F/Data_Prep/g2f_2-17_Weather_Data.csv')

Environ_Params$Date <- as.Date(with(Environ_Params,paste(Year,Month,Day,sep="-")),"%Y-%m-%d")
Environ_Params <- Environ_Params[,-c(2:5)]

#Environ_Params <- aggregate(cbind(Environ_Params$Temperature_C, Environ_Params$Dew_Point_C, Environ_Params$Relative_Humidity), by=list(Date=Environ_Params$Date), FUN= mean)

library(dplyr)
Environ_Params %>%
  group_by(Location, Date) %>% 
  summarise_each(funs(mean))

Environ_Params <- Environ_Params %>% 
  group_by(Location, Date) %>%
  summarise(across(c("Temperature_C","Dew_Point_C","Relative_Humidity"), mean))


names(Environ_Params)[c(1,3,4,5)] <- c("env_code","Temp","DewPoint","Humidity")

names(Env_meta_table)[c(1,2,3,4)]<- c("env_code","lat","lon","PlantingDate")

names(Trait_records)[c(1,2)]<- c("line_code","env_code")


#Export Files

write.csv(Trait_records, './G2F/Trait_records.csv', row.names = FALSE)
write.csv(Env_meta_table, './G2F/Env_meta_table.csv', row.names = FALSE)
write.csv(Environ_Params, './G2F/Environ_Params.csv', row.names = FALSE)

