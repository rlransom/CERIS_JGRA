###### Two sections: the 1st section is to identify the environmental index (by Xianran Li lixr@iastate.edu); 
######               the 2nd section is for performance prediction combined with genomic prediction(by Tingting Guo, tguo@iastate.edu) 
############################################################
##### Option (line 32-34: select the crop & trait
##### Line 95-97: modify the value based on the searching results

#################################################################################################
## install required packages
if (!require(dplyr)) { install.packages("dplyr", repos = "https://cloud.r-project.org");}
if (!require(corrgram)) { install.packages("corrgram", repos = "https://cloud.r-project.org");}
if (!require(data.table)) { install.packages("data.table", repos = "https://cloud.r-project.org");}
if (!require(colorspace)) { install.packages("colorspace", repos = "https://cloud.r-project.org");}
if (!require(RColorBrewer)) { install.packages("RColorBrewer", repos = "https://cloud.r-project.org");}

if (!require(rrBLUP)) { install.packages("rrBLUP", repos = "https://cloud.r-project.org");}
if (!require(BGLR)) { install.packages("BGLR", repos = "https://cloud.r-project.org");}
if (!require(yarrr)) { install.packages("yarrr", repos = "https://cloud.r-project.org");}
if (!require(openxlsx)) { install.packages("openxlsx", repos = "https://cloud.r-project.org");}

col_wdw <- 25;
col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1)
t_base <- 50; t_max1 <- 100; t_max2 <- 1000; Haun_threshold <- 0; p <- 1 

###
Top_dir <- 'C:/Users/schamart/Box/CSSA virtual meeting/Workshops/Phenotype plasticity/Workshop_2020/'

###### If you modify some functions in this file, please run Line 27 each time to reload the updated functions
subfunction_file <- paste(Top_dir, 'Sub_functions_Workshop.r', sep = '');
source(subfunction_file);

######################################
experiment <- 'Maize'; ## Options: Maize, Wheat, Oat; 
##### For Maize, the traits are DTA, PH, and Yield; for Sorghum: the trait is called as FTgdd; for Rice, the trait is FTdap; ; 
trait <- 'DTA'; ##### Options: for Maize: FTdap, PH,; Wheat, FTdap, PH, or GY; for oat: FTdap, PH, GY
######################################

exp_dir <- paste(Top_dir, experiment, '/', sep = '')
env_meta_file <- paste(exp_dir, 'Env_meta_table.txt', sep = ''); ## make sure the PlantingDate is formated as 'YYYY-MM-DD'
env_meta_info_0 <- read.table(env_meta_file, header = T, sep = "\t", stringsAsFactors = F);

searching_daps <- 100; #### For wheat, for oat if (experiment == '1Sorghum') { searching_daps <- 122};
ptt_ptr_file <- paste(exp_dir, nrow(env_meta_info_0), 'Envs_envParas_DAP',  searching_daps, '.txt', sep = ''); ##  
PTT_PTR <- read.table(ptt_ptr_file, header = T , sep = "\t");
Paras <- c('DL', 'GDD', 'PTT', 'PTR', 'PTS');

exp_traits_file <- paste(exp_dir, 'Trait_records.txt', sep = '');
exp_traits <- read.table(exp_traits_file, sep = "\t", header = T, stringsAsFactors = F, na.string = 'NA');

if(!('FTdap' %in% colnames(exp_traits))) {exp_traits$FTdap <- exp_traits$DTA};

all_env_codes <- unique(exp_traits$env_code);
env_cols <- rainbow_hcl(length(all_env_codes), c = 80, l = 60, start = 0, end = 300, fixup = TRUE, alpha = 0.75);

##########################################################
  lInd <- which(colnames(exp_traits) == 'line_code'); eInd <- which(colnames(exp_traits) == 'env_code'); tInd <- which(colnames(exp_traits) == trait);
  exp_trait_dir <- paste(exp_dir, trait,  '/',  sep = ''); if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir, recursive= T)};
  exp_trait <- exp_traits[,c(lInd, eInd, tInd)]; 
  
  colnames(exp_trait)[3] <- 'Yobs';
  exp_trait <- aggregate(Yobs ~  line_code + env_code, data = exp_trait, mean) ## To make sure only one phenotype record per each line in each environment
  exp_trait <- exp_trait[!is.na(exp_trait$Yobs),];

  line_codes <- unique(exp_trait$line_code); 
  env_mean_trait_0 <- na.omit(aggregate(x = exp_trait$Yobs, by = list(env_code = exp_trait$env_code), mean, na.rm = T));
  colnames(env_mean_trait_0)[2] <- 'meanY';
  env_mean_trait <- env_mean_trait_0[order(env_mean_trait_0$meanY),];

### pairwise correlations among environments; trait distribution across environments;
### two figures and the correspondent output files will be saved in the trait directory;
  try(Pairwise_trait_env_distribution_plot(exp_trait, exp_trait_dir, trait, all_env_codes, env_meta_info_0));

  ##### searching the critical window having the strongest correlation with environmental mean
  ##### the window can be adjusted based on biological meaning
  ##### 'FTgdd_7Envs_PTTPTR_0LOO_cor.txt' stores all correlations from all the tested windows and environmental parameters;
  ##### 'MaxR_FTgdd_7Envs_0LOO.png' is the visualization 
  pop_cor_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_', 0, 'LOO_cor.txt', sep = '');
  Exhaustive_search(env_mean_trait, PTT_PTR, searching_daps, exp_trait_dir, exp_traits$FTdap, trait, 1, searching_daps, searching_daps, 0, Paras, pop_cor_file)#; searching_daps, searching_daps);

#######################################################
##### Contributed from Laura Cortes (ltibbs@iastate.edu) ##########
##### To help navigate the strongest correlations identified for each environment parameter. 
##### View the results with 'FTgdd_7Envs_PTTPTR_0LOO_cor.txt' and 'MaxR_FTgdd_7Envs_0LOO.png'
##### R_PTT means the original correlation between PTT and population mean; while nR_PTT means the correlation of 0 - R_PTT. 
  search.results <- read.table(pop_cor_file, header = T)
  search.results <- search.results %>%
    tidyr::gather(key="Parameter", value="Corr", -Day_x, -Day_y, -window, -pop_code) %>%
    arrange(-Corr)
  search.results <- search.results %>%
    group_by(Parameter) %>%
    top_n(5, Corr)    
  View(search.results);
  
##########################################################
### Change the following three parameters for the window and environmental parameter with the strongest correlation
  maxR_dap1 <- 23;
  maxR_dap2 <- 43;
  kPara_Name <- 'GDD';
##### #################################################### 

  PTT_PTR_ind <-  which(colnames(PTT_PTR) == kPara_Name); 
#### Visualization of the relationships between environmental mean and environmental parameters from the selected window.  
  Plot_Trait_mean_envParas(env_mean_trait, PTT_PTR, maxR_dap1, maxR_dap2, trait, exp_trait_dir, env_cols, Paras);  

#### Output intercept and slope estimation for each line based on environmental mean and environmental parameter
  Slope_Intercept(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, exp_trait, line_codes, exp_trait_dir);
  #### LOOCV function for 1 -> 2
  obs_prd_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_', kPara_Name, 'D', maxR_dap1, '_', maxR_dap2, '.txt', sep = '');
  LOO_pdf_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Env_LOO_by_Lines_', kPara_Name, 'D', maxR_dap1, '_', maxR_dap2, '.png', sep = '');
  if (!file.exists(obs_prd_file)) { 
    prdM <- LOOCV(maxR_dap1, maxR_dap2, env_mean_trait, PTT_PTR, PTT_PTR_ind, exp_trait, obs_prd_file, p)
  }
  Plot_prediction_result(obs_prd_file, all_env_code, prdM, kPara_Name, LOO_pdf_file,env_cols);


####################################################################################################################
############ From Tingting Guo (tguo@iastate.edu) #################
#### 3 prediction scenarios: 1->2; 1->3; 1->4
LOC=read.table(paste(exp_dir,"Env_meta_table.txt",sep=""),header=T,sep="\t");
geno=read.table(paste(exp_dir,"Genotype.txt",sep=""),header=T,sep="\t");
pheno=read.table(paste(exp_trait_dir,"LbE_table.txt",sep=""),header=F,sep="\t");
envir=read.table(paste(exp_trait_dir,trait,'_envMeanPara_', maxR_dap1, '_', maxR_dap2, '.txt',sep=""),header=T,sep="\t");

tt.line=nrow(pheno)*0.5; ##remove the environment for which the number of missing line > tt.line 
tt.e=c(ncol(pheno)-1)-3; ##remove the line for which the number of missing environment > tt.e
enp=which(colnames(envir) == kPara_Name); 
fold=10;
reshuffle=2;

###Through reaction norm parameter
out1.2=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.E",fold,reshuffle)  ## 1->2 prediction
write.xlsx(out1.2, file = paste(exp_trait_dir,"Result_Norm_1.2.xlsx",sep=""))

out1.3=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.G",fold,reshuffle)  ## 1->3 prediction
write.xlsx(out1.3, file = paste(exp_trait_dir,"Result_Norm_1.3.xlsx",sep=""))

out1.4=JGRA(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.GE",fold,reshuffle) ## 1->4 prediction
write.xlsx(out1.4, file = paste(exp_trait_dir,"Result_Norm_1.4.xlsx",sep=""))

### Through marker effect
out1.2=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.E",fold,reshuffle)
write.xlsx(out1.2, file = paste(exp_trait_dir,"Result_Marker_1.2.xlsx",sep=""))

out1.3=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.G",fold,reshuffle)
write.xlsx(out1.3, file = paste(exp_trait_dir,"Result_Marker_1.3.xlsx",sep=""))

out1.4=JGRA.marker(LOC,pheno,geno,envir,enp,tt.line,tt.e,mets="RM.GE",fold,reshuffle) #### Need a longer time to run than others; be patient
write.xlsx(out1.4, file = paste(exp_trait_dir,"Result_Marker_1.4.xlsx",sep=""))


