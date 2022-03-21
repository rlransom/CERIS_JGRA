rm(list = ls())


## install required packages
if (!require(dplyr)) { install.packages("dplyr", repos = "https://cloud.r-project.org");}
if (!require(corrgram)) { install.packages("corrgram", repos = "https://cloud.r-project.org");}
if (!require(data.table)) { install.packages("data.table", repos = "https://cloud.r-project.org");}
if (!require(colorspace)) { install.packages("colorspace", repos = "https://cloud.r-project.org");}
if (!require(RColorBrewer)) { install.packages("RColorBrewer", repos = "https://cloud.r-project.org");}

if (!require(rrBLUP)) { install.packages("rrBLUP", repos = "https://cloud.r-project.org");}
if (!require(BGLR)) { install.packages("BGLR", repos = "https://cloud.r-project.org");}
if (!require(yarrr)) { install.packages("yarrr", repos = "https://cloud.r-project.org");}
#if (!require(openxlsx)) { install.packages("openxlsx", repos = "https://cloud.r-project.org");}
if (!require(lme4)) { install.packages("lme4", repos = "https://cloud.r-project.org");}
if (!require(emmeans)) { install.packages("emmeans", repos = "https://cloud.r-project.org");}



col_wdw <- 25;
col_palette <- diverge_hcl(col_wdw + 1, h = c(260, 0), c = 100, l = c(50, 90), power = 1)
t_base <- 50; t_max1 <- 100; t_max2 <- 1000; Haun_threshold <- 0; p <- 1 

###
#Top_dir <- 'D:/0GbE/CERES_MS/MP/GitHub_files/'
Top_dir <- '~/School/Graduate/Spring 2022/Envirotyping/CERIS_JGRA/'

###### If you modify some functions in this file, please run Line 27 each time to reload the updated functions
subfunction_file <- paste(Top_dir, 'G2F_Sub_functions.r', sep = '');
source(subfunction_file);

######################################
experiment <- 'G2F'; ## Options: Maize, Wheat, Oat; 
##### For Maize, the traits are DTA, PH, and Yield; for Sorghum: the trait is called as FTgdd; for Rice, the trait is FTdap; ; 
trait <- 'Yield'; ##### Options: for Maize: FTdap, PH,; Wheat, FTdap, PH, or GY; for oat: FTdap, PH, GY
######################################

exp_dir <- paste(Top_dir, experiment, '/', sep = '')
env_meta_file <- paste(exp_dir, 'Env_meta_table.csv', sep = ''); ## make sure the PlantingDate is formated as 'YYYY-MM-DD'
env_meta_info_0 <- read.table(env_meta_file, header = T, sep = ",", stringsAsFactors = F);

searching_daps <- 100; #### For wheat, for oat if (experiment == '1Sorghum') { searching_daps <- 122};
#ptt_ptr_file <- paste(exp_dir, nrow(env_meta_info_0), 'Envs_envParas_DAP',  searching_daps, '.txt', sep = ''); ##  
ptt_ptr_file <- '~/School/Graduate/Spring 2022/Envirotyping/CERIS_JGRA/G2F/Environ_Params.csv'
#PTT_PTR <- read.table(ptt_ptr_file, header = T , sep = "\t");
PTT_PTR <- read.table(ptt_ptr_file, header = T , sep = ",");
#Paras <- c('DL', 'GDD', 'PTT', 'PTR', 'PTS');
Paras <- c('Temp', 'DewPoint', 'Humidity');

exp_traits_file <- paste(exp_dir, 'Trait_records.csv', sep = '');
exp_traits <- read.table(exp_traits_file, sep = ",", header = T, stringsAsFactors = F, na.string = 'NA');

all_env_codes <- unique(exp_traits$env_code);
env_cols <- rainbow_hcl(length(all_env_codes), c = 80, l = 60, start = 0, end = 300, fixup = TRUE, alpha = 0.75);

##########################################################
lInd <- which(colnames(exp_traits) == 'line_code'); eInd <- which(colnames(exp_traits) == 'env_code'); tInd <- which(colnames(exp_traits) == trait);
exp_trait_dir <- paste(exp_dir, trait,  '/',  sep = ''); if (!dir.exists(exp_trait_dir))  { dir.create(exp_trait_dir, recursive= T)};
exp_trait <- exp_traits[,c(lInd, eInd, tInd)]; 

colnames(exp_trait)[3] <- 'Yobs';
exp_trait <- aggregate(Yobs ~  line_code + env_code, data = exp_trait, mean) ## To make sure only one phenotype record per each line in each environment
exp_trait <- exp_trait[!is.na(exp_trait$Yobs),];

# ### remove outlier environments, such as one with high missing rate
# if (trait == 'FT') {
#   env_outliers <- c('03NY06', '08MO07');
#   exp_trait <- exp_trait[!(exp_trait$env_code %in% env_outliers),];
# }

line_codes <- unique(exp_trait$line_code); 
##  env_mean_trait_0 <- na.omit(aggregate(x = exp_trait$Yobs, by = list(env_code = exp_trait$env_code), mean, na.rm = T));
##  colnames(env_mean_trait_0)[2] <- 'meanY';

################# modified to add the opitions for estimating environmental mean by arithmetic mean (ari), mixed linear model (mlm), or least-square (emm)
env_mean_method <- 'mlm'; ### 'mlm', 'emm'. Modifying this based on the method to calculate environmental mean
exp_trait_m <- exp_trait
exp_trait_m$env_code <- as.factor(exp_trait_m$env_code)
exp_trait_m$line_code <- as.factor(exp_trait_m$line_code)

# if (env_mean_method == 'ari') {
#  env_mean_trait_0 <- na.omit(aggregate(x = exp_trait$Yobs, by = list(env_code = exp_trait$env_code), mean, na.rm = T));
#  colnames(env_mean_trait_0)[2] <- 'meanY';
# } else if (env_mean_method == 'mlm') {
#  lm_ <- lmer(Yobs ~ env_code + (1|line_code), data = exp_trait_m)
#  mlm_col <- c(fixef(lm_)[1], fixef(lm_)[1] + fixef(lm_)[-1]) ## BLUE for environnment
#  env_mean_trait_0 <- data.frame(env_code = as.vector(unique(exp_trait_m$env_code), meanY = mlm_col)
#  } else if (env_mean_method == 'emm') {
#   env_n <- length(as.vector(unique(exp_trait_m$env_code)))
#   trait_ori_lm1 <- lm(Yobs ~ env_code + line_code, data = exp_trait_m)
#   trait_ori.pred1 <- matrix(predict(ref_grid(trait_ori_lm1)), nrow = env_n)
#   emm_col <- apply(trait_ori.pred1, 1, mean) ### marginal mean for environments
#   env_mean_trait_0 <- data.frame(env_code = as.vector(unique(exp_trait_m$env_code), meanY = emm_col)
#  }

#########################################
#mlm
lm_ <- lmer(Yobs ~ env_code + (1|line_code), data = exp_trait_m)
mlm_col <- c(fixef(lm_)[1], fixef(lm_)[1] + fixef(lm_)[-1]) ## BLUE for environment
env_mean_trait_0 <- data.frame(env_code = as.vector(unique(exp_trait_m$env_code)), meanY = mlm_col)


#########################################


env_mean_trait <- env_mean_trait_0[order(env_mean_trait_0$meanY),];

### pairwise correlations among environments; trait distribution across environments;
### two figures and the correspondent output files will be saved in the trait directory;
try(Pairwise_trait_env_distribution_plot(exp_trait, exp_trait_dir, trait, all_env_codes, env_meta_info_0));

##### searching the critical window having the strongest correlation with environmental mean
##### the window can be adjusted based on biological meaning
##### 'FT_9Envs_PTTPTR_0LOO_cor.txt' stores all correlations from all the tested windows and environmental parameters;
##### 'MaxR_FTgdd_9Envs_0LOO.png' is the visualization 
pop_cor_file <- paste(exp_trait_dir, trait, '_', nrow(env_mean_trait), 'Envs_PTTPTR_', 0, 'LOO_cor.csv', sep = '');
Exhaustive_search(env_mean_trait, PTT_PTR, searching_daps, exp_trait_dir, exp_traits$Yield, trait, 1, 25, 50, 0, Paras, pop_cor_file)#; searching_daps, searching_daps);

##########################################################
### Change the following three parameters for the window and environmental parameter with the strongest correlation
maxR_dap1 <- 22;
maxR_dap2 <- 37;
kPara_Name <- 'Temp';
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


