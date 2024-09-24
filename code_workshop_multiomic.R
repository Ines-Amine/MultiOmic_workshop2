library(tidyverse)
library(RGCCA)
devtools::install_github(repo="https://github.com/rgcca-factory/RGCCA.git")

ncore.use <- parallel::detectCores() - 1
set.seed(2000)

### Step 0: prepare the dataset #####

# Omic data
proteomic_data <- readRDS("2_Data/2_final_data/proteomic_adjusted.Rds") |>
  arrange(HelixID) %>% select(-HelixID)
celltype_data  <- readRDS("2_Data/2_final_data/cell_type_adjusted.Rds")|>
  arrange(HelixID) %>% select(-HelixID)
methylome_data  <- readRDS("2_Data/2_final_data/methylome_reduced_adjusted_celltypes_killed.Rds") %>%
  select(-HelixID)
methylome_data <- methylome_data[,1:100]

# Outcome
Y <- readRDS("2_Data/2_final_data/health_scores_adjusted.Rds")%>% 
  arrange(HelixID) %>% dplyr::select(mean_health_score) |> as.matrix()

# Divide it in train test
id_train <- sample(c(TRUE, FALSE), nrow(Y), replace=TRUE, prob=c(0.7,0.3))

X_train <- list(proteins = proteomic_data[id_train,],
          celltype=celltype_data[id_train,],
          dnam = methylome_data[id_train,],
          Y=Y[id_train,])


X_test <- list(proteins = proteomic_data[!id_train,],
               celltype=celltype_data[!id_train,],
               dnam = methylome_data[!id_train,],
               Y=Y[!id_train,])


# Defining the connection matrix
connection = 1 - diag(4)


### Step 1: tuning parameters using cross validation for number of components  #

# Number of components  
print("- Tuning of number of components per blocks")
cv_sgcca_ncomp  <- rgcca_cv(
  blocks = X_train,
  response = 4,
  method = "rgcca",
  par_value = expand.grid(proteins = 1:3, 
                          celltype =1:3,  
                          dnam = 1:3, Y = 1),
  par_type = "ncomp",
  n_cores = ncore.use,
  connection = connection,
  n_run = 1
)
plot(cv_sgcca_ncomp )
ncomp = cv_sgcca_ncomp$best_params


## Step 2: RGCCA with optimized parameters ##
rgcca_res <- rgcca(cv_sgcca_ncomp)
summary(rgcca_res)

## Step 3: Interpretation of the model and performance ##

# Significance: bootstrap (Time estimated: 2 min)
bootstrap <- rgcca_bootstrap(rgcca_res,n_cores = ncore.use)
View(bootstrap$stats)

# Plots loading for all
plot(bootstrap,comp=1,n_mark=20,display_order=T,type="weight")

# Plots loading for specific components
plot(bootstrap,comp=1,block=1,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=2,block=1,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=1,block=2,n_mark=30,display_order=T,type="weight")
plot(bootstrap,comp=1,block=3,n_mark=30,display_order=T,type="weight")


### Step 4: Prediction ability on test set 

# With rgcca_predict, you can get standard indicators (RMSE, MAE, R2) using different predictor models
rgcca_predict_res <- rgcca_predict(rgcca_res,blocks_test=X_test,prediction_model="lm")
rgcca_predict_res$metric$test

# Projection of the variables on the latent components for the test set
latent_variables <- rgcca_predict_res$projection %>% purrr::reduce(cbind)%>% as.data.frame()
name_components <- c("Prot_LC1","Prot_LC2","Prot_LC3",
                               "WBC_LC1","WBC_LC2","WBC_LC3",
                               "DNAm_LC1")
colnames(latent_variables)<- name_components
  
# Estimate R2 of each latent variable by runing a linear model with only this component
for (comp in name_components){
  data_r2 <- cbind(latent_variables,X_test$Y) %>% as.data.frame()
  form = paste0(" X_test$Y ~", comp)
  res <- summary(lm(as.formula(form),data=data_r2))
  r2 <- res$r.squared
  print(paste0("R2 of ",comp," = ", 100*signif(r2,2),"%"))
}
  
# Estimate correlation between components
cor(latent_variables) %>% round(2)


  


  
