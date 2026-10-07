#--- Load libraries
library(reticulate)

# To create a Python environment containing 'pycox': change to TRUE if needed
if (FALSE) {
  
  cnd <- conda_binary("auto")
  
  if ("pycox-clean" %in% conda_list(conda = cnd)$name)
    conda_remove("pycox-clean", conda = cnd)
  
  conda_create(
    "pycox-clean",
    c("python=3.10", "pip", "numpy=1.23.5", "pandas=1.4.4"),#"zlib",
    channel = "conda-forge",
    conda = cnd
  )
  
  conda_install(
    "pycox-clean",
    c("pycox==0.3.0", "torchtuples==0.2.2", "torch==2.0.1"),
    pip = TRUE,
    conda = cnd
  )

}

reticulate::use_condaenv("pycox-clean", required = TRUE)
#py_config()
reticulate::import("pycox", convert = TRUE) 
library(survival)
library(truncnorm)
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)
library(randomForestSRC)
library(xgboost)
library(survXgboost)
library(survivalmodels)
library(pec)
set_seed(47)

 
