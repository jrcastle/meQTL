library(MatrixEQTL)
require(GEM)
setwd('/home/jcastle/VitaminUse_Study/eQTL')
DATADIR <- "/mnt/DATA/meqtl/seq/"


##### RUN cis-meQTL ANALYSIS USING GEM #####
covariate_file = "../data/VDUse_vs_methylation/cov_VDUse_EUR.txt"
methylation_file = "../data/VDUse_vs_methylation/meth_VDUse_EUR.txt"
env_file = '../data/VDUse_vs_methylation/env_VDUse_EUR.txt'
Emodel_pv = 1
Emodel_result_file_name = "Result_Emodel_EUR.txt"
Emodel_qqplot_file_name = "QQplot_Emodel_EUR.jpg"
GEM_Emodel(env_file, 
           covariate_file, 
           methylation_file, 
           Emodel_pv, 
           Emodel_result_file_name, 
           Emodel_qqplot_file_name, 
           savePlot=TRUE
           )

