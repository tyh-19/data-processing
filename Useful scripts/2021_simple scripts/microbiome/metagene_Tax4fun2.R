# install
install.packages(pkgs = "C:/Users/Tao/Desktop/Tax4Fun2_1.1.5.tar.gz", repos = NULL, source = TRUE)
library(Tax4Fun2)
setwd("C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/Tax4Fun2/")
??Tax4Fun2

buildReferenceData(path_to_working_directory="C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/Tax4Fun2/reference",use_force=T)
testReferenceData(path_to_reference_data = "C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/Tax4Fun2/reference/Tax4Fun2_ReferenceData_v2")

makeFunctionalPrediction(path_to_otu_table = 'C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/counts-G.txt', path_to_reference_data = 'C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/Tax4Fun2/reference/Tax4Fun2_ReferenceData_v2', 
                         path_to_temp_folder = 'C:/Users/Tao/Desktop/Lulab/2020/02.Project/02. Data/metagene/Tax4Fun2/Ref100NR', database_mode = 'Ref100NR', 
                         normalize_by_copy_number = TRUE, min_identity_to_reference = 0.97, normalize_pathways = FALSE)


makeFunctionalPrediction()