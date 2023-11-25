getwd()
list.files()
install.packages('readxl')
library(readxl)
data_label <- read_excel("Merge_file.xlsx")
head(data_label)
tail(data_label)

##########################################
# feature table
data_feature <- read_excel("all_data.xlsx")
head(data_feature)
