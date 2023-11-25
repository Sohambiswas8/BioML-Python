library(caret)

# assume pos_data and neg_data are your matrix variables for positive and negative data

# create labels for your data
pos_labels <- rep("positive", nrow(pos_data))
neg_labels <- rep("negative", nrow(neg_data))

# combine data and labels
all_data <- rbind(pos_data, neg_data)
all_labels <- c(pos_labels, neg_labels)

# set the seed for reproducibility
set.seed(123)

# split data into 80% training and 20% testing
train_index <- createDataPartition(all_labels, p = 0.8, list = FALSE)
train_data <- all_data[train_index, ]
train_labels <- all_labels[train_index]
test_data <- all_data[-train_index, ]
test_labels <- all_labels[-train_index]
