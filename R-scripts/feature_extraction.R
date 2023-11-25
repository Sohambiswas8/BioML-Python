# This is the R-script for feature selection in plant Resistance (R) Protein

# Available positive and negative data had been divided into Training and Test Data
# # with 80% and 20% ration in each data set

# Test and Train data file access for Positive Data

library(BiocManager)
library(protr)

### Load FASTA file
### system.file is for accessing example file in 'protr' package
library(seqinr)
positive_train <- read.fasta(file = 'pos_train.fasta')
positive_test <- read.fasta(file = 'pos_test.fasta')
length(positive_train) # No. of training sequence
length(positive_test)  # No. of test sequence

#............................................................................................
#............................................................................................

#########
## UPDATE: 
## First feature extraction from combined Positive file
## Secondly, same way feature extraction from the combined negative file

### Afterwards, divide them into training and test data


#............................................................................................
# POSITIVE DATASET

# Load fasta file > "combined_positive_dataset.fasta"

library(seqinr)
pos_data <- read.fasta(file = "combined_positive_dataset.fasta")
length(pos_data) # check the no. of protein sequences -> 98

# checking unusual Amino acid composition

pos_data_check <- pos_data[(sapply(pos_data, protcheck))]
length(pos_data_check)

# prot_data_check yield no sequences as having unusual Amino acid sequence.
## Hence, prot_data have 98 'perfect' sequences to go with !!

# DESCRIPTORS:
# 1. Amino Acid Composition (AAC) #extractAAC() function

positive_data <- readFASTA("combined_positive_dataset.fasta")
aac_pos_check <- extractAAC(positive_data[[1]])
aac_pos_check

# Since positive_data is a list, the variable can't pass the extractAAC() function

aac_pos <- t(sapply(positive_data, extractAAC))
aac_pos

# Dipeptide with extractDC()

dc_pos_check <- extractDC(postive_data[[1]])
dc_pos_check

dc_pos <- t(sapply(positive_data, extractDC))
head(dc_pos)


# Tripeptide with extractTC()

tc_pos <- t(sapply(positive_data, extractTC))
head(tc_pos, n = 36L)
head(tc_pos)
tc_pos[[1]]

#............................................................................................
## Autocorrelation descriptors

# -> MORAEU-BOROTO CORRELATION DESCRIPTOR

mbc_pos <- t(sapply(positive_data, extractMoreauBroto))
head(mbc_pos)

# -> Moran Correlation

mc_pos <- t(sapply(positive_data, extractMoran))
head(mc_pos, n = 36L)

# -> Geary Autocorrelation

geary_pos <- t(sapply(positive_data, extractGeary))
head(geary_pos, n = 36L)


##..........................................................................................
#CTD - Composition

c_pos <- t(sapply(positive_data, extractCTDC))
head(c_pos)

##CTD - Transition
t_pos <- t(sapply(positive_data, extractCTDT))
head(t_pos)

## CTD - Distribution
d_pos <- t(sapply(positive_data, extractCTDD))
head(d_pos)

##..................................................................................
# Conjoint-triad

conjoint_triad_pos <- t(sapply(positive_data, extractCTriad))
head(conjoint_triad_pos, n = 65L)

##..................................................................................
## Sequence Order Coupling Number

socp_pos <- t(sapply(positive_data, extractSOCN))
head(socp_pos)


##..................................................................................
#Quasi-sequence order

qso_pos <- t(sapply(positive_data, extractQSO))
head(qso_pos)

#...................................................................................
## Pseudo-amino acid composition

paac_pos <- t(sapply(positive_data, extractPAAC))
head(paac_pos)

#...................................................................................
## Amphi-philic PAAC
apaac_pos <- t(sapply(positive_data, extractAPAAC))
head(apaac_pos)

#...................................................................................
## Profile-based descriptors

# PSSM-based descriptors - extractPSSM(), extractPSSMAcc() & extractPSSMFeature

?extractPSSM()
class(positive_data)
# NO feature extraction for PSSM

#...................................................................................

## BLOSUM Matrix descriptor: extractBLOSUM()

?extractBLOSUM()
blosum_pos <- t(sapply(positive_data, extractBLOSUM, k = 5, lag = 7))
head(blosum_pos)
length(blosum_pos)


#...................................................................................
## scale-based descriptors
?extractDescScales()
descscales_pos <- t(sapply(positive_data, extractDescScales, 
                           propmat = "AATopo", pc = 20 
                           , lag = 7))
head(descscales_pos)
descscales_aamoe2d_pos <- t(sapply(positive_data, extractDescScales,
                                   propmat = "AAMOE2D", pc = 20,
                                   lag = 7))
descscales_aamoe3d_pos <- t(sapply(positive_data, extractDescScales, 
                                   propmat = "AAMOE3D", pc = 20,
                                   lag = 7))
View(descscales_aamoe2d_pos)
View(descscales_aamoe3d_pos)

descscales_aamolprop_pos <- t(sapply(positive_data, extractDescScales,
                                     propmat = "AAMolProp", pc = 20,
                                     lag = 7))
length(descscales_aamolprop_pos)
View(descscales_aamolprop_pos)

print(paste0("row: ",nrow(descscales_aamolprop_pos), 
             ", column: ", ncol(descscales_aamolprop_pos)))

print(paste0("row: ",nrow(descscales_aamoe3d_pos), 
             ", column: ", ncol(descscales_aamoe3d_pos)))

print(paste0("row: ",nrow(descscales_aamoe2d_pos), 
             ", column: ", ncol(descscales_aamoe2d_pos)))
## We calculated 4 different descriptors with scaled-based properties. MOE2D,
## MOE3D, Topology and MolProp
## combine all the 4 matrices [ by using reduce() func.]

# merge all matrices
scale_based_des_pos <- Reduce(function(x, y) cbind(x, y), 
                       list(descscales_aamolprop_pos, 
                            descscales_aamoe2d_pos, 
                            descscales_aamoe3d_pos, 
                            descscales_pos))

# check dimensions [should be 98 * 9408]
dim(scale_based_des_pos)
head(scale_based_des_pos, n = 60L)




































