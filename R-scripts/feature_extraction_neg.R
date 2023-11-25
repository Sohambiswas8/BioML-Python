###
###
###
# THIS IS THE NEGATIVE DATASET

library(BiocManager)
library(protr)

### Load FASTA file
### system.file is for accessing example file in 'protr' package
library(seqinr)
negative_train <- read.fasta(file = 'neg_train.fasta')
negative_test <- read.fasta(file = 'neg_test.fasta')
length(negative_train) # No. of training sequence
length(negative_test)  # No. of test sequence

#............................................................................................
#............................................................................................

#########
## UPDATE: 

## Same way feature extraction from the combined negative file

### Afterwards, divide them into training and test data


#............................................................................................
# NEGATIVE DATASET

# Load fasta file > "combined_negative_dataset.fasta"

library(seqinr)
neg_data <- read.fasta(file = "combined_negative_dataset.fasta")
length(neg_data) # check the no. of protein sequences -> 198

View(neg_data)

# checking unusual Amino acid composition

neg_data_check <- neg_data[(sapply(neg_data, protcheck))]
length(neg_data_check)

# neg_data_check yield no sequences as having unusual Amino acid sequence.
## Hence, neg_data have 198 'perfect' sequences to go with !!

# DESCRIPTORS:
# 1. Amino Acid Composition (AAC) #extractAAC() function

negative_data <- readFASTA("combined_negative_dataset.fasta")
aac_neg_check <- extractAAC(negative_data[[1]])
aac_neg_check

# Since positive_data is a list, the variable can't pass the extractAAC() function

aac_neg <- t(sapply(negative_data, extractAAC))
aac_neg

# Dipeptide with extractDC()

#dc_neg_check <- extractDC(negative_data[[1]])
#dc_neg_check

dc_neg <- t(sapply(negative_data, extractDC))
head(dc_neg)


# Tripeptide with extractTC()

tc_neg <- t(sapply(negative_data, extractTC))
head(tc_neg, n = 36L)
head(tc_neg)
tc_neg[[1]]

#............................................................................................
## Autocorrelation descriptors

# -> MORAEU-BOROTO CORRELATION DESCRIPTOR

mbc_neg <- t(sapply(negative_data, extractMoreauBroto))
head(mbc_neg)

# -> Moran Correlation

mc_neg <- t(sapply(negative_data, extractMoran))
head(mc_neg, n = 36L)

# -> Geary Autocorrelation

geary_neg <- t(sapply(negative_data, extractGeary))
head(geary_neg, n = 36L)


##..........................................................................................
#CTD - Composition

c_neg <- t(sapply(negative_data, extractCTDC))
head(c_neg)

##CTD - Transition
t_neg <- t(sapply(negative_data, extractCTDT))
head(t_neg)

## CTD - Distribution
d_neg <- t(sapply(negative_data, extractCTDD))
head(d_neg)

##..................................................................................
# Conjoint-triad

conjoint_triad_neg <- t(sapply(negative_data, extractCTriad))
head(conjoint_triad_neg, n = 65L)

##..................................................................................
## Sequence Order Coupling Number

socp_neg <- t(sapply(negative_data, extractSOCN))
head(socp_neg)


##..................................................................................
#Quasi-sequence order

qso_neg <- t(sapply(negative_data, extractQSO))
head(qso_neg)

#...................................................................................
## Pseudo-amino acid composition

paac_neg <- t(sapply(negative_data, extractPAAC))
head(paac_neg)

#...................................................................................
## Amphi-philic PAAC
apaac_neg <- t(sapply(negative_data, extractAPAAC))
head(apaac_neg)

#...................................................................................
## Profile-based descriptors

# PSSM-based descriptors - extractPSSM(), extractPSSMAcc() & extractPSSMFeature

?extractPSSM()
class(positive_data)
# NO feature extraction for PSSM

#...................................................................................

## BLOSUM Matrix descriptor: extractBLOSUM()

?extractBLOSUM()
blosum_neg <- t(sapply(negative_data, extractBLOSUM, k = 5, lag = 7))
head(blosum_neg)
length(blosum_neg)


#...................................................................................
## scale-based descriptors
?extractDescScales()
descscales_neg <- t(sapply(negative_data, extractDescScales, 
                           propmat = "AATopo", pc = 20 
                           , lag = 7))
head(descscales_neg)
descscales_aamoe2d_neg <- t(sapply(negative_data, extractDescScales,
                                   propmat = "AAMOE2D", pc = 20,
                                   lag = 7))
descscales_aamoe3d_neg <- t(sapply(negative_data, extractDescScales, 
                                   propmat = "AAMOE3D", pc = 20,
                                   lag = 7))
View(descscales_aamoe2d_neg)
View(descscales_aamoe3d_neg)

descscales_aamolprop_neg <- t(sapply(negative_data, extractDescScales,
                                     propmat = "AAMolProp", pc = 20,
                                     lag = 7))
length(descscales_aamolprop_neg)
View(descscales_aamolprop_neg)

print(paste0("row: ",nrow(descscales_aamolprop_neg), 
             ", column: ", ncol(descscales_aamolprop_neg)))

print(paste0("row: ",nrow(descscales_aamoe3d_neg), 
             ", column: ", ncol(descscales_aamoe3d_neg)))

print(paste0("row: ",nrow(descscales_aamoe2d_neg), 
             ", column: ", ncol(descscales_aamoe2d_neg)))

print(paste0("row: ",nrow(descscales_neg), 
             ", column: ", ncol(descscales_neg)))

## We calculated 4 different descriptors with scaled-based properties. MOE2D,
## MOE3D, Topology and MolProp
## combine all the 4 matrices [ by using reduce() func.]

# merge all matrices
scale_based_des_neg <- Reduce(function(x, y) cbind(x, y), 
                              list(descscales_aamolprop_neg, 
                                   descscales_aamoe2d_neg, 
                                   descscales_aamoe3d_neg, 
                                   descscales_neg))

# check dimensions [should be 198 * 9408]
dim(scale_based_des_neg)
head(scale_based_des_neg, n = 60L)




































