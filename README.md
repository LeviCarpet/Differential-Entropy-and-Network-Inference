# Differential-Entropy-and-Network-Inference
Contains the necessary R codes associated to the article: "Identification of Nanog-Helper Genes in Early Mouse Embryo Differentiation Using Differential Entropy and Network Inference "

# Data-preprocessing 
Contains the R code used to pre-process the 5 single-cell datasets, including the mRNA transformations for the scRT-qPCR datasets.

# Standard Analysis
Contains the R code used for the correlation and PCA analysis of the 5 single-cell datasets.

# Entropy Generation Bootstrap
Contains the R code used for the bootstrapped computation of the differential entropy, and includes the verification of the fittings using the KS statistic.

# CARDAMOM 
Contains the R code to create the necessary datafiles that are used in CARDAMOM. Once the files are created, run CARDAMOM on python as described in Ventre et al. 2023, and then use the second part of this code to create the network graphs as in the manuscript.
