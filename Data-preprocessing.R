######################### Pre-Processing File #########################
## In this file, we pre-process the all datasets into mRNA count matrices
## where rows are different cells, and columns are 'Code','Population' and Genes. 

######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### Setting the directory for loading data ######### ######### 
setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Entropy")

######### ######### ######### ######### ######### ######### ######### #########  
######### ######### ######### Load necessary libraries ######### #########  #########  
library(dplyr)
library(stringr) # str_detect function, super important
library('org.Mm.eg.db') # to get gene annotations for the scRNA-seq datasets 
library(data.table) # to read bigger csv faster
######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### 



######### ######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### #########  Define the functions  ######### ######### ######### ######### ######### 
Ct.Transform.To.mRNA.Counts_fun <- function(Dataset, Threshold, Preamplification, Wells, Start_Column){
  Transient_Dataset <- NULL
  Transient_Dataset <- 2**(Threshold  - Dataset[,Start_Column:ncol(Dataset)]-Preamplification)*Wells # 18 is for the pre-amplification, 96 is for the number of wells.
  Final_Dataset <- cbind(Dataset[,1:(Start_Column-1)], Transient_Dataset)
  
  return(Final_Dataset)
} # Transforms Ct values to mRNA counts or proportional to it 
TPM_fun <- function(Dataset,Genes,bool){
  
  Final <- NULL
  if(bool){
    for (names in rownames(Genes)){
      #print(names)
      vec <- Dataset[,names]/as.numeric(Genes[names,'GeneLength'])
      Final <- cbind(Final,vec)
    }
  }
  
  if(!bool){
    for (names in rownames(Genes)){
      #print(names)
      vec <- Dataset[,names]
      Final <- cbind(Final,vec)
    }
  }
  
  for (i in 1:nrow(Final)){
    #print(i)
    Final[i,] <- as.numeric(Final[i,]/sum(Final[i,],na.rm = TRUE))
  }
  colnames(Final) <- rownames(GeneLength)
  Final <- cbind(Dataset[,c(1,2)],Final)
  colnames(Final)[1] <- 'Code'
  colnames(Final)[2] <- 'Population'
  return(Final)
} # TPM normalisation for scRNA-seq
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### #########


######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### Pre-Process the data ######### ######### ######### 
{
  ######### scRT-qPCR datasets ########
  {
    # load the data
    Mutants.and.Ref16 <- as.data.frame(read.csv("DataYoan/MouseDataUsedData/MutantsAndRef16.csv",stringsAsFactors = FALSE)) # Mutants and Ref16, 48 genes
    References <- read.csv("DataYoan/MouseDataUsedData/References.csv",stringsAsFactors = FALSE) # Refs 32,64 and 90, 96 genes
    Gata6 <- read.csv("DataYoan/MouseDataUsedData/Gata6.csv",stringsAsFactors = FALSE) # Ref 16, 32, 64, 90, 3 genes (Gata6, Rps17, Rpl30 because of probes)
    
    # rename the first column for simplicity
    colnames(Mutants.and.Ref16)[1] <- 'Code' # We rename the first column to 'Code'
    colnames(References)[1] <- 'Code'
    colnames(Gata6)[1] <- 'Code'
    
    # We rename the Codes given to the cells for clarity
    Mutants.and.Ref16[str_detect(Mutants.and.Ref16$Code,'16C'),'Code'] <- as.character('C16') # we rename the cell type to C16 for 16C stage cells.
    Mutants.and.Ref16[str_detect(Mutants.and.Ref16$Code,'DKO'),'Code'] <- as.character('DKO') # we rename the cell type to DKO for Double KO Stage cells.
    Mutants.and.Ref16[str_detect(Mutants.and.Ref16$Code,'Gata6'),'Code'] <- as.character('G-/-') # we rename the cell type to G-/- for Gata6 knocked down stage cells.
    Mutants.and.Ref16[str_detect(Mutants.and.Ref16$Code,'Nanog'),'Code'] <- as.character('N-/-') # we rename the cell type to N-/- for Nanog Knocked down stage cells.
    Mutants.and.Ref16[str_detect(Mutants.and.Ref16$Code,'WT'),'Code'] <- as.character('WT') # we rename the cell type to C16 for 16 Cells.
    
    References[str_detect(References$Code,'32C'),'Code'] <- as.character('C32') # we rename the cell type to C32 for 32C stage cells.
    References[str_detect(References$Code,'64C'),'Code'] <- as.character('C64') # we to rename the cell type to C64 for 64C stage cells.
    References[str_detect(References$Code,'90C'),'Code'] <- as.character('C90') # we to rename the cell type to C90 for 90C stage cells.
    Gata6[str_detect(Gata6$Code,'16C'),'Code'] <- as.character('C16') # we rename the cell type to C16 for 16C stage cells.
    Gata6[str_detect(Gata6$Code,'32C'),'Code'] <- as.character('C32') # we rename the cell type to C32 for 32C stage cells.
    Gata6[str_detect(Gata6$Code,'64C'),'Code'] <- as.character('C64') # we rename the cell type to C64 for 64C stage cells.
    Gata6[str_detect(Gata6$Code,'90C'),'Code'] <- as.character('C90') # we rename the cell type to C90 for 90C stage cells.
    
    # This is specific to this dataset, the probe interacted with Gata6, so we rename Gata6 to Gata6_2 and inversely so as to be careful about the results obtained here.
    colnames(Mutants.and.Ref16)[which(colnames(Mutants.and.Ref16) == 'Gata6')] <- 'Gata6_2' 
    colnames(References)[which(colnames(References) == 'Gata6')] <- 'Gata6_2' 
    colnames(Gata6)[which(colnames(Gata6) == 'Gata6_2')] <- 'Gata6' 
    
    # Create a new column called Population with potential Epi or PrE progenitors
    Mutants.and.Ref16[Mutants.and.Ref16[,'Fgf4'] < 34.9,'EpiOrPrE'] <- as.character('Epi') # If Fgf4 < threshold then it becomes Epi
    Mutants.and.Ref16[Mutants.and.Ref16[,'Fgf4'] > 34.9,'EpiOrPrE'] <- as.character('PrE') # otherwise it becomes PrE
    Mutants.and.Ref16[,'Population'] <- paste(Mutants.and.Ref16$Code,Mutants.and.Ref16$EpiOrPrE,sep="")
    
    References[References[,'Fgf4'] <34.9,'EpiOrPrE'] <- as.character('Epi')
    References[References[,'Fgf4'] >34.9,'EpiOrPrE'] <- as.character('PrE')
    References[,'Population'] <- paste(References$Code,References$EpiOrPrE,sep="")
    
    # Create new datasets, one for the mutants, one for the references
    not.taken.indices <- c(which(colnames(Mutants.and.Ref16) == 'Code'), which(colnames(Mutants.and.Ref16) == 'Population'), which(colnames(Mutants.and.Ref16) == 'Name'),
                           which(colnames(Mutants.and.Ref16) == 'EpiOrPrE'))
    Mutants.and.Ref16.With.Population <- cbind(Mutants.and.Ref16[,'Code'],Mutants.and.Ref16[,'Population'],Mutants.and.Ref16[,-not.taken.indices])
    not.taken.indices <- c(which(colnames(References) == 'Code'), which(colnames(References) == 'Population'), which(colnames(References) == 'Name'),
                           which(colnames(References) == 'EpiOrPrE'))
    References.With.Population <- cbind(References[,'Code'],References[,'Population'],References[,-not.taken.indices])
    
    # give the appropriate names to the first two columns of these datasets
    colnames(Mutants.and.Ref16.With.Population)[1] <- 'Code'
    colnames(Mutants.and.Ref16.With.Population)[2] <- 'Population'
    colnames(References.With.Population)[1] <- 'Code'
    colnames(References.With.Population)[2] <- 'Population'
    
    ######### Create different Ct datasets for mutants and refs #########
    All.NoGata6 <- bind_rows(Mutants.and.Ref16.With.Population[str_detect(Mutants.and.Ref16.With.Population$Code,'C16'),],References.With.Population)
    Allegre.Ct <- bind_cols(All.NoGata6[str_detect(All.NoGata6$Code,'C'),],Gata6[,3:5])
    Mutants.Allegre.Ct <- Mutants.and.Ref16.With.Population[!str_detect(Mutants.and.Ref16.With.Population$Code,'C16'),]
    
    write.csv(Allegre.Ct,'FinalData/PreProcessedData/Allegre.Ct.csv', row.names = FALSE)
    write.csv(Mutants.Allegre.Ct,'FinalData/PreProcessedData/Mutants.Allegre.Ct.csv', row.names = FALSE)
    
    TwoCtValues <- 2**Allegre.Ct[,3:ncol(Allegre.Ct)] 
    Code <- Allegre.Ct[,'Code']
    Population <- Allegre.Ct[,'Population']
    Allegre.TwoCt <- cbind(Code,Population,TwoCtValues)
    
    TwoCtValues <- 2**Mutants.Allegre.Ct[,3:ncol(Mutants.Allegre.Ct)] 
    Code <- Mutants.Allegre.Ct[,'Code']
    Population <- Mutants.Allegre.Ct[,'Population']
    Mutants.Allegre.TwoCt <- cbind(Code,Population,TwoCtValues)
    
    write.csv(Allegre.TwoCt,'FinalData/PreProcessedData/Allegre.TwoCt.csv', row.names = FALSE)
    write.csv(Mutants.Allegre.TwoCt,'FinalData/PreProcessedData/Mutants.Allegre.TwoCt.csv', row.names = FALSE)
    
    rm(TwoCtValues,Code,Population)
    
    ######### Create the mRNA datasets #########
    Allegre.mRNA <- Ct.Transform.To.mRNA.Counts_fun(Dataset = Allegre.Ct, Threshold = 35, Preamplification = 18, Wells = 96, Start_Column =  3) 
    Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C16'),3:ncol(Allegre.mRNA)] <- Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C16'),3:ncol(Allegre.mRNA)]/2 # stage 16 has only 48 wells
    Mutants.Allegre.mRNA <- Ct.Transform.To.mRNA.Counts_fun(Dataset = Mutants.Allegre.Ct, Threshold = 35, Preamplification = 18, Wells = 48, Start_Column = 3)
    
    # delete variables that are not going to be used again
    rm(All.NoGata6,Gata6)
    rm(Mutants.and.Ref16,Mutants.and.Ref16.With.Population)
    rm(References,References.With.Population)
    rm(not.taken.indices)
    
    # mRNA counts datasets 
    write.csv(Allegre.mRNA,'FinalData/PreProcessedData/Allegre.mRNA.csv', row.names = FALSE)
    
    # Supplementary Datasets if needed
    write.csv(Mutants.Allegre.mRNA,'FinalData/PreProcessedData/Mutants.Allegre.mRNA.csv', row.names = FALSE)
    
    rm(Allegre.TwoCt,Mutants.Allegre.Ct,Mutants.Allegre.mRNA,Mutants.Allegre.TwoCt)
    print('Finished Allègre')
  } # Allègre Dataset - create mRNA counts.
  
  {
    setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Entropy") # first put ourselves in the directory where the data is, to load it.
    
    # load the data
    Guo.Before.Process <- as.data.frame(read.csv("DataGuo/GuoData/Guo.csv")) # actual Ct data from 1C to 64C
    
    # rename the first column for simplicity
    colnames(Guo.Before.Process)[1] <- 'Code'
    
    # We rename the Codes given to the cells for clarity.
    Guo.Before.Process[str_detect(Guo.Before.Process$Code,'1C'),'Code'] <- as.character('1C') # we rename the cell type to 1C for 1C stage cells.
    Guo.Before.Process[str_detect(Guo.Before.Process$Code,'^2C'),'Code'] <- as.character('2C') # we rename the cell type to 2C for 2C stage cells.
    Guo.Before.Process[str_detect(Guo.Before.Process$Code,'^4C'),'Code'] <- as.character('4C') # we rename the cell type to 4C for 4C stage cells.
    Guo.Before.Process[str_detect(Guo.Before.Process$Code,'8C'),'Code'] <- as.character('8C') # we rename the cell type to 8C for 8C stage cells.
    Guo.Before.Process[str_detect(Guo.Before.Process$Code,'16C'),'Code'] <- as.character('C16') # we rename the cell type to 16C for C16 stage cells.
    Guo.Before.Process[str_detect(Guo.Before.Process$Code,'32C'),'Code'] <- as.character('C32') # we rename the cell type to 32C for C32 stage cells.
    Guo.Before.Process[str_detect(Guo.Before.Process$Code,'64C'),'Code'] <- as.character('C64') # we rename the cell type to 64C for C64 stage cells.
    
    
    # Create a new column called Population with potential Epi or PrE progenitors
    Guo.Before.Process[(Guo.Before.Process[,'Fgf4'] < 28),'EpiOrPrE'] <- as.character('Epi')
    Guo.Before.Process[!(Guo.Before.Process[,'Fgf4'] < 28),'EpiOrPrE'] <- as.character('PrE')
    Guo.Before.Process[,'Population'] <- paste(Guo.Before.Process$Code,Guo.Before.Process$EpiOrPrE,sep="")
    
    # threshold chosen so we don't have TE cells when analysing, just to be sure. For more info, check supplementary material 'Guo.Id2.R'.
    threshold <- 24 
    Guo.In <- Guo.Before.Process[Guo.Before.Process[,'Id2'] > threshold,] # take only cells that are supposed to be the ICM cells and not TE
    
    # Create new datasets, one with TE cells, the other without
    not.taken.indices <- c(which(colnames(Guo.In) == 'Code'), which(colnames(Guo.In) == 'Population'), which(colnames(Guo.In) == 'EpiOrPrE'))
    Guo.Ct <- cbind(Guo.In[,'Code'],Guo.In[,'Population'],Guo.In[,-not.taken.indices])
    not.taken.indices <- c(which(colnames(Guo.Before.Process) == 'Code'), which(colnames(Guo.Before.Process) == 'Population'), which(colnames(Guo.Before.Process) == 'EpiOrPrE'))
    Guo.With.TE.Ct <- cbind(Guo.Before.Process[,'Code'],Guo.Before.Process[,'Population'],Guo.Before.Process[,-not.taken.indices])
    
    
    # give the appropriate names to the first two columns of these datasets
    colnames(Guo.Ct)[1] <- 'Code'
    colnames(Guo.Ct)[2] <- 'Population'
    
    colnames(Guo.With.TE.Ct)[1] <- 'Code'
    colnames(Guo.With.TE.Ct)[2] <- 'Population'
    
    write.csv(Guo.Ct,'FinalData/PreProcessedData/Guo.Ct.csv', row.names = FALSE)
    write.csv(Guo.With.TE.Ct,'FinalData/PreProcessedData/Guo.With.TE.Ct.csv', row.names = FALSE)
    
    
    TwoCtValues <- 2**Guo.Ct[,3:ncol(Guo.Ct)] 
    Code <- Guo.Ct[,'Code']
    Population <- Guo.Ct[,'Population']
    Guo.TwoCt <- cbind(Code,Population,TwoCtValues)
    
    
    TwoCtValues <- 2**Guo.With.TE.Ct[,3:ncol(Guo.With.TE.Ct)]
    Code <- Guo.With.TE.Ct[,'Code']
    Population <- Guo.With.TE.Ct[,'Population']
    Guo.With.TE.TwoCt <- cbind(Code,Population,TwoCtValues)
    
    write.csv(Guo.TwoCt,'FinalData/PreProcessedData/Guo.TwoCt.csv', row.names = FALSE)
    write.csv(Guo.With.TE.TwoCt,'FinalData/PreProcessedData/Guo.With.TE.TwoCt.csv', row.names = FALSE)
    
    
    #View(TwoCt.Guo.With.Population)
    rm(TwoCtValues,Code,Population)
    
    ######### Create the mRNA datasets #########
    Guo.mRNA.for.entropy <- Ct.Transform.To.mRNA.Counts_fun(Dataset = Guo.Ct, Threshold = 28, Preamplification = 18, Wells = 48, Start_Column = 3)
    Guo.With.TE.mRNA.for.entropy <- Ct.Transform.To.mRNA.Counts_fun(Dataset = Guo.With.TE.Ct, Threshold = 28, Preamplification = 18, Wells = 48, Start_Column = 3) 
    Guo.mRNA.for.CARDAMOM <- Ct.Transform.To.mRNA.Counts_fun(Dataset = Guo.Ct, Threshold = 28, Preamplification = 0, Wells = 1, Start_Column = 3)
    Guo.With.TE.mRNA.for.CARDAMOM <- Ct.Transform.To.mRNA.Counts_fun(Dataset = Guo.With.TE.Ct, Threshold = 28, Preamplification = 0, Wells = 1, Start_Column = 3) 
    
    # delete variables that are not going to be used again
    rm(Guo.In,Guo.Before.Process)
    rm(threshold,not.taken.indices)
    
    # mRNA counts datasets 
    write.csv(Guo.mRNA.for.entropy,'FinalData/PreProcessedData/Guo.mRNA.for.entropy.csv', row.names = FALSE)
    write.csv(Guo.mRNA.for.CARDAMOM,'FinalData/PreProcessedData/Guo.mRNA.for.CARDAMOM.csv', row.names = FALSE)
    
    # Supplementary Datasets if needed
    write.csv(Guo.With.TE.mRNA.for.entropy,'FinalData/PreProcessedData/Guo.mRNA.for.entropy.csv', row.names = FALSE)
    write.csv(Guo.With.TE.mRNA.for.CARDAMOM,'FinalData/PreProcessedData/Guo.With.TE.mRNA.for.CARDAMOM.csv', row.names = FALSE)
    
    rm(Guo.TwoCt,Guo.With.TE.Ct,Guo.With.TE.TwoCt)
    print('Finished Guo')
  } # Guo Dataset - create mRNA counts with TE and without TE.
  
  
  ######### scRNA-seq datasets ########
  {
    # load the data
    Goolam.Counts <- as.data.frame(fread('DataGoolam/Counts.csv', stringsAsFactors = FALSE, sep=';'))
    
    # Some pre-processing so we have the same format where rows are cell types and columns are genes
    Transposed.Goolam.Counts <- t(Goolam.Counts)
    Global_Stage <- rownames(Transposed.Goolam.Counts)
    Goolam.Counts.With.Code <- cbind(Global_Stage,Transposed.Goolam.Counts)
    rownames(Goolam.Counts.With.Code)[1] <- 'gene'
    colnames(Goolam.Counts.With.Code)  <- Goolam.Counts.With.Code['gene',]
    colnames(Goolam.Counts.With.Code)[1] <- 'Code'
    Goolam.Counts.With.Code <- as.data.frame(Goolam.Counts.With.Code[-c(1),])
    
    # We rename the Codes given to the cells for clarity.
    Goolam.Counts.With.Code[str_detect(Goolam.Counts.With.Code$Code,'Zygote'),'Code'] <- as.character('1C')
    Goolam.Counts.With.Code[str_detect(Goolam.Counts.With.Code$Code,'X2cell'),'Code'] <- as.character('2C')
    Goolam.Counts.With.Code[str_detect(Goolam.Counts.With.Code$Code,'4cell'),'Code'] <- as.character('4C')
    Goolam.Counts.With.Code[str_detect(Goolam.Counts.With.Code$Code,'8cell'),'Code'] <- as.character('8C')
    Goolam.Counts.With.Code[str_detect(Goolam.Counts.With.Code$Code,'16cell'),'Code'] <- as.character('C16')
    Goolam.Counts.With.Code[str_detect(Goolam.Counts.With.Code$Code,'32cell'),'Code'] <- as.character('C32')
    
    # Map the gene codes to the gene names in terms of ENSEMBL since they are all like this in other datasets
    Colnames <- as.data.frame(colnames(Goolam.Counts.With.Code[2:ncol(Goolam.Counts.With.Code)]))
    colnames(Colnames)[1] <- 'GeneCodes'
    Colnames$Symbol <- mapIds(org.Mm.eg.db, keys = Colnames$GeneCodes,keytype = 'ENSEMBL', column = 'SYMBOL')
    colnames(Goolam.Counts.With.Code)[2:ncol(Goolam.Counts.With.Code)] <- Colnames$Symbol
    
    # Create a new column called Population with potential Epi or PrE progenitors
    Goolam.Counts.With.Code[as.numeric(Goolam.Counts.With.Code[,'Fgf4']) > 0,'EpiOrPrE'] <- as.character('Epi')
    Goolam.Counts.With.Code[as.numeric(Goolam.Counts.With.Code[,'Fgf4']) == 0,'EpiOrPrE'] <- as.character('PrE')
    Goolam.Counts.With.Code[,'Population'] <- paste(Goolam.Counts.With.Code$Code,Goolam.Counts.With.Code$EpiOrPrE,sep="")
    
    # Create new dataset
    not.taken.indices <- c(which(colnames(Goolam.Counts.With.Code) == 'Code'), which(colnames(Goolam.Counts.With.Code) == 'Population'),
                           which(colnames(Goolam.Counts.With.Code) == 'EpiOrPrE'))
    Goolam.mRNA <- cbind(Goolam.Counts.With.Code[,'Code'], Goolam.Counts.With.Code[,'Population'], Goolam.Counts.With.Code[,-not.taken.indices])
    
    # give the appropriate names to the first two columns of the dataset
    colnames(Goolam.mRNA)[1] <- 'Code'
    colnames(Goolam.mRNA)[2] <- 'Population'
    
    # The counts table is written in 'str' format, and to make it numeric format we chose to write it and load it again. 
    # There might be other better options for this that are unknown to me
    fwrite(Goolam.mRNA,'FinalData/PreProcessedData/Goolam.mRNA.csv', row.names = FALSE)
    
    Goolam.mRNA <- as.data.frame(fread('FinalData/PreProcessedData/Goolam.mRNA.csv',stringsAsFactors = TRUE))
    # Gamma distributions cannot have 0 therefore, we do a pseudo-count by putting all zeros to one.
    # Goolam.mRNA[,3:ncol(Goolam.mRNA)] <- Goolam.mRNA[,3:ncol(Goolam.mRNA)] + 1
    
    # delete variables that are not going to be used again
    rm(Goolam.Counts,Goolam.Counts.With.Code)
    rm(Transposed.Goolam.Counts)
    rm(not.taken.indices,Global_Stage)
    rm(Colnames)
    
    colnames(Goolam.mRNA)[1] <- 'Code'
    colnames(Goolam.mRNA)[2] <- 'Population'
    
    # mRNA counts datasets 
    fwrite(Goolam.mRNA,'FinalData/PreProcessedData/Goolam.mRNA.csv', row.names = FALSE)
    
  } # Goolam Dataset: Here it's smartseq2, so don't need to divide by the gene length for coverage when doing TPM.
  { 
    # first write the gene lengths
    Bmp4 <- -(46620977-46628126)
    Cdx2 <- -(147237710-147244059)
    Dab2 <- -(6329269-6470196)
    Esrrb <- -(86407891-86568402)
    Fgf4 <- -(144415123-144418982)
    Fgfr2 <- -(129764181-129868538)
    Gata3 <- -(9861889-9892762)
    Gata4 <- -(63436371-63509141)
    Gata6 <- -(11052510-11085636)
    Hnf4a <- -(163348731-163414827)
    Id2 <- -(25143798-25146091)
    Klf2 <- -(73072906-73075498)
    Klf4 <- -(55527143-55532466)
    Klf5 <- -(99536127-99550848)
    Nanog <- -(122684448-122691592)
    Pdgfa <- -(138961769-138983125)
    Pdgfra <- -(75312953-75358876)
    Pecam1 <- -(106545039-106606107)
    Pou5f1 <- -(35816929-35821674)
    Sox2 <- -(34704554-34706610)
    Sox17 <- -(4561154-4567577)
    
    # create a dataframe with the length of the genes for normalization
    GeneLength <- as.data.frame(c(Bmp4,Cdx2,Dab2,Esrrb,Fgf4,Fgfr2,Gata3,Gata4,Gata6,Hnf4a,Id2,Klf2,Klf4,Klf5,Nanog,Pdgfa,Pdgfra,Pecam1,Pou5f1,Sox2,Sox17))
    rownames(GeneLength) <- c('Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6','Hnf4a','Id2','Klf2','Klf4','Klf5','Nanog','Pdgfa',
                              'Pdgfra','Pecam1','Pou5f1','Sox2','Sox17')
    colnames(GeneLength) <- c('GeneLength')
    GeneLength$GeneLength <- GeneLength/1000 # transform to kb
    
    # take only the genes that intersect between the RT-qPCR datasets (since they posess the lesser number of genes)
    genes.32 <- intersect(names(Allegre.mRNA),names(Guo.mRNA.for.entropy))
    
    # Generate the TPMs using the TPM_fun
    Transient.Goolam <- Goolam.mRNA[,genes.32]
    Goolam.TPM <- TPM_fun(Dataset = Transient.Goolam, Genes = GeneLength, bool = FALSE)
    
    # delete variables that are not going to be used again
    rm(Bmp4,Cdx2,Dab2,Esrrb,Fgf4,Fgfr2,Gata3,Gata4,Gata6,Hnf4a,Id2,Klf2,Klf4,Klf5,Nanog,Pdgfa,Pdgfra,Pecam1,Pou5f1,Sox2,Sox17)
    rm(Transient.Goolam)
    rm(intersection.RT.qPCR)
    rm(genes.32)
    
    # TPM datasets
    fwrite(Goolam.TPM,'FinalData/PreProcessedData/Goolam.TPM.csv', row.names = FALSE)
    
    print('Finished Goolam')
    
  } # TPM Normalization for Goolam's scRNA-seq dataset, no division by gene length here
  
  {
    setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Entropy")
    
    
    # load the data
    Posfai.Annotations <- as.data.frame(fread("DataPosfai/Annotation.csv", stringsAsFactors = FALSE)) # This is the annotation data for the stages
    Posfai.Counts <- as.data.frame(fread("DataPosfai/Count.csv", stringsAsFactors = TRUE)) # These are the actual count numbers
    
    # Some pre-processing so we have the same format where rows are cell types and columns are genes
    Transposed.Posfai.Counts <- t(Posfai.Counts)
    colnames(Transposed.Posfai.Counts) <- Transposed.Posfai.Counts[1,]
    Transposed.Posfai.Counts <- Transposed.Posfai.Counts[-c(1),]
    
    # little code to change the cell name to its cellular stage (16C, 32C or 64C using the annotations)
    vec <- c()
    for(i in 1:nrow(Transposed.Posfai.Counts)){
      for(j in 1:nrow(Posfai.Annotations)){
        if(rownames(Transposed.Posfai.Counts)[i] == Posfai.Annotations[j,'cells']){
          vec <- c(vec,Posfai.Annotations[j,'global_Stage'])
        }
      }
    }
    Posfai.Counts.With.Code <- as.data.frame(cbind(vec,Transposed.Posfai.Counts))
    colnames(Posfai.Counts.With.Code)[1] <- 'Code'
    
    # We rename the Codes given to the cells for clarity.
    Posfai.Counts.With.Code[str_detect(Posfai.Counts.With.Code$Code,'E16'),'Code'] <- as.character('C16')
    Posfai.Counts.With.Code[str_detect(Posfai.Counts.With.Code$Code,'E32'),'Code'] <- as.character('C32')
    Posfai.Counts.With.Code[str_detect(Posfai.Counts.With.Code$Code,'E64'),'Code'] <- as.character('C64')
    
    custom_order <- c("C16", "C32", "C64")
    Posfai.Counts.With.Code$Code <- factor(Posfai.Counts.With.Code$Code, levels = custom_order)
    Posfai.Counts.With.Code <- Posfai.Counts.With.Code[order(Posfai.Counts.With.Code$Code), ]
    
    
    # Create a new column called Population with potential Epi or PrE progenitors
    Posfai.Counts.With.Code[as.numeric(Posfai.Counts.With.Code[,'Fgf4']) > 0,'EpiOrPrE'] <- as.character('Epi')
    Posfai.Counts.With.Code[as.numeric(Posfai.Counts.With.Code[,'Fgf4']) == 0,'EpiOrPrE'] <- as.character('PrE')
    Posfai.Counts.With.Code[,'Population'] <- paste(Posfai.Counts.With.Code$Code,Posfai.Counts.With.Code$EpiOrPrE,sep="")
    
    # Create new dataset
    not.taken.indices <- c(which(colnames(Posfai.Counts.With.Code) == 'Code'), which(colnames(Posfai.Counts.With.Code) == 'Population'),
                           which(colnames(Posfai.Counts.With.Code) == 'EpiOrPrE'))
    Posfai.mRNA <- cbind(Posfai.Counts.With.Code[,'Code'], Posfai.Counts.With.Code[,'Population'], Posfai.Counts.With.Code[,-not.taken.indices])
    
    # give the appropriate names to the first two columns of the dataset
    colnames(Posfai.mRNA)[1] <- 'Code'
    colnames(Posfai.mRNA)[2] <- 'Population'
    
    # The counts table is written in 'str' format, and to make it numeric format we chose to write it and load it again. 
    # There might be other better options for this that are unknown to me
    fwrite(Posfai.mRNA,'FinalData/PreProcessedData/Posfai.mRNA.csv', row.names = FALSE)
    
    Posfai.mRNA <- as.data.frame(fread('FinalData/PreProcessedData/Posfai.mRNA.csv',stringsAsFactors = TRUE))
    # Gamma distributions cannot have 0 therefore, we do a pseudo-count by putting all zeros to one.
    # Posfai.mRNA[,3:ncol(Posfai.mRNA)] <- Posfai.mRNA[,3:ncol(Posfai.mRNA)] + 1
    
    # delete variables that are not going to be used again
    rm(Posfai.Annotations,Posfai.Counts)
    rm(Transposed.Posfai.Counts,Posfai.Counts.With.Code)
    rm(i,j,not.taken.indices,vec)
    
    # mRNA counts datasets 
    fwrite(Posfai.mRNA,'FinalData/PreProcessedData/Posfai.mRNA.csv', row.names = FALSE)
    
  } # Posfai Dataset: Here it's smartseq2, so don't need to divide by the gene length for coverage when doing TPM.
  { 
    # first write the gene lengths
    Bmp4 <- -(46620977-46628126)
    Cdx2 <- -(147237710-147244059)
    Dab2 <- -(6329269-6470196)
    Esrrb <- -(86407891-86568402)
    Fgf4 <- -(144415123-144418982)
    Fgfr2 <- -(129764181-129868538)
    Gata3 <- -(9861889-9892762)
    Gata4 <- -(63436371-63509141)
    Gata6 <- -(11052510-11085636)
    Hnf4a <- -(163348731-163414827)
    Id2 <- -(25143798-25146091)
    Klf2 <- -(73072906-73075498)
    Klf4 <- -(55527143-55532466)
    Klf5 <- -(99536127-99550848)
    Nanog <- -(122684448-122691592)
    Pdgfa <- -(138961769-138983125)
    Pdgfra <- -(75312953-75358876)
    Pecam1 <- -(106545039-106606107)
    Pou5f1 <- -(35816929-35821674)
    Sox2 <- -(34704554-34706610)
    Sox17 <- -(4561154-4567577)
    
    # create a dataframe with the length of the genes for normalization
    GeneLength <- as.data.frame(c(Bmp4,Cdx2,Dab2,Esrrb,Fgf4,Fgfr2,Gata3,Gata4,Gata6,Hnf4a,Id2,Klf2,Klf4,Klf5,Nanog,Pdgfa,Pdgfra,Pecam1,Pou5f1,Sox2,Sox17))
    rownames(GeneLength) <- c('Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6','Hnf4a','Id2','Klf2','Klf4','Klf5','Nanog','Pdgfa',
                              'Pdgfra','Pecam1','Pou5f1','Sox2','Sox17')
    colnames(GeneLength) <- c('GeneLength')
    GeneLength$GeneLength <- GeneLength/1000 # transform to kb
    
    # take only the genes that intersect between the RT-qPCR datasets (since they posess the lesser number of genes)
    genes.32 <- intersect(names(Allegre.mRNA),names(Guo.mRNA.for.entropy))
    
    # Generate the TPMs using the TPM_fun
    Transient.Posfai <- Posfai.mRNA[,genes.32]
    Posfai.TPM <- TPM_fun(Dataset = Transient.Posfai, Genes = GeneLength, bool = FALSE) # bool = False because no need to divide by gene length for coverage
    
    # delete variables that are not going to be used again
    rm(Bmp4,Cdx2,Dab2,Esrrb,Fgf4,Fgfr2,Gata3,Gata4,Gata6,Hnf4a,Id2,Klf2,Klf4,Klf5,Nanog,Pdgfa,Pdgfra,Pecam1,Pou5f1,Sox2,Sox17)
    rm(Transient.Posfai)
    rm(genes.32)
    rm(GeneLength)
    
    # TPM datasets
    fwrite(Posfai.TPM,'FinalData/PreProcessedData/Posfai.TPM.csv', row.names = FALSE)
    
    print('Finished Posfai')
    
  } # TPM Normalization for Posfai's scRNA-seq dataset, no division by gene length here
  
  { 
    
    # load the data
    Wang.Counts <- as.data.frame(fread('DataWang/CountTableWang.csv'))

    # Some pre-processing so we have the same format where rows are cell types and columns are genes
    Transposed.Wang.Counts <- t(Wang.Counts)
    Global_Stage <- rownames(Transposed.Wang.Counts)
    Wang.Counts.With.Code <- cbind(Global_Stage,Transposed.Wang.Counts)
    colnames(Wang.Counts.With.Code)  <- Wang.Counts.With.Code['gene',]
    colnames(Wang.Counts.With.Code)[1] <- 'Code'
    Wang.Counts.With.Code <- as.data.frame(Wang.Counts.With.Code[-c(1),])
    
    # We rename the Codes given to the cells for clarity.
    Wang.Counts.With.Code[str_detect(Wang.Counts.With.Code$Code,'Zygote'),'Code'] <- as.character('1C')
    Wang.Counts.With.Code[str_detect(Wang.Counts.With.Code$Code,'4cell'),'Code'] <- as.character('4C')
    Wang.Counts.With.Code[str_detect(Wang.Counts.With.Code$Code,'Late4cell'),'Code'] <- as.character('4LC')
    Wang.Counts.With.Code[str_detect(Wang.Counts.With.Code$Code,'8cell'),'Code'] <- as.character('8C')
    Wang.Counts.With.Code[str_detect(Wang.Counts.With.Code$Code,'16cell'),'Code'] <- as.character('C16')
    Wang.Counts.With.Code[str_detect(Wang.Counts.With.Code$Code,'32cell'),'Code'] <- as.character('C32')
    Wang.Counts.With.Code[str_detect(Wang.Counts.With.Code$Code,'2cell_embryo'),'Code'] <- as.character('2C')
    
    # Map the gene codes to the gene names in terms of ENSEMBL since they are all like this in other datasets
    Colnames <- as.data.frame(colnames(Wang.Counts.With.Code[2:ncol(Wang.Counts.With.Code)]))
    colnames(Colnames)[1] <- 'GeneCodes'
    Colnames$Symbol <- mapIds(org.Mm.eg.db, keys = Colnames$GeneCodes,keytype = 'ENSEMBL', column = 'SYMBOL')
    colnames(Wang.Counts.With.Code)[2:ncol(Wang.Counts.With.Code)] <- Colnames$Symbol
    
    # Create a new column called Population with potential Epi or PrE progenitors
    Wang.Counts.With.Code[as.numeric(Wang.Counts.With.Code[,'Fgf4']) > 0,'EpiOrPrE'] <- as.character('Epi')
    Wang.Counts.With.Code[as.numeric(Wang.Counts.With.Code[,'Fgf4']) == 0,'EpiOrPrE'] <- as.character('PrE')
    Wang.Counts.With.Code[,'Population'] <- paste(Wang.Counts.With.Code$Code,Wang.Counts.With.Code$EpiOrPrE,sep="")
    
    # Create new dataset
    not.taken.indices <- c(which(colnames(Wang.Counts.With.Code) == 'Code'), which(colnames(Wang.Counts.With.Code) == 'Population'),
                           which(colnames(Wang.Counts.With.Code) == 'EpiOrPrE'))
    Wang.mRNA <- cbind(Wang.Counts.With.Code[,'Code'], Wang.Counts.With.Code[,'Population'], Wang.Counts.With.Code[,-not.taken.indices])
    
    # give the appropriate names to the first two columns of the dataset
    colnames(Wang.mRNA)[1] <- 'Code'
    colnames(Wang.mRNA)[2] <- 'Population'
    
    # The counts table is written in 'str' format, and to make it numeric format we chose to write it and load it again. 
    # There might be other better options for this that are unknown to me
    fwrite(Wang.mRNA,'FinalData/PreProcessedData/Wang.mRNA.csv', row.names = FALSE)
    
    Wang.mRNA <- as.data.frame(fread('FinalData/PreProcessedData/Wang.mRNA.csv',stringsAsFactors = TRUE))
    
    # Gamma distributions cannot have 0 therefore, we do a pseudo-count by putting all zeros to 1
    # Wang.mRNA[,3:ncol(Wang.mRNA)] <- Wang.mRNA[,3:ncol(Wang.mRNA)] + 1
    
    # delete variables that are not going to be used again
    rm(Wang.Counts,Wang.Counts.With.Code)
    rm(Transposed.Wang.Counts)
    rm(not.taken.indices,Global_Stage)
    rm(Colnames)
    
    # mRNA counts datasets 
    fwrite(Wang.mRNA,'FinalData/PreProcessedData/Wang.mRNA.csv', row.names = FALSE)
    
  } # Wang Dataset: Here no smartseq2, so should need to divide by the gene length for coverage when doing TPM.
  {
    # first write the gene lengths
    Bmp4 <- -(46620977-46628126)
    Cdx2 <- -(147237710-147244059)
    Dab2 <- -(6329269-6470196)
    Esrrb <- -(86407891-86568402)
    Fgf4 <- -(144415123-144418982)
    Fgfr2 <- -(129764181-129868538)
    Gata3 <- -(9861889-9892762)
    Gata4 <- -(63436371-63509141)
    Gata6 <- -(11052510-11085636)
    Hnf4a <- -(163348731-163414827)
    Id2 <- -(25143798-25146091)
    Klf2 <- -(73072906-73075498)
    Klf4 <- -(55527143-55532466)
    Klf5 <- -(99536127-99550848)
    Nanog <- -(122684448-122691592)
    Pdgfa <- -(138961769-138983125)
    Pdgfra <- -(75312953-75358876)
    Pecam1 <- -(106545039-106606107)
    Pou5f1 <- -(35816929-35821674)
    Sox2 <- -(34704554-34706610)
    Sox17 <- -(4561154-4567577)
    
    # create a dataframe with the length of the genes for normalization
    GeneLength <- as.data.frame(c(Bmp4,Cdx2,Dab2,Esrrb,Fgf4,Fgfr2,Gata3,Gata4,Gata6,Hnf4a,Id2,Klf2,Klf4,Klf5,Nanog,Pdgfa,Pdgfra,Pecam1,Pou5f1,Sox2,Sox17))
    rownames(GeneLength) <- c('Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6','Hnf4a','Id2','Klf2','Klf4','Klf5','Nanog','Pdgfa',
                              'Pdgfra','Pecam1','Pou5f1','Sox2','Sox17')
    colnames(GeneLength) <- c('GeneLength')
    GeneLength$GeneLength <- GeneLength/1000 # transform to kb
    
    # take only the genes that intersect between the RT-qPCR datasets (since they posess the lesser number of genes)
    intersection.RT.qPCR <- intersect(names(Allegre.mRNA),names(Guo.mRNA.for.entropy))
    
    # Generate the TPMs using the TPM_fun
    Transient.Wang <- Wang.mRNA[,intersection.RT.qPCR]
    Wang.TPM <- TPM_fun(Dataset = Transient.Wang, Genes = GeneLength, bool = TRUE) # bool = TRUE because we need to divide by gene length
    
    # delete variables that are not going to be used again
    rm(Bmp4,Cdx2,Dab2,Esrrb,Fgf4,Fgfr2,Gata3,Gata4,Gata6,Hnf4a,Id2,Klf2,Klf4,Klf5,Nanog,Pdgfa,Pdgfra,Pecam1,Pou5f1,Sox2,Sox17)
    rm(Transient.Wang)
    rm(intersection.RT.qPCR)
    rm(GeneLength)
    
    # TPM datasets
    fwrite(Wang.TPM,'FinalData/PreProcessedData/Wang.TPM.csv', row.names = FALSE)
    
    print('Finished Wang')
    
  } # TPM Normalization for Wang's scRNA-seq dataset, division by gene length to take care of coverage
  
  rm(custom_order,Ct.Transform.To.mRNA.Counts_fun,TPM_fun)
  
} # All Pre-processing of the datasets
######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### #########