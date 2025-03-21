######################### Standard Analysis File #########################
## In this file, we perform the number of cell barplots, PCAs, correlations analysis of the different datasets

######### setting the directory for loading data ######### 
setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Entropy")
load('FinalData/InterCellularEntropy/EntropyGamma/PreProcessedData.RData') 


######### loading of libraries #########
library(dplyr)
library(stringr) # str_detect function, super important
library(ggbiplot) # PCA plots
library(gridExtra) # Extra refining for 2x2 plots
library(pheatmap) # Pheatmap package for nice heatmap production
library(RColorBrewer) # RColorBrewer package for nice colors on the heatmap
library(tidyverse) # for pivot longer
library(ggplot2) # for ggplots
library(viridis) # for inferno

######### Functions ######### 
PCA_fun <- function(Dataset,genes,title,values_used,BOOL_axes, position = 'none'){
  
  PCA <- prcomp(Dataset[,genes], center = TRUE, scale = TRUE)
  
  p <- ggbiplot(PCA, ellipse = FALSE, var.axes = BOOL_axes, varname.size = 10) +
    geom_point(aes(color = factor(x = Dataset$Code, levels = c('8C','C16','C32','C64','C90'), labels = c('8C','16C','32C','64C','90C'))), size = 4)+
    scale_colour_manual(name = "Cellular Stage", values = values_used) +
    ggtitle(title) +
    theme_minimal() +
    xlim(-3.5,3.5) + 
    ylim(-3.5,3.5) +
    theme(legend.position = position, legend.title=element_text(size=30),legend.text=element_text(size=40),
          plot.title = element_text(hjust = 0.5, size = 45),
          axis.title = element_text(size = 25), axis.text = element_text(size = 25))
  
  return(p)
} # Plots of PCA
PCA_fun_with_Fgf4_index <- function(Dataset, genes, title, values_used, BOOL_axes, position = 'none') {
  
  PCA <- prcomp(Dataset[, genes], center = TRUE, scale = TRUE)
  
  Dataset$Index <- ifelse(str_detect(Dataset$Population,"Epi"), "Fgf4+", "Fgf4-")
  
  
  p <- ggbiplot(PCA, ellipse = FALSE, var.axes = BOOL_axes, varname.size = 10) +
    geom_point(aes(
      color = factor(Dataset$Code, levels = c('8C','C16','C32','C64','C90'), 
                     labels = c('8C','16C','32C','64C','90C')),
      shape = factor(Dataset$Index)  # Shape represents cell type
    ), size = 4) +
    scale_colour_manual(name = "Cellular Stage", values = values_used) +
    scale_shape_manual(name = "Cell Type", values = c(16, 17)) +  # Adjust shapes as needed
    ggtitle(title) +
    theme_minimal() +
    xlim(-3.5, 3.5) + 
    ylim(-3.5, 3.5) +
    theme(legend.position = position, legend.title = element_text(size = 30),
          legend.text = element_text(size = 40),
          plot.title = element_text(hjust = 0.5, size = 45),
          axis.title = element_text(size = 25), axis.text = element_text(size = 25))
  
  return(p)
}
Subset_fun <- function(Dataset,number_of_subsets,Population_Or_Code){
  data <- Dataset
  
  # Create a list to store the subsets
  subsets_list <- list()
  
  # Loop over each unique value of Code and create subsets of that value
  for (pop_code_value in unique(data[,Population_Or_Code])) {
    # Subset the data for the current code
    subset_data <- data[data[,Population_Or_Code] == pop_code_value, ]
    
    # Divide the subset into three parts
    subset_list <- split(subset_data, rep(1:number_of_subsets, each = nrow(subset_data)/number_of_subsets))
    
    # Store the subsets in the subsets_list
    subsets_list[[pop_code_value]] <- subset_list
  }
  return(subsets_list)
} # Subsetting function to create list of based on Code or Population
Expression.Correlation_fun <- function(Dataset,genes,threshold,threshold_value){
  
  data.frame <- as.data.frame(cor(dplyr::select(Dataset, genes),method = 'spearman'))
  
  if(threshold == TRUE){
    data.frame[data.frame <= threshold_value & data.frame >= -threshold_value] <- 0
  }
  
  return(data.frame)
} # Creating Correlation Heatmaps 
Barplots_fun <- function(Dataset,genes,Cellular.Stage,threshold){
  threshold_bool <- FALSE
  threshold_used <- NA
  df <- NULL
  for (names in Cellular.Stage){
    a <- Expression.Correlation_fun(Dataset[str_detect(Dataset$Code,names),], genes,threshold = threshold_bool, threshold_value = threshold_used)
    Correlation.Counts <- as.data.frame(table(a > threshold))
    Anti.Correlation.Counts <- as.data.frame(table(a < -threshold))
    
    row <- rbind((Correlation.Counts[2,'Freq']-length(genes))/2,(Anti.Correlation.Counts[2,'Freq'])/2)
    df <- cbind(df,row)
    
  }
  df <- as.matrix(df)
  colnames(df) <- c('8C','16C','32C','64C','90C') 
  return(df)
} # Function to create the barplots of correlations/anti-correlations with a threshold


######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### Genes chosen and degradation rates ######### ######### 
genes.16 <- c('Bmp4','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6',
              'Id2','Nanog','Pdgfra','Pecam1','Pou5f1','Sox17','Sox2')
genes.21 <- c('Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6',
              'Hnf4a','Id2','Klf2', 'Klf4','Klf5','Nanog','Pdgfa','Pdgfra','Pecam1','Pou5f1','Sox17','Sox2')
genes.21.with.code.and.pop <- c('Code','Population','Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6',
                                'Hnf4a','Id2','Klf2', 'Klf4','Klf5','Nanog','Pdgfa','Pdgfra','Pecam1','Pou5f1','Sox17','Sox2')
######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### 


######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### PCAs and Correlations ######### ######### ######### 
{
  genes <- c('Nanog','Pou5f1','Actb','Gapdh')
  df <- Guo.Ct[!(str_detect(Guo.Ct$Code,'4C') |
               str_detect(Guo.Ct$Code,'2C')),] # Allegre.Ct
  # Exclude the 'Code' and 'Population' columns for plotting
  columns_to_plot <- colnames(df)[colnames(df) %in% genes]
  
  # Reshape the dataset to a long format for easier plotting with ggplot2
  df_long <- df %>%
    pivot_longer(
      cols = all_of(columns_to_plot),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    mutate(Gene = factor(Gene, levels = genes))
  
  # Create the boxplot
  ggplot(df_long, aes(x = Code, y = Expression, fill = Code)) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~ Gene, scales = "free_y", nrow = 1) +  # Separate plot for each gene
    theme_minimal() +
    scale_x_discrete(labels = c('8C','16C', '32C', '64C'), drop = TRUE) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5)
    ) +
    labs(
      title = "Guo Dataset",
      x = "Cellular Stage",
      y = "Ct Value"
    ) +
    ylim(c(12,28))
  
  ggsave(filename = 'GGplotsFigs/Barplots/Ct_Values_Guo.pdf', width = 5, height = 2.5)
} # Boxplots of Ct values

{
  ######### PCA genes in common at 16C ######### 
  {
    genes.used <- genes.16
    p1 <- PCA_fun(Dataset = Allegre.mRNA,genes = genes.used,
                  title = 'Allègre',values_used =  c("purple", "red","orange",'green'), BOOL_axes = FALSE)
    p2 <- PCA_fun(Dataset = Goolam.TPM[str_detect(Goolam.TPM$Code,'^C') | str_detect(Goolam.TPM$Code,'8C') ,],genes = genes.used,
                  title = 'Goolam',values_used = c("darkblue", "purple", "red"), BOOL_axes = FALSE)
    p3 <- PCA_fun(Dataset = Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'^C') | str_detect(Guo.mRNA.for.entropy$Code,'8C'),],genes = genes.used,
                  title = 'Guo',values_used =  c("darkblue", "purple", "red","orange"), BOOL_axes = FALSE)
    p4 <- PCA_fun(Dataset = Posfai.TPM,genes = genes.used,
                  title = 'Posfai', values_used = c("purple", "red","orange"), BOOL_axes = FALSE)
    p5 <- PCA_fun(Dataset = Wang.TPM[str_detect(Wang.TPM$Code,'^C') | str_detect(Wang.TPM$Code,'8C'),],genes = genes.used,
                  title = 'Wang',values_used = c("darkblue", "purple", "red"), BOOL_axes = FALSE)
    
    p6 <- PCA_fun(Dataset = Allegre.mRNA,genes = genes.used,
                  title = 'Allègre', values_used =  c("purple", "red","orange",'green'), BOOL_axes = TRUE)
    
    p <- grid.arrange(p1, p6, p2, p3, p4, p5, nrow = 3,
                      layout_matrix = rbind(c(1, 2, 3), c(4, 5, 6)))
    p
    ggsave(paste('GGplotsFigs/EntropyGamma/PCA/PCA8C.pdf', sep = ''), p, width = 19, height = 25)
    rm(genes.used)
    rm(p,p1,p2,p3,p4,p5,p6)
    
  } # PCA with genes in common at 16C + PCA with var_axis on Allègre
  
  {
    genes.used <- genes.21
    
    p1 <- PCA_fun(Dataset = Allegre.mRNA[!str_detect(Allegre.mRNA$Code,'C16'),], genes = genes.used,
                  title = 'Allègre', values_used = c("red","orange",'green'), BOOL_axes = FALSE)
    p2 <- PCA_fun(Dataset = Goolam.TPM[str_detect(Goolam.TPM$Code,'^C') | str_detect(Goolam.TPM$Code,'8C') ,],genes = genes.used,
                  title = 'Goolam', c("darkblue", "purple", "red"), BOOL_axes = FALSE)
    p3 <- PCA_fun(Dataset = Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'^C') | str_detect(Guo.mRNA.for.entropy$Code,'8C') ,], genes = genes.used,
                  title = 'Guo', values_used = c("darkblue", "purple", "red","orange"), BOOL_axes = FALSE)
    p4 <- PCA_fun(Dataset = Posfai.TPM[str_detect(Posfai.TPM$Code,'^C'),], genes = genes.used,
                  title = 'Posfai', values_used = c("purple", "red","orange"), BOOL_axes = FALSE)
    p5 <- PCA_fun(Dataset = Wang.TPM[str_detect(Wang.TPM$Code,'^C') | str_detect(Wang.TPM$Code,'8C') ,], genes = genes.used,
                  title = 'Wang', values_used = c("darkblue", "purple", "red"), BOOL_axes = FALSE)
    
    p6 <- PCA_fun(Dataset = Allegre.mRNA[!str_detect(Allegre.mRNA$Code,'C16'),], genes = genes.used,
                  title = 'Allègre', values_used =  c("red","orange",'green'), BOOL_axes = TRUE)
    
    p <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2)
    
    ggsave(paste('GGplotsFigs/EntropyGamma/PCA/PCAC32.pdf', sep = ''), p, width = 24, height = 18)
    rm(genes.used)
    rm(p,p1,p2,p3,p4,p5)
    
  } # PCA with genes in common at 32C 
  
  ######### PCA genes in common at 32C ######### 
  {
    genes.used <- genes.16
    
    p1 <- PCA_fun(Dataset = Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C16') |
                                         str_detect(Allegre.mRNA$Code,'C90'),], genes = genes.used,
                  title = 'Allègre',values_used = c("purple","green"), BOOL_axes = FALSE)
    p2 <- PCA_fun(Dataset = Goolam.TPM[str_detect(Goolam.TPM$Code,'8C') |
                                         str_detect(Goolam.TPM$Code,'C32'),] , genes = genes.used,
                  title = 'Goolam', values_used =  c("darkblue", "red"), BOOL_axes = FALSE)
    p3 <- PCA_fun(Dataset = Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'8C') |
                                       str_detect(Guo.mRNA.for.entropy$Code,'C64'),], genes = genes.used,
                  title = 'Guo',values_used = c("darkblue","orange"), BOOL_axes = FALSE)
    p4 <- PCA_fun(Dataset = Posfai.TPM[str_detect(Posfai.TPM$Code,'C16') |
                                         str_detect(Posfai.TPM$Code,'C64'),], genes = genes.used,
                  title = 'Posfai',values_used = c("purple","orange"), BOOL_axes = FALSE)
    p5 <- PCA_fun(Dataset = Wang.TPM[str_detect(Wang.TPM$Code,'8C') |
                                       str_detect(Wang.TPM$Code,'C32'),], genes = genes.used,
                  title = 'Wang', values_used = c("darkblue","red"), BOOL_axes = FALSE)
    
    
    p6 <- PCA_fun(Dataset = Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C16') |
                                           str_detect(Allegre.mRNA$Code,'C90'),], genes = genes.used,
                  title = 'Allègre', values_used =  c("purple",'green'), BOOL_axes = TRUE)
    
    
    p <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2)
    
    ggsave(paste('GGplotsFigs/EntropyGamma/PCA/PCAEP16.pdf', sep = ''), p, width = 30, height= 25)
    rm(genes.used)
    rm(p,p1,p2,p3,p4,p5) # PCA for the extreme time points
    
  } # PCA ETP with genes in common at 16C
  
  {
    genes.used <- genes.21
    
    p1 <- PCA_fun(Dataset = Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C32') |
                                           str_detect(Allegre.mRNA$Code,'C90'),], genes = genes.used,
                  title = 'Allègre', values_used = c("red","green"), BOOL_axes = FALSE)
    p2 <- PCA_fun(Dataset = Goolam.TPM[str_detect(Goolam.TPM$Code,'8C') |
                                         str_detect(Goolam.TPM$Code,'C32'),] , genes = genes.used,
                  title = 'Goolam', values_used =  c("darkblue", "red"), BOOL_axes = FALSE)
    p3 <- PCA_fun(Dataset = Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'8C') |
                                       str_detect(Guo.mRNA.for.entropy$Code,'C64'),], genes = genes.used,
                  title = 'Guo', values_used = c("darkblue","orange"), BOOL_axes = FALSE)
    p4 <- PCA_fun(Dataset = Posfai.TPM[str_detect(Posfai.TPM$Code,'C16') |
                                         str_detect(Posfai.TPM$Code,'C64'),], genes = genes.used,
                  title = 'Posfai', values_used = c("purple","orange"), BOOL_axes = FALSE)
    p5 <- PCA_fun(Dataset = Wang.TPM[str_detect(Wang.TPM$Code,'8C') |
                                       str_detect(Wang.TPM$Code,'C32'),], genes = genes.used,
                  title = 'Wang', values_used = c("darkblue","red"), BOOL_axes = FALSE)
    
    p6 <- PCA_fun(Dataset = Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C32') |
                                           str_detect(Allegre.mRNA$Code,'C90'),], genes = genes.used,
                  title = 'Allègre', values_used =  c("red",'green'), BOOL_axes = TRUE)
    
    
    p <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2)
    ggsave(paste('GGplotsFigs/EntropyGamma/PCA/PCAEP32.pdf', sep = ''), p, width = 30, height= 25)
    rm(genes.used)
    rm(p,p1,p2,p3,p4,p5) # PCA for the extreme time points
    
  } # PCA ETP with genes in common at 32C
  
  {
    genes.used <- genes.21
    
    p1 <- PCA_fun_with_Fgf4_index(Dataset = Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C32') |
                                           str_detect(Allegre.mRNA$Code,'C90'),], genes = genes.used,
                                  title = 'Allègre', values_used = c("red","green"), BOOL_axes = FALSE)
    p2 <- PCA_fun_with_Fgf4_index(Dataset = Goolam.TPM[str_detect(Goolam.TPM$Code,'8C') |
                                         str_detect(Goolam.TPM$Code,'C32'),] , genes = genes.used,
                                  title = 'Goolam', values_used =  c("darkblue", "red"), BOOL_axes = FALSE)
    p3 <- PCA_fun_with_Fgf4_index(Dataset = Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'8C') |
                                                   str_detect(Guo.mRNA.for.entropy$Code,'C64'),], genes = genes.used,
                                  title = 'Guo', values_used = c("darkblue","orange"), BOOL_axes = FALSE)
    p4 <- PCA_fun_with_Fgf4_index(Dataset = Posfai.TPM[str_detect(Posfai.TPM$Code,'C16') |
                                         str_detect(Posfai.TPM$Code,'C64'),], genes = genes.used,
                                  title = 'Posfai', values_used = c("purple","orange"), BOOL_axes = FALSE)
    p5 <- PCA_fun_with_Fgf4_index(Dataset = Wang.TPM[str_detect(Wang.TPM$Code,'8C') |
                                       str_detect(Wang.TPM$Code,'C32'),], genes = genes.used,
                                  title = 'Wang', values_used = c("darkblue","red"), BOOL_axes = FALSE)
    
    p6 <- PCA_fun_with_Fgf4_index(Dataset = Posfai.TPM[str_detect(Posfai.TPM$Code,'C16') |
                                                       str_detect(Posfai.TPM$Code,'C64'),], genes = genes.used,
                                  title = 'Posfai', values_used =  c("purple",'orange'), BOOL_axes = TRUE)
    
    
    p <- grid.arrange(p1,p2,p3,p4,p5,p6, nrow = 2)
    ggsave(paste('GGplotsFigs/EntropyGamma/PCA/PCAEP32WithIndexes.pdf', sep = ''), p, width = 24, height= 18)
    rm(genes.used)
    rm(p,p1,p2,p3,p4,p5) # PCA for the extreme time points
    
  } # PCA ETP with genes in common at 32C with indexes
  
} # PCAs of all and ETP on with genes in common at the 16C and 32C stage

{
  dummy_df <- data.frame(
    Code = factor(c("8C", "C16", "C32", "C64", "C90"), levels = c("8C", "C16", "C32", "C64", "C90")),
    PC1 = runif(n = 100, max = 100),  # Dummy values
    PC2 = runif(n = 100, max = 100)
  )
  
  p <- PCA_fun(Dataset = dummy_df, genes = c('PC1','PC2'), title= 'test', values_used = c("darkblue", "purple", "red",'orange','green'),
               BOOL_axes = FALSE, position = 'bottom')
  p
  ggsave(filename = 'GGplotsFigs/EntropyGamma/PCA/For_Legend.pdf', width = 20)
} # For legend, don't forget to but legend = 'bottom' in PCA_fun

{
  {
    threshold_bool <- TRUE
    threshold_used <- 0.4
    breaksList = seq(-1, 1, by = 0.1)
    
    ######### Allègre ######### 
    p1 <- pheatmap(Expression.Correlation_fun(Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C16'),], genes.21,
                                              threshold_bool, threshold_used),
                   cluster_rows = FALSE, cluster_cols = FALSE, main = 'Allègre 16C', silent = TRUE,
                   breaks = breaksList, color =  inferno(20),
                   fontsize = 25,
                   legend = FALSE)
    p1
    
    p2 <- pheatmap(Expression.Correlation_fun(Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C32'),], genes.21,
                                              threshold_bool, threshold_used),
                   cluster_rows = FALSE, cluster_cols = FALSE, main = 'Allègre 32C', silent = TRUE,
                   breaks = breaksList, color =  inferno(20),
                   fontsize = 25,
                   legend = FALSE)
    
    p3 <- pheatmap(Expression.Correlation_fun(Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C64'),], genes.21,
                                              threshold_bool, threshold_used),
                   cluster_rows = FALSE, cluster_cols = FALSE, main = 'Allègre 64C', silent = TRUE,
                   breaks = breaksList, color = inferno(20),
                   fontsize = 25,
                   legend = FALSE)
    
    p4 <- pheatmap(Expression.Correlation_fun(Allegre.mRNA[str_detect(Allegre.mRNA$Code,'C90'),], genes.21,
                                              threshold_bool, threshold_used),
                   cluster_rows = FALSE, cluster_cols = FALSE, main = 'Allègre 90C', silent = TRUE,
                   breaks = breaksList, color =  inferno(20),
                   fontsize = 25,
                   legend = FALSE)
    
    p <- grid.arrange(p1$gtable,p2$gtable,p3$gtable,p4$gtable,
                      nrow = 2)
    ggsave(paste('GGplotsFigs/EntropyGamma/Heatmaps/Correlations/ExpressionCorrelationAllegre.jpeg', sep = ''), p, width = 15, height = 15)
    rm(p)
    
  } # Allegre Correlations 
  
  {
    threshold_bool <- TRUE
    threshold_used <- 0.4
    breaksList = seq(-1, 1, by = 0.1)
    
    p5 <- pheatmap(Expression.Correlation_fun(Goolam.TPM[str_detect(Goolam.TPM$Code,'8C'),], genes.21,
                                              threshold_bool, threshold_used),
                   cluster_rows = FALSE, cluster_cols = FALSE, main = 'Goolam 8C', silent = TRUE,
                   breaks = breaksList, color = inferno(20),
                   fontsize = 25)
    
    p6 <- pheatmap(Expression.Correlation_fun(Goolam.TPM[str_detect(Goolam.TPM$Code,'C16'),], genes.21,
                                              threshold_bool, threshold_used),
                   cluster_rows = FALSE, cluster_cols = FALSE, main = 'Goolam 16C', silent = TRUE,
                   breaks = breaksList, color = inferno(20),
                   fontsize = 25)
    
    p7 <- pheatmap(Expression.Correlation_fun(Goolam.TPM[str_detect(Goolam.TPM$Code,'C32'),], genes.21,
                                              threshold_bool, threshold_used),
                   cluster_rows = FALSE, cluster_cols = FALSE, main = 'Goolam 32C', silent = TRUE,
                   breaks = breaksList, color = inferno(20),
                   fontsize = 25)
    
    p <- grid.arrange(p5$gtable,p6$gtable,p7$gtable,
                      nrow = 2)
    
    ggsave(paste('GGplotsFigs/EntropyGamma/Heatmaps/Correlations/ExpressionCorrelationGoolam.jpeg', sep = ''), p, width = 15, height = 15)
    rm(p)
    
  } # Goolam Correlations
  
  {
    threshold_bool <- TRUE
    threshold_used <- 0.4
    breaksList = seq(-1, 1, by = 0.1)
    
    p8 <- pheatmap(Expression.Correlation_fun(Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'8C'),], genes.21,
                                              threshold_bool, threshold_used),
                   cluster_rows = FALSE, cluster_cols = FALSE, main = 'Guo 8C', silent = TRUE,
                   breaks = breaksList, color = inferno(20),
                   fontsize = 25,
                   legend = FALSE)
    
    p9 <- pheatmap(Expression.Correlation_fun(Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'C16'),], genes.21,
                                              threshold_bool, threshold_used),
                   cluster_rows = FALSE, cluster_cols = FALSE, main = 'Guo 16C', silent = TRUE,
                   breaks = breaksList, color = inferno(20),
                   fontsize = 25,
                   legend = FALSE)
    
    p10 <- pheatmap(Expression.Correlation_fun(Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'C32'),], genes.21,
                                               threshold_bool, threshold_used),
                    cluster_rows = FALSE, cluster_cols = FALSE, main = 'Guo 32C', silent = TRUE,
                    breaks = breaksList, color = inferno(20),
                    fontsize = 25,
                    legend = FALSE)
    
    p11 <- pheatmap(Expression.Correlation_fun(Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'C64'),], genes.21,
                                               threshold_bool, threshold_used),
                    cluster_rows = FALSE, cluster_cols = FALSE, main = 'Guo 64C', silent = TRUE,
                    breaks = breaksList, color = inferno(20),
                    fontsize = 25,
                    legend = FALSE)
    
    
    p <- grid.arrange(p8$gtable,p9$gtable,p10$gtable,p11$gtable,
                      nrow = 2)
    ggsave(paste('GGplotsFigs/EntropyGamma/Heatmaps/Correlations/ExpressionCorrelationGuo.jpeg', sep = ''), p, width = 15, height = 15)
    rm(threshold_bool,threshold_used,breaksList)
    rm(p)
    
  } # Guo Correlations
  
  {
    threshold_bool <- TRUE
    threshold_used <- 0.4
    breaksList = seq(-1, 1, by = 0.1)
    
    
    p12 <- pheatmap(Expression.Correlation_fun(Posfai.TPM[str_detect(Posfai.TPM$Code,'C16'),], genes.21,
                                               threshold_bool, threshold_used),
                    cluster_rows = FALSE, cluster_cols = FALSE, main = 'Posfai 16C', silent = TRUE,
                    breaks = breaksList, color = inferno(20),
                    fontsize = 25)
    
    p13 <- pheatmap(Expression.Correlation_fun(Posfai.TPM[str_detect(Posfai.TPM$Code,'C32'),], genes.21,
                                               threshold_bool, threshold_used),
                    cluster_rows = FALSE, cluster_cols = FALSE, main = 'Posfai 32C', silent = TRUE,
                    breaks = breaksList, color = inferno(20),
                    fontsize = 25)
    
    p14 <- pheatmap(Expression.Correlation_fun(Posfai.TPM[str_detect(Posfai.TPM$Code,'C64'),], genes.21,
                                               threshold_bool, threshold_used),
                    cluster_rows = FALSE, cluster_cols = FALSE, main = 'Posfai 64C', silent = TRUE,
                    breaks = breaksList, color = inferno(20),
                    fontsize = 25)
    
    p <- grid.arrange(p12$gtable,p13$gtable,p14$gtable,
                      nrow = 2)
    ggsave(paste('GGplotsFigs/EntropyGamma/Heatmaps/Correlations/ExpressionCorrelationPosfai.jpeg', sep = ''), p, width = 15, height = 15)
    rm(threshold_bool,threshold_used,breaksList)
    rm(p)
    
  } # Posfai Correlations
  
  {
    threshold_bool <- TRUE
    threshold_used <- 0.4
    breaksList = seq(-1, 1, by = 0.1)
    
    p15 <- pheatmap(Expression.Correlation_fun(Wang.TPM[str_detect(Wang.TPM$Code,'8C'),], genes.21,
                                               threshold_bool, threshold_used),
                    cluster_rows = FALSE, cluster_cols = FALSE, main = 'Wang 8C', silent = TRUE,
                    breaks = breaksList, color = inferno(20),
                    fontsize = 25)
    
    p16 <- pheatmap(Expression.Correlation_fun(Wang.TPM[str_detect(Wang.TPM$Code,'C16'),], genes.21,
                                               threshold_bool, threshold_used),
                    cluster_rows = FALSE, cluster_cols = FALSE, main = 'Wang 16C', silent = TRUE,
                    breaks = breaksList, color = inferno(20),
                    fontsize = 25)
    
    p17 <- pheatmap(Expression.Correlation_fun(Wang.TPM[str_detect(Wang.TPM$Code,'C32'),], genes.21,
                                               threshold_bool, threshold_used),
                    cluster_rows = FALSE, cluster_cols = FALSE, main = 'Wang 32C', silent = TRUE,
                    breaks = breaksList, color = inferno(20),
                    fontsize = 25)
    
    p <- grid.arrange(p15$gtable,p16$gtable,p17$gtable,
                      nrow = 2)
    
    ggsave(paste('GGplotsFigs/EntropyGamma/Heatmaps/Correlations/ExpressionCorrelationWang.jpeg', sep = ''), p, width = 15, height = 15)
    rm(p)
    
  } # Wang Correlations
  
  {
    p_all_final <- grid.arrange(NULL,p1$gtable,p2$gtable,p3$gtable,p4$gtable,
                      p5$gtable,p6$gtable,p7$gtable,NULL,NULL,
                      p8$gtable,p9$gtable,p10$gtable,p11$gtable,NULL,
                      NULL,p12$gtable,p13$gtable,p14$gtable,NULL,
                      p15$gtable,p16$gtable,p17$gtable,NULL,NULL,
                      nrow = 5)
    ggsave(paste('GGplotsFigs/EntropyGamma/Heatmaps/Correlations/ExpressionCorrelationAll.pdf', sep = ''), p_all_final, width = 35, height = 35)
    
    p_All_Guo <- grid.arrange(p1$gtable,p2$gtable,p3$gtable,p4$gtable,
                      p8$gtable,p9$gtable,p10$gtable,p11$gtable, nrow = 2)
    
    ggsave(paste('GGplotsFigs/EntropyGamma/Heatmaps/Correlations/ExpressionCorrelationAllegreGuo.pdf', sep = ''), p_All_Guo, width = 29, height = 15)
    
    
    rm(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17)
    rm(p_all_final,p_All_Guo)
  } # All Correlations Graph
  
} # Heatmaps of Correlations  

{
  
  {
    threshold_cor <- 0.5
    genes <- genes.32
    Cellular.Stage <- c('8C','C16','C32','C64','C90')
    list.of.datasets <- list(Allegre.mRNA,Guo.mRNA.for.entropy,Posfai.mRNA,Wang.mRNA,Goolam.mRNA)
    
    list.of.df <- lapply(list.of.datasets, FUN = Barplots_fun, genes = genes.21, Cellular.Stage = Cellular.Stage, threshold = threshold_cor)
    
    
    par(mfrow = c(1,5))
    barplot(list.of.df[[1]], main=paste('Allegre'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[5]], main=paste('Goolam'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[2]], main=paste('Guo'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[3]], main=paste('Posfai'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[4]], main=paste('Wang'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    
    
    rm(list.of.datasets,list.of.df,threshold_cor,genes,Cellular.Stage)
    
    
  } # Correlation Barplots with threshold 0.5
  
  {
    threshold_cor <- 0.3
    genes <- genes.32
    Cellular.Stage <- c('8C','C16','C32','C64','C90')
    list.of.datasets <- list(Allegre.mRNA,Guo.mRNA.for.entropy,Posfai.mRNA,Wang.mRNA,Goolam.mRNA)
    
    list.of.df <- lapply(list.of.datasets, FUN = Barplots_fun, genes = genes.21, Cellular.Stage = Cellular.Stage, threshold = threshold_cor)
    
    
    par(mfrow = c(1,5))
    barplot(list.of.df[[1]], main=paste('Allegre'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[5]], main=paste('Goolam'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[2]], main=paste('Guo'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[3]], main=paste('Posfai'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[4]], main=paste('Wang'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    
    rm(list.of.datasets,list.of.df,threshold_cor,genes,Cellular.Stage)
    
  } # Correlation Barplots with threshold 0.3
  
  {
    threshold_cor <- 0.1
    genes <- genes.32
    Cellular.Stage <- c('8C','C16','C32','C64','C90')
    list.of.datasets <- list(Allegre.mRNA,Guo.mRNA.for.entropy,Posfai.mRNA,Wang.mRNA,Goolam.mRNA)
    
    list.of.df <- lapply(list.of.datasets, FUN = Barplots_fun, genes = genes.21, Cellular.Stage = Cellular.Stage, threshold = threshold_cor)
    
    
    par(mfrow = c(1,5))
    barplot(list.of.df[[1]], main=paste('Allegre'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[5]], main=paste('Goolam'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[2]], main=paste('Guo'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[3]], main=paste('Posfai'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    barplot(list.of.df[[4]], main=paste('Wang'),
            col=c("red","darkblue"),ylim = c(0,150),
            beside = TRUE)
    
    
    legend('topright',fill = c('red','darkblue'), c('Correlations','Anti-Correlations'))
    
    rm(list.of.datasets,list.of.df,threshold_cor,genes,Cellular.Stage)
  } # Correlation Barplots with threshold 0.1
  
} # Correlation Barplots


######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### 

