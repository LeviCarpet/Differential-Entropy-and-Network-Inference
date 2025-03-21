######################## CARDAMOM #########################
## In this file, we pre-process data for CARDAMOM and create the qgraphs of the inferred GRNs for the article.

######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### Setting the directory for loading data ######### ######### 
setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Entropy")
load('FinalData/InterCellularEntropy/EntropyGamma/PreProcessedData.RData') 
# Note 8C = E2.5 - 0
# Note 16C = E3.0 - 12h
# Note 32C = E3.25 - 18h
# Note 64C = E3.5 - 24h
# Note 90C = E3.75 - 30h
######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### 





######### ######### ######### ######### ######### ######### ######### #########  
######### ######### Load necessary libraries ######### ######### 
library(dplyr)  # For data manipulation
library(stringr) # str_detect function, super important
library(ggbiplot) # PCA plots
library(gridExtra) # To generate the PCAs together
library(qgraph) # To generate the regulation graphs that look the same
library(uwot) # for umaps
######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### 

######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### Genes chosen and degradation rates ######### #########  
genes.20 <- c('Code','Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6',
              'Hnf4a','Id2','Klf2', 'Klf4','Klf5','Nanog','Pdgfra','Pecam1','Pou5f1','Sox17','Sox2')

genes.21 <- c('Code','Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6',
              'Hnf4a','Id2','Klf2', 'Klf4','Klf5','Nanog','Pdgfa','Pdgfra','Pecam1','Pou5f1','Sox17','Sox2')
genes.21.with.code.and.pop <- c('Code','Population','Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6',
                                'Hnf4a','Id2','Klf2', 'Klf4','Klf5','Nanog','Pdgfa','Pdgfra','Pecam1','Pou5f1','Sox17','Sox2')

{
  degradation_rates_df_21_genes <- data.frame(GeneName=character(), 
                                              mRNA_deg = numeric(), 
                                              prot_deg = numeric(),
                                              stringsAsFactors = FALSE) 
  degradation_rates_df_21_genes[1,] <- c('Stimulus', 1, 0.2)
  degradation_rates_df_21_genes[2,] <- c('Bmp4', log(2)/2.47, log(2)/0.25) # use a protein life time of 7-16 min which is 0.25 hours at best
  degradation_rates_df_21_genes[3,] <- c('Cdx2', log(2)/2.09, log(2)/3) # THIS VALUE IS HIGHLY SPECULATIVE, BASED ON ARTICLE ABOUT HUMAN TUMOURS
  degradation_rates_df_21_genes[4,] <- c('Dab2', log(2)/7.58, log(2)/14.4) # This value is for monocytes/hepatocytes, not on mice embryonic cells
  degradation_rates_df_21_genes[5,] <- c('Esrrb', log(2)/4.56, log(2)/2) # in article by Knusden Bipartite
  degradation_rates_df_21_genes[6,] <- c('Fgf4',  log(2)/2.17, log(2)/72)
  degradation_rates_df_21_genes[7,] <- c('Fgfr2', log(2)/2.98, log(2)/35.921) 
  degradation_rates_df_21_genes[8,] <- c('Gata3', log(2)/1.7, log(2)/2.5) # based on decay of GATA3 in 6 hours and another article that says it has a half life of 2h, so median of 3 and 2 makes 2.5, and also information saying it is similar half life to GATA6 apparently from litterature
  degradation_rates_df_21_genes[9,] <- c('Gata4', log(2)/2.35, log(2)/6) # half life of GATA4 is larger than 6hours based on article
  degradation_rates_df_21_genes[10,] <- c('Gata6', log(2)/1.78, log(2)/4) # based on many articles
  degradation_rates_df_21_genes[11,] <- c('Hnf4a', log(2)/10.3, log(2)/4.5) # fairly speculative based on article Fig6C of an article looking at interaction between HNF4A and HSP90b, also done on humans
  degradation_rates_df_21_genes[12,] <- c('Id2', log(2)/1.87, log(2)/0.25) # based on article in human
  degradation_rates_df_21_genes[13,] <- c('Klf2', log(2)/1.87, log(2)/2) # based on 2 articles, but quite speculative
  degradation_rates_df_21_genes[14,] <- c('Klf4', log(2)/1.91, log(2)/1.46) # ASK Olivier, have an article about it. For now i'm using the half life when differentiated
  degradation_rates_df_21_genes[15,] <- c('Klf5', log(2)/2.92, log(2)/1.5) # Based on one article that states this verbatim that it's around 1.5h, and two others that say it has a half life of 2hs (in tumours though) 
  degradation_rates_df_21_genes[16,] <- c('Nanog', log(2)/5.21, log(2)/2.3) # Based on article that computed OCT4 and NANOG half life
  degradation_rates_df_21_genes[17,] <- c('Pdgfa', log(2)/4.52, (log(2)/5.42)/5) # demi vie fois 5, ou taux de degradation divise par 5, en gros faux que la proteine soit environ 5 fois plus stable que l'ARNm
  degradation_rates_df_21_genes[18,] <- c('Pdgfra', log(2)/3.03, log(2)/47.99) # mouse neurons
  degradation_rates_df_21_genes[19,] <- c('Pecam1', log(2)/16.7, log(2)/36) # Rough approximation using article that says it dissapears from cell surface after 3 days
  degradation_rates_df_21_genes[20,] <- c('Pou5f1', log(2)/7.4, log(2)/8) # Using 8h because MANY articles found the range >6h, and <24h (but not for WT, for which it was 8h), but the best article found around 8h. 
  degradation_rates_df_21_genes[21,] <- c('Sox17', log(2)/1.18, log(2)/3.5) # based on couple of articles and a nice figure in one of them 
  degradation_rates_df_21_genes[22,] <- c('Sox2', log(2)/1.09, log(2)/8.2)
  #View(degradation_rates_df_21_genes)
  
  setwd("C:/Users/franc/CARDAMOM_OG/cardamom/DONE/")
  write.csv(degradation_rates_df_21_genes, 'degradation_rates_df_21_genes.csv' ,row.names = FALSE)
} # degradation rates df for 21 genes in common 
######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### 


######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### Define the functions ######### ######### ######### ######### ######### #########
automate_analysis <- function(mRNA_data, levels, hours, degradation_rates, saving_name,bootstrap = FALSE,n_bootstraps = 1000) {
  
  # Ensure the inputs are data frames
  if (!is.data.frame(mRNA_data) || !is.data.frame(degradation_rates)) {
    stop("Both mRNA_data and degradation_rates must be data frames.")
  }
  
  new_data <- mRNA_data
  if (bootstrap == TRUE){
    # Step 1: Filter mRNA_data for the genes of interest and relevant time points (dynamically detect them)
    unique_timepoints <- unique(mRNA_data$Code)
    
    # Step 2: Split the mRNA_data by time points and perform bootstrapping
    bootstrapped_list <- lapply(unique_timepoints, function(tp) {
      tp_data <- mRNA_data %>% filter(Code == tp)
      sample_n(tp_data, n_bootstraps, replace = TRUE)
    })
    
    
    # Combine all bootstrapped samples
    bootstrapped_df <- bind_rows(bootstrapped_list)
    new_data <- bootstrapped_df
  }
  
  
  new_data[,2:ncol(new_data)] <- new_data[,2:ncol(new_data)]
  new_data <- new_data %>%
    mutate_if(is.numeric, as.integer)
  transposed <- as.data.frame(t(new_data))
  
  for (i in 1:length(levels)){
    #print(levels[i])
    #print(hours[i])
    transposed['Code',str_detect(transposed['Code',],levels[i])] <- hours[i]
  }
  
  transposed['Stimulus',] <- 1
  transposed['Stimulus',transposed['Code',] == '0'] <- 0
  panel_almost_real <- rbind(transposed['Code',],transposed['Stimulus',])
  panel_real <- rbind(panel_almost_real,transposed[2:(nrow(transposed)-1),])

  genes <- rownames(panel_real)[2:dim(panel_real)[1]]
  N <- dim(degradation_rates)[1]
  number.vector <- 0:(N-1)
  panel_genes <- cbind(number.vector,genes)
  Genenames <- panel_genes[1:N,2]
  
  degradation.rates_df_auto <- degradation_rates[,2:3]
  
  dir.create(paste('C:/Users/franc/CARDAMOM_OG/cardamom/',saving_name, sep = ''))
  dir.create(paste('C:/Users/franc/CARDAMOM_OG/cardamom/',saving_name,'/Data', sep = ''))
  dir.create(paste('C:/Users/franc/CARDAMOM_OG/cardamom/',saving_name,'/Rates', sep = ''))
  dir.create(paste('C:/Users/franc/CARDAMOM_OG/cardamom/',saving_name,'/cardamom', sep = ''))
  dir.create(paste('C:/Users/franc/CARDAMOM_OG/cardamom/',saving_name,'/Results', sep = ''))
  
  write.table(panel_real, file = paste('C:/Users/franc/CARDAMOM_OG/cardamom/',saving_name,'/Data/panel_real.txt',sep =''), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(panel_genes, file = paste('C:/Users/franc/CARDAMOM_OG/cardamom/',saving_name,'/Data/panel_genes.txt',sep =''), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(Genenames, file = paste('C:/Users/franc/CARDAMOM_OG/cardamom/',saving_name,'/Data/Genenames.txt',sep =''), sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(degradation.rates_df_auto, file = paste('C:/Users/franc/CARDAMOM_OG/cardamom/',saving_name,'/Rates/degradation_rates.txt',sep ='') , sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Return the final results (modify as per your requirements)
  return(invisible(NULL))
  
}
matrix.for.qgraphs_fun <- function(path,matrix_name,genes_rows_columns){
  df <- as.data.frame(read.csv(paste(path,matrix_name,sep = ''), header = FALSE))
  colnames(df) <- genes_rows_columns
  rownames(df) <- genes_rows_columns
  colnames(df)[1] <- 'Stimulus'
  rownames(df)[1] <- 'Stimulus'
  final.matrix <- df
  return(final.matrix)
}
Qgraph.MergedData  <- function(load_path,threshold_load_path,saving_path,
                               genes_modeled,genes_displayed,layout_matrix,
                               graph_title,threshold_vec,disposition_full = c(2,2),
                               nbr_of_stage_transitions_or_columns = 4){
  
  path <- paste('DONE/',load_path,'/cardamom/',threshold_load_path,'/',sep = '')
  print(path)
  
  setwd("C:/Users/franc/CARDAMOM_OG/cardamom/")
  dataset <- matrix.for.qgraphs_fun(path, 'inter.csv', genes_modeled)
  qgraph0 <- matrix.for.qgraphs_fun(path, 'inter_0.csv', genes_modeled)
  qgraph1 <- matrix.for.qgraphs_fun(path, 'inter_1.csv', genes_modeled)
  qgraph2 <- matrix.for.qgraphs_fun(path, 'inter_2.csv', genes_modeled)
  qgraph3 <- matrix.for.qgraphs_fun(path, 'inter_3.csv', genes_modeled)
  qgraph4 <- matrix.for.qgraphs_fun(path, 'inter_4.csv', genes_modeled)
  genes_displayed[1] <-  'Stimulus'
  
  dataset <- dataset[genes_displayed,genes_displayed]
  qgraph0 <- qgraph0[genes_displayed,genes_displayed]
  qgraph1 <- qgraph1[genes_displayed,genes_displayed]
  qgraph2 <- qgraph2[genes_displayed,genes_displayed]
  qgraph3 <- qgraph3[genes_displayed,genes_displayed]
  qgraph4 <- qgraph4[genes_displayed,genes_displayed]
  setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Gene Network Inference/")
  print(paste("qgraphs/MergedData/",saving_path,"_GRN.pdf", sep = ''))
  pdf(paste("qgraphs/MergedData/",saving_path,"_GRN.pdf", sep = ''))
  par(mfrow = disposition_full) 
  for (i in 1:length(threshold_vec)){
    qgraph(input = as.matrix(dataset), directed = TRUE, layout = layout_matrix,
           #title = paste(graph_title,' Starting at 8C - Thresh ', threshold_vec[i] , sep = ''),
           threshold = threshold_vec[i], cut = threshold_vec[i]*5)
  }
  dev.off()
  
  
  print(paste("qgraphs/MergedData/",saving_path,"_GRN_By_Stage.pdf",sep = ''))
  pdf(paste("qgraphs/MergedData/",saving_path,"_GRN_By_Stage.pdf",sep = ''))
  par(mfrow = c(length(threshold_vec), nbr_of_stage_transitions_or_columns))
  for (i in 1:length(threshold_vec)){
    qgraph(input = as.matrix(qgraph0), directed = TRUE, layout = layout_matrix,
           #title = paste('8C'),# - Thresh ', threshold_vec[i] , sep = ''), 
           threshold = threshold_vec[i], cut = threshold_vec[i]*5)
    mtext("8C", side = 3, line = 1, cex = 3, font = 2, adj = 0.2)  # cex controls font size, font=2 makes it bold
    qgraph(input = as.matrix(qgraph1), directed = TRUE, layout = layout_matrix,
           #title = paste('16C'),# - Thresh ', threshold_vec[i] , sep = ''), 
           threshold = threshold_vec[i], cut = threshold_vec[i]*5)
    mtext("16C", side = 3, line = 1, cex = 3, font = 2, adj = 0.2)  # cex controls font size, font=2 makes it bold
    qgraph(input = as.matrix(qgraph2), directed = TRUE, layout = layout_matrix,
           #title = paste('32C'),# - Thresh ', threshold_vec[i] , sep = ''), 
           threshold = threshold_vec[i], cut = threshold_vec[i]*5)
    mtext("32C", side = 3, line = 1, cex = 3, font = 2, adj = 0.2)  # cex controls font size, font=2 makes it bold
    qgraph(input = as.matrix(qgraph3), directed = TRUE, layout = layout_matrix,
           #title = paste('64C'),# - Thresh ', threshold_vec[i] , sep = ''), 
           threshold = threshold_vec[i], cut = threshold_vec[i]*5)
    mtext("64C", side = 3, line = 1, cex = 3, font = 2, adj = 0.2)  # cex controls font size, font=2 makes it bold
    qgraph(input = as.matrix(qgraph4), directed = TRUE, layout = layout_matrix,
           #title = paste('90C'), - Thresh ', threshold_vec[i] , sep = ''), 
           threshold = threshold_vec[i], cut = threshold_vec[i]*5)
    mtext("90C", side = 3, line = 1, cex = 3, font = 2, adj = 0.2)  # cex controls font size, font=2 makes it bold
  }
  dev.off()
}
Qgraph.Only.One.Gene  <- function(load_path, threshold_load_path, saving_path,
                                  genes_modeled, genes_displayed, gene_name, layout_matrix,
                                  graph_title, threshold_vec, disposition_full = c(2,2),
                                  nbr_of_stage_transitions_or_columns = 4){
  
  path <- paste('DONE/', load_path, '/cardamom/', threshold_load_path, '/', sep = '')
  print(path)
  
  setwd("C:/Users/franc/CARDAMOM_OG/cardamom/")
  
  # Load datasets
  dataset <- matrix.for.qgraphs_fun(path, 'inter.csv', genes_modeled)
  qgraph0 <- matrix.for.qgraphs_fun(path, 'inter_0.csv', genes_modeled)
  qgraph1 <- matrix.for.qgraphs_fun(path, 'inter_1.csv', genes_modeled)
  qgraph2 <- matrix.for.qgraphs_fun(path, 'inter_2.csv', genes_modeled)
  qgraph3 <- matrix.for.qgraphs_fun(path, 'inter_3.csv', genes_modeled)
  qgraph4 <- matrix.for.qgraphs_fun(path, 'inter_4.csv', genes_modeled)
  
  dataset <- as.matrix(dataset)
  qgraph0 <- as.matrix(qgraph0)
  qgraph1 <- as.matrix(qgraph1)
  qgraph2 <- as.matrix(qgraph2)
  qgraph3 <- as.matrix(qgraph3)
  qgraph4 <- as.matrix(qgraph4)
  genes_displayed[1] <- 'Stimulus'
  
  # Filter to keep only displayed genes
  
  dataset <- dataset[genes_displayed, genes_displayed]
  qgraph0 <- qgraph0[genes_displayed, genes_displayed]
  qgraph1 <- qgraph1[genes_displayed, genes_displayed]
  qgraph2 <- qgraph2[genes_displayed, genes_displayed]
  qgraph3 <- qgraph3[genes_displayed, genes_displayed]
  qgraph4 <- qgraph4[genes_displayed, genes_displayed]
  
  filter_interactions <- function(mat, gene_name) {
    # Ensure mat is a matrix
    if (!is.matrix(mat)) {
      stop("Error: Input is not a matrix. Check your data format.")
    }
    
    # Check if gene exists in both rows and columns
    if (!(gene_name %in% rownames(mat))) {
      stop(paste("Error: Gene", gene_name, "not found in rownames. Available genes:", paste(rownames(mat), collapse = ", ")))
    }
    if (!(gene_name %in% colnames(mat))) {
      stop(paste("Error: Gene", gene_name, "not found in colnames. Available genes:", paste(colnames(mat), collapse = ", ")))
    }
    
    # Create an empty matrix with the same dimensions
    filtered_mat <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
    
    # Keep only the interactions related to gene_name
    filtered_mat[gene_name, ] <- mat[gene_name, , drop = FALSE]  # Row interactions
    filtered_mat[, gene_name] <- mat[, gene_name, drop = FALSE]  # Column interactions
    
    return(filtered_mat)
  }
  
  
  
  # Apply filtering to keep only interactions of the 16th gene
  dataset_filtered <- filter_interactions(dataset,gene_name)
  qgraph0_filtered <- filter_interactions(qgraph0,gene_name)
  qgraph1_filtered <- filter_interactions(qgraph1,gene_name)
  qgraph2_filtered <- filter_interactions(qgraph2,gene_name)
  qgraph3_filtered <- filter_interactions(qgraph3,gene_name)
  qgraph4_filtered <- filter_interactions(qgraph4,gene_name)
  
  setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Gene Network Inference/")
  
  # **Save Merged Network**
  print(paste("qgraphs/MergedData/", saving_path, "_GRN.pdf", sep = ''))
  pdf(paste("qgraphs/MergedData/", saving_path, "_GRN.pdf", sep = ''))
  par(mfrow = disposition_full)
  
  for (i in seq_along(threshold_vec)) {
    qgraph(input = as.matrix(dataset_filtered), directed = TRUE, layout = layout_matrix,
           threshold = threshold_vec[i], cut = threshold_vec[i] * 5)
  }
  dev.off()
  
  # **Save Network by Stage**
  print(paste("qgraphs/MergedData/", saving_path, "_GRN_By_Stage.pdf", sep = ''))
  pdf(paste("qgraphs/MergedData/", saving_path, "_GRN_By_Stage.pdf", sep = ''))
  par(mfrow = c(length(threshold_vec), nbr_of_stage_transitions_or_columns))
  
  for (i in seq_along(threshold_vec)) {
    qgraph(input = as.matrix(qgraph0_filtered), directed = TRUE, layout = layout_matrix,
          threshold = threshold_vec[i], cut = threshold_vec[i] * 5)
    mtext("8C", side = 3, line = 2, cex = 3, font = 2)  # cex controls font size, font=2 makes it bold
    qgraph(input = as.matrix(qgraph1_filtered), directed = TRUE, layout = layout_matrix,
           threshold = threshold_vec[i], cut = threshold_vec[i] * 5)
    mtext("16C", side = 3, line = 2, cex = 3, font = 2)  # cex controls font size, font=2 makes it bold
    qgraph(input = as.matrix(qgraph2_filtered), directed = TRUE, layout = layout_matrix,
           threshold = threshold_vec[i], cut = threshold_vec[i] * 5)
    mtext("32C", side = 3, line = 2, cex = 3, font = 2)  # cex controls font size, font=2 makes it bold
    qgraph(input = as.matrix(qgraph3_filtered), directed = TRUE, layout = layout_matrix,
           threshold = threshold_vec[i], cut = threshold_vec[i] * 5)
    mtext("64C", side = 3, line = 2, cex = 3, font = 2)  # cex controls font size, font=2 makes it bold
    qgraph(input = as.matrix(qgraph4_filtered), directed = TRUE, layout = layout_matrix,
           threshold = threshold_vec[i], cut = threshold_vec[i] * 5)
    mtext("90C", side = 3, line = 2, cex = 3, font = 2)  # cex controls font size, font=2 makes it bold
  }
  dev.off()
}
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 


######### ######### ######### ######### ######### ######### ######### ######### ######### #########  
######### ######### Pre-Processing of datasets and creation of Full Datasets ######### ######### 
{ # Merging multiple datasets with a new column indicating their source
  # Processing Allegre dataset
  Allegre.data <- cbind(
    Allegre.mRNA[!str_detect(Allegre.mRNA$Code, 'C16'), genes.21], # Filter rows excluding 'C16'
    Dataset = 'Allegre' # Add 'Dataset' column
  )
  
  # Processing Goolam dataset
  Goolam.data <- cbind(
    Goolam.mRNA[
      str_detect(Goolam.mRNA[, 'Code'], '8C') | 
        str_detect(Goolam.mRNA[, 'Code'], 'C16') | 
        str_detect(Goolam.mRNA[, 'Code'], 'C32'),
      genes.21
    ],
    Dataset = 'Goolam'
  )
  
  # Processing Guo dataset
  Guo.data <- cbind(
    Guo.mRNA.for.CARDAMOM[!(
      str_detect(Guo.mRNA.for.CARDAMOM[, 'Code'], '2C') | 
        str_detect(Guo.mRNA.for.CARDAMOM[, 'Code'], '4C')
    ), genes.21],
    Dataset = 'Guo'
  )
  
  # Processing Posfai dataset
  Posfai.data <- cbind(
    Posfai.mRNA[, genes.21], # Include all rows
    Dataset = 'Posfai'
  )
  
  # Processing Wang dataset
  Wang.data  <- cbind(
    Wang.mRNA[
      str_detect(Wang.mRNA[, 'Code'], '8C') | 
        str_detect(Wang.mRNA[, 'Code'], 'C16') | 
        str_detect(Wang.mRNA[, 'Code'], 'C32'),
      genes.21
    ],
    Dataset = 'Wang'
  )
  
  # Combine all processed datasets into one
  Five_Datasets_21_genes <- rbind(Allegre.data,Goolam.data,Guo.data,Posfai.data,Wang.data)
  
  # Cleanup intermediate variables
  rm(Allegre.data,Goolam.data,Guo.data,Posfai.data,Wang.data)
  
} # Creation of Five_Datasets_21_genes
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 




######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### PCAs to check that the datasets are indeed comparable ######### ######### 
setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco")
{
  PCA.8C <- prcomp(Five_Datasets_21_genes[,genes.21[2:length(genes.21)]], center = TRUE, scale = TRUE)
  
  
  p1 <- ggbiplot(PCA.8C, var.axes = FALSE, varname.size = 8) +
    geom_point(aes(color = factor(x = Five_Datasets_21_genes$Code, levels = c('8C','C16','C32','C64','C90'), labels = c('8C','16C','32C','64C','90C'))), size = 4) +
    scale_colour_manual(name = "Cellular Stage", values = c('blue',"purple", "red","orange",'green') )+
    ggtitle('PCA Colored By Cellular Stage') +
    theme_minimal() +
    theme(legend.position = "bottom", legend.title=element_text(size=22),legend.text=element_text(size=22),
          plot.title = element_text(hjust = 0.5, size = 35),
          axis.title = element_text(size = 25), axis.text = element_text(size = 25)) +
    guides(color = guide_legend(override.aes = list(size = 6)))  # Increase legend point size
  
  p2 <- ggbiplot(PCA.8C, var.axes = FALSE, varname.size = 8) +
    geom_point(aes(color = factor(x = Five_Datasets_21_genes$Dataset, levels = c('Allegre','Goolam','Guo','Posfai','Wang'), labels = c('Allegre','Goolam','Guo','Posfai','Wang'))), size = 4) +
    ggtitle('PCA Colored by Dataset') +
    theme_minimal() +
    scale_colour_manual(name = "Dataset", values = c('darkred','blue','darkgreen','orange','purple'))+
    theme(legend.position = "bottom", legend.title=element_text(size=22),legend.text=element_text(size=22),
          plot.title = element_text(hjust = 0.5, size = 35),
          axis.title = element_text(size = 25), axis.text = element_text(size = 25)) +
    guides(color = guide_legend(override.aes = list(size = 6)))  # Increase legend point size
  
  p <- grid.arrange(p1,p2, nrow = 1)
  p
  ggsave(paste('Gene Network Inference/PCA/PCA_comparable_20_genes.pdf', sep = ''), p, width = 20, height = 11)
  rm(p,p1,p2,PCA.8C)
} # PCAs to check that the 5 datasets are comparable 20 genes

{   
  map <- umap2(Five_Datasets_21_genes[,genes.21[2:length(genes.21)]])
  
  # `umap2` might return a matrix, so convert it to a data frame for ggplot compatibility.
  map_df <- as.data.frame(map)
  
  # Rename columns for easy plotting
  colnames(map_df) <- c("UMAP1", "UMAP2")
  
  # Add the desired color column from `Five_Datasets_21_genes`, replace "your_column_name" with the actual column name.
  map_df$Code <- Five_Datasets_21_genes$Code
  map_df$Dataset <- Five_Datasets_21_genes$Dataset
  
  # Plotting
  p1 <- ggplot(map_df, aes(x = UMAP1, y = UMAP2, color = Code)) +
    geom_point() +
    scale_color_manual(values = c('blue',"purple", "red","orange",'green')) + 
    theme_minimal() +
    theme(legend.position = "bottom", legend.title=element_text(size=20),legend.text=element_text(size=20),
          plot.title = element_text(hjust = 0.5, size = 35),
          axis.title = element_text(size = 25), axis.text = element_text(size = 25)) +
    guides(color = guide_legend(override.aes = list(size = 6))) + # Increase legend point size
    labs(title = "UMAP Plot Colored by Cellular Stage", color = "Cellular Stage")
  
  p2 <- ggplot(map_df, aes(x = UMAP1, y = UMAP2, color = Dataset)) +
    geom_point() +
    scale_color_manual(values = c('darkred','blue','darkgreen','orange','purple')) + 
    
    theme_minimal() +
    theme(legend.position = "bottom", legend.title=element_text(size=20),legend.text=element_text(size=20),
          plot.title = element_text(hjust = 0.5, size = 35),
          axis.title = element_text(size = 25), axis.text = element_text(size = 25)) +
    guides(color = guide_legend(override.aes = list(size = 6))) + # Increase legend point size
    labs(title = "UMAP Plot Colored by Dataset", color = "Dataset")
  
  p <- grid.arrange(p1,p2,nrow = 1)
  ggsave(paste('Gene Network Inference/UMAP/UMAP_comparable_20_genes.pdf', sep = ''), p, width = 20, height = 11)
  rm(p1,p2,p,map_df,map)
} # UMAP to check that the 5 datasets are comparable 20 genes
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 





######### ######### ######### ######### ######### ######### ######### ######### 
######### Raw datasets, genes in common at the 32C stage for CARDAMOM ######### 
{
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.21],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = degradation_rates_df_21_genes,
                    saving_name = '5_Datasets_21_genes',
                    bootstrap = FALSE, n_bootstraps = 0)
} # Pre-processing into CARDAMOM 5 Datasets 21 genes

{
  genes.main.components <- c('Code','Fgf4','Fgfr2','Gata6','Nanog')
  deg_df_main_components <- degradation_rates_df_21_genes[c(1,6,7,10,16),]
  deg_df_main_components
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components,
                    saving_name = '5_Datasets_maincomponents',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components,deg_df_main_components)
} # Pre-processing into CARDAMOM 5 Datasets Main components

{
  genes.main.components.with.Bmp4 <- c('Code','Bmp4','Fgf4','Fgfr2','Gata6','Nanog')
  deg_df_main_components_with_Bmp4 <- degradation_rates_df_21_genes[c(1,2,6,7,10,16),]
  deg_df_main_components_with_Bmp4
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Bmp4],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Bmp4,
                    saving_name = '5_Datasets_maincomponents_with_Bmp4',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Bmp4,deg_df_main_components_with_Bmp4)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Bmp4

{
  genes.main.components.with.Cdx2 <- c('Code','Cdx2','Fgf4','Fgfr2','Gata6','Nanog')
  deg_df_main_components_with_Cdx2 <- degradation_rates_df_21_genes[c(1,3,6,7,10,16),]
  deg_df_main_components_with_Cdx2
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Cdx2],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Cdx2,
                    saving_name = '5_Datasets_maincomponents_with_Cdx2',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Cdx2,deg_df_main_components_with_Cdx2)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Cdx2

{
  genes.main.components.with.Dab2 <- c('Code','Dab2','Fgf4','Fgfr2','Gata6','Nanog')
  deg_df_main_components_with_Dab2 <- degradation_rates_df_21_genes[c(1,4,6,7,10,16),]
  deg_df_main_components_with_Dab2
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Dab2],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Dab2,
                    saving_name = '5_Datasets_maincomponents_with_Dab2',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Dab2,deg_df_main_components_with_Dab2)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Dab2

{
  genes.main.components.with.Esrrb <- c('Code','Esrrb','Fgf4','Fgfr2','Gata6','Nanog')
  deg_df_main_components_with_Esrrb <- degradation_rates_df_21_genes[c(1,5,6,7,10,16),]
  deg_df_main_components_with_Esrrb
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Esrrb],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Esrrb,
                    saving_name = '5_Datasets_maincomponents_with_Esrrb',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Esrrb,deg_df_main_components_with_Esrrb)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Esrrb

{
  genes.main.components.with.Gata3 <- c('Code','Fgf4','Fgfr2','Gata3','Gata6','Nanog')
  deg_df_main_components_with_Gata3 <- degradation_rates_df_21_genes[c(1,6,7,8,10,16),]
  deg_df_main_components_with_Gata3
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Gata3],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Gata3,
                    saving_name = '5_Datasets_maincomponents_with_Gata3',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Gata3,deg_df_main_components_with_Gata3)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Gata3

{
  genes.main.components.with.Gata4 <- c('Code','Fgf4','Fgfr2','Gata4','Gata6','Nanog')
  deg_df_main_components_with_Gata4 <- degradation_rates_df_21_genes[c(1,6,7,9,10,16),]
  deg_df_main_components_with_Gata4
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Gata4],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Gata4,
                    saving_name = '5_Datasets_maincomponents_with_Gata4',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(deg_df_main_components_with_Gata4,deg_df_main_components_with_Gata4)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Gata4

{
  genes.main.components.with.Hnf4a <- c('Code','Fgf4','Fgfr2','Gata6','Hnf4a','Nanog')
  deg_df_main_components_with_Hnf4a <- degradation_rates_df_21_genes[c(1,6,7,10,11,16),]
  deg_df_main_components_with_Hnf4a
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Hnf4a],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Hnf4a,
                    saving_name = '5_Datasets_maincomponents_with_Hnf4a',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Hnf4a,deg_df_main_components_with_Hnf4a)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Hnf4a

{
  genes.main.components.with.Id2 <- c('Code','Fgf4','Fgfr2','Gata6','Id2','Nanog')
  deg_df_main_components_with_Id2 <- degradation_rates_df_21_genes[c(1,6,7,10,12,16),]
  deg_df_main_components_with_Id2
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Id2],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Id2,
                    saving_name = '5_Datasets_maincomponents_with_Id2',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Id2,deg_df_main_components_with_Id2)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Id2

{
  genes.main.components.with.Klf2 <- c('Code','Fgf4','Fgfr2','Gata6','Klf2','Nanog')
  deg_df_main_components_with_Klf2 <- degradation_rates_df_21_genes[c(1,6,7,10,13,16),]
  deg_df_main_components_with_Klf2
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Klf2],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Klf2,
                    saving_name = '5_Datasets_maincomponents_with_Klf2',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Klf2,deg_df_main_components_with_Klf2)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Klf2

{
  genes.main.components.with.Klf4 <- c('Code','Fgf4','Fgfr2','Gata6','Klf4','Nanog')
  deg_df_main_components_with_Klf4 <- degradation_rates_df_21_genes[c(1,6,7,10,14,16),]
  deg_df_main_components_with_Klf4
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Klf4],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Klf4,
                    saving_name = '5_Datasets_maincomponents_with_Klf4',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Klf4,deg_df_main_components_with_Klf4)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Klf4

{
  genes.main.components.with.Klf5 <- c('Code','Fgf4','Fgfr2','Gata6','Klf5','Nanog')
  deg_df_main_components_with_Klf5 <- degradation_rates_df_21_genes[c(1,6,7,10,15,16),]
  deg_df_main_components_with_Klf5
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Klf5],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Klf5,
                    saving_name = '5_Datasets_maincomponents_with_Klf5',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Klf5,deg_df_main_components_with_Klf5)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Klf5

{
  genes.main.components.with.Pdgfa <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pdgfa')
  deg_df_main_components_with_Pdgfa <- degradation_rates_df_21_genes[c(1,6,7,10,16,17),]
  deg_df_main_components_with_Pdgfa
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Pdgfa],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Pdgfa,
                    saving_name = '5_Datasets_maincomponents_with_Pdgfa',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Pdgfa,deg_df_main_components_with_Pdgfa)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Pdgfa

{
  genes.main.components.with.Pdgfra <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pdgfra')
  deg_df_main_components_with_Pdgfra <- degradation_rates_df_21_genes[c(1,6,7,10,16,18),]
  deg_df_main_components_with_Pdgfra
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Pdgfra],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Pdgfra,
                    saving_name = '5_Datasets_maincomponents_with_Pdgfra',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Pdgfra,deg_df_main_components_with_Pdgfra)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Pdgfra

{
  genes.main.components.with.Pecam1 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pecam1')
  deg_df_main_components_with_Pecam1 <- degradation_rates_df_21_genes[c(1,6,7,10,16,19),]
  deg_df_main_components_with_Pecam1
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Pecam1],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Pecam1,
                    saving_name = '5_Datasets_maincomponents_with_Pecam1',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Pecam1,deg_df_main_components_with_Pecam1)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Pecam1

{
  genes.main.components.with.Pou5f1 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pou5f1')
  deg_df_main_components_with_Pou5f1 <- degradation_rates_df_21_genes[c(1,6,7,10,16,20),]
  deg_df_main_components_with_Pou5f1
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Pou5f1],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Pou5f1,
                    saving_name = '5_Datasets_maincomponents_with_Pou5f1',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Pou5f1,deg_df_main_components_with_Pou5f1)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Pou5f1

{
  genes.main.components.with.Sox17 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Sox17')
  deg_df_main_components_with_Sox17 <- degradation_rates_df_21_genes[c(1,6,7,10,16,21),]
  deg_df_main_components_with_Sox17
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Sox17],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Sox17,
                    saving_name = '5_Datasets_maincomponents_with_Sox17',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Sox17,deg_df_main_components_with_Sox17)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Sox17

{
  genes.main.components.with.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Sox2')
  deg_df_main_components_with_Sox2 <- degradation_rates_df_21_genes[c(1,6,7,10,16,22),]
  deg_df_main_components_with_Sox2
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Sox2],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Sox2,
                    saving_name = '5_Datasets_maincomponents_with_Sox2',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Sox2,deg_df_main_components_with_Sox2)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Sox2

{
  genes.main.components.with.Bmp4.Pecam1 <- c('Code','Bmp4','Fgf4','Fgfr2','Gata6','Nanog','Pecam1')
  deg_df_main_components_with_Bmp4_Pecam1 <- degradation_rates_df_21_genes[c(1,2,6,7,10,16,19),]
  deg_df_main_components_with_Bmp4_Pecam1
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Bmp4.Pecam1],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Bmp4_Pecam1,
                    saving_name = '5_Datasets_maincomponents_with_Bmp4_Pecam1',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Bmp4.Pecam1,deg_df_main_components_with_Bmp4_Pecam1)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Bmp4-Pecam1

{
  genes.main.components.with.Pecam1.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pecam1','Sox2')
  deg_df_main_components_with_Pecam1_Sox2 <- degradation_rates_df_21_genes[c(1,6,7,10,16,19,22),]
  deg_df_main_components_with_Pecam1_Sox2
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Pecam1.Sox2],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Pecam1_Sox2,
                    saving_name = '5_Datasets_maincomponents_with_Pecam1_Sox2',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Pecam1.Sox2,deg_df_main_components_with_Pecam1_Sox2)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Pecam1-Sox2

{
  genes.main.components.with.Hnf4a.Pecam1.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Hnf4a','Nanog','Pecam1','Sox2')
  deg_df_main_components_with_Hnf4a_Pecam1_Sox2 <- degradation_rates_df_21_genes[c(1,6,7,10,11,16,19,22),]
  deg_df_main_components_with_Hnf4a_Pecam1_Sox2
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Hnf4a.Pecam1.Sox2],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Hnf4a_Pecam1_Sox2,
                    saving_name = '5_Datasets_maincomponents_with_Hnf4a_Pecam1_Sox2',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Hnf4a.Pecam1.Sox2,deg_df_main_components_with_Hnf4a_Pecam1_Sox2)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Hnf4a-Pecam1-Sox2

{
  genes.main.components.with.Klf2.Pecam1.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Klf2','Nanog','Pecam1','Sox2')
  deg_df_main_components_with_Klf2_Pecam1_Sox2 <- degradation_rates_df_21_genes[c(1,6,7,10,13,16,19,22),]
  deg_df_main_components_with_Klf2_Pecam1_Sox2
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Klf2.Pecam1.Sox2],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Klf2_Pecam1_Sox2,
                    saving_name = '5_Datasets_maincomponents_with_Klf2_Pecam1_Sox2',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Klf2.Pecam1.Sox2,deg_df_main_components_with_Klf2_Pecam1_Sox2)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Klf2-Pecam1-Sox2

{
  genes.main.components.with.Hnf4a.Klf2.Pecam1.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Hnf4a','Klf2','Nanog','Pecam1','Sox2')
  deg_df_main_components_with_Hnf4a_Klf2_Pecam1_Sox2 <- degradation_rates_df_21_genes[c(1,6,7,10,11,13,16,19,22),]
  deg_df_main_components_with_Hnf4a_Klf2_Pecam1_Sox2
  
  automate_analysis(mRNA_data = Five_Datasets_21_genes[,genes.main.components.with.Hnf4a.Klf2.Pecam1.Sox2],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Hnf4a_Klf2_Pecam1_Sox2,
                    saving_name = '5_Datasets_maincomponents_with_Hnf4a_Klf2_Pecam1_Sox2',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Hnf4a.Klf2.Pecam1.Sox2,deg_df_main_components_with_Hnf4a_Klf2_Pecam1_Sox2)
} # Pre-processing into CARDAMOM 5 Datasets Main components with Hnf4a-Klf2-Pecam1-Sox2
######### #########  ######### #########  ######### ######### ######### ######### 
######### #########  ######### #########  ######### ######### ######### ######### 



######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### Generate QGraphs ######### ######### #########
{
  {
    layout_maincomponents <- matrix(c(
        # X-coordinates:
        -1, #Stimulus
        0, #Fgf4 
        0, #Fgfr2
        -1,#Gata6
        1, #Nanog
        # Y-coordinates:
        1.5,#Stimulus
        2,#Fgf4 
        0,#Fgfr2
        1,#Gata6
        1 #Nanog
      ), ncol=2, byrow=FALSE) # for main components
      
      
      layout_5genes_BFgf4 <- matrix(c(
        # X-coordinates:
        -1, #Stimulus
        1, #Helper Gene
        0, #Fgf4 
        0, #Fgfr2
        -1,#Gata6
        1, #Nanog
        # Y-coordinates:
        1.5,#Stimulus
        2,#Helper Gene
        2,#Fgf4 
        0,#Fgfr2
        1,#Gata6
        1 #Nanog
      ), ncol=2, byrow=FALSE) # for MC + Helper Gene before Fgf4
      
      layout_5genes_AFgfr2BGata6 <- matrix(c(
        # X-coordinates:
        -1, #Stimulus
        0, #Fgf4 
        0, #Fgfr2
        1, #Helper Gene
        -1,#Gata6
        1, #Nanog
        # Y-coordinates:
        1.5,#Stimulus
        2,#Fgf4 
        0,#Fgfr2
        2,#Helper Gene
        1,#Gata6
        1 #Nanog
      ), ncol=2, byrow=FALSE) # for MC + Helper Gene after Fgf4/Fgfr2 before Gata6
      
      
      layout_5genes_AGata6BNanog <- matrix(c(
        # X-coordinates:
        -1, #Stimulus
        0, #Fgf4 
        0, #Fgfr2
        -1,#Gata6
        1, #Helper Gene
        1, #Nanog
        # Y-coordinates:
        1.5,#Stimulus
        2,#Fgf4 
        0,#Fgfr2
        1,#Gata6
        2,#Helper Gene
        1 #Nanog
      ), ncol=2, byrow=FALSE) # for MC + Helper Gene after Gata6 before Nanog
      
      
      layout_5genes_ANanog <- matrix(c(
        # X-coordinates:
        -1, #Stimulus
        0, #Fgf4 
        0, #Fgfr2
        -1,#Gata6
        1, #Nanog
        1, #Helper Gene
        # Y-coordinates:
        1.5,#Stimulus
        2,#Fgf4 
        0,#Fgfr2
        1,#Gata6
        1,#Nanog
        2 #Helper Gene
      ), ncol=2, byrow=FALSE) # for MC + Helper Gene after Nanog
      
      
      layout_6genes_BGata6ANanog <- matrix(c(
        # X-coordinates:
        -1, #Stimulus
        1, #Helper Gene 1
        0, #Fgf4 
        0, #Fgfr2
        -1,#Gata6
        1, #Nanog
        0.5, #Helper Gene 2
        # Y-coordinates:
        1.5,#Stimulus
        2,#Helper Gene 1
        2,#Fgf4 
        0,#Fgfr2
        1,#Gata6
        1,#Nanog
        1.5#Helper Gene 2
      ), ncol=2, byrow=FALSE) # for MC + Helper Gene 1 and 2 after Nanog
      
      layout_6genes_ANanog <- matrix(c(
        # X-coordinates:
        -1, #Stimulus
        0, #Fgf4 
        0, #Fgfr2
        -1,#Gata6
        1, #Nanog
        1, #Helper Gene 1
        0.5, #Helper Gene 2
        # Y-coordinates:
        1.5,#Stimulus
        2,#Fgf4 
        0,#Fgfr2
        1,#Gata6
        1,#Nanog
        2,#Helper Gene 1
        1.5#Helper Gene 2
      ), ncol=2, byrow=FALSE) # for MC + Helper Gene 1 and 2 after Nanog
      
      
      layout_7genes_AGata6BNanog <- matrix(c(
        # X-coordinates:
        -1, #Stimulus
        0, #Fgf4 
        0, #Fgfr2
        -1,#Gata6
        0.5, #Helper Gene 1
        1, #Nanog
        1, #Helper Gene 2
        0.5, #Helper Gene 3
        # Y-coordinates:
        1.5,#Stimulus
        2,#Fgf4 
        0,#Fgfr2
        1,#Gata6
        2.2,#Helper Gene 1
        1, #Nanog
        2,#Helper Gene 2 
        1.5#Helper Gene 3
      ), ncol=2, byrow=FALSE) # for MC + Helper Gene 1, 2, 3 Hnf4a-Pecam1-Sox2
      
      layout_8genes_AGata6BNanog <- matrix(c(
        # X-coordinates:
        -1, #Stimulus
        0, #Fgf4 
        0, #Fgfr2
        -1,#Gata6
        0.5, #Helper Gene 1
        0.8, #Helper Gene 2
        1, #Nanog
        1, #Helper Gene 3
        0.5, #Helper Gene 4
        # Y-coordinates:
        1.5,#Stimulus
        2,#Fgf4 
        0,#Fgfr2
        1,#Gata6
        2.2,#Helper Gene 1
        2.5, #Helper Gene 2
        1, #Nanog
        2,#Helper Gene 3
        1.5#Helper Gene 4
      ), ncol=2, byrow=FALSE) # for MC + Helper Gene 1, 2, 3 Hnf4a-Klf2-Pecam1-Sox2
      
      # Define a layout matrix for 21 nodes
      layout_matrix_20_genes <- matrix(nrow = 21, ncol = 2)
      
      # Excentered node (Node 1)
      layout_matrix_20_genes[1, ] <- c(-1.5, 1.2)  # Positioned separately
      
      # Centered square (Nodes 6, 7, 10, 16)
      layout_matrix_20_genes[6, ]  <- c(-0.5,  0.5)
      layout_matrix_20_genes[7, ]  <- c( 0.5,  0.5)
      layout_matrix_20_genes[10, ] <- c(-0.5, -0.5)
      layout_matrix_20_genes[16, ] <- c( 0.5, -0.5)
      
      # Arrange remaining nodes around the square
      angles <- seq(0, 2 * pi, length.out = 17)[-1]  # Generate 16 angles in a circle
      radius <- 1.5  # Radius of the surrounding nodes
      
      counter <- 1
      for (i in 2:21) {
        if (i %in% c(6, 7, 10, 16)) next  # Skip square nodes
        layout_matrix_20_genes[i, ] <- c(radius * cos(angles[counter]), radius * sin(angles[counter]))
        counter <- counter + 1
      }
      
    } # Create the layout matrices
    
  {
      genes.main.components <- c('Code','Fgf4','Fgfr2','Gata6','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents',
                        genes_modeled = genes.main.components, genes_displayed = genes.main.components, layout_matrix = layout_maincomponents,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components)
      
    } # 5 Datasets MC
    
  {
      genes.main.components.Bmp4 <- c('Code','Bmp4','Fgf4','Fgfr2','Gata6','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Bmp4', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Bmp4',
                        genes_modeled = genes.main.components.Bmp4, genes_displayed = genes.main.components.Bmp4, layout_matrix = layout_5genes_BFgf4,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Bmp4)
      
    } # 5 Datasets MC + Bmp4
    
  {
      genes.main.components.Cdx2 <- c('Code','Cdx2','Fgf4','Fgfr2','Gata6','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Cdx2', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Cdx2',
                        genes_modeled = genes.main.components.Cdx2, genes_displayed = genes.main.components.Cdx2, layout_matrix = layout_5genes_BFgf4,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Cdx2)
      
    } # 5 Datasets MC + Cdx2
    
  {
      genes.main.components.Dab2 <- c('Code','Dab2','Fgf4','Fgfr2','Gata6','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Dab2', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Dab2',
                        genes_modeled = genes.main.components.Dab2, genes_displayed = genes.main.components.Dab2, layout_matrix = layout_5genes_BFgf4,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Dab2)
      
    } # 5 Datasets MC + Dab2
    
  {
      genes.main.components.Esrrb <- c('Code','Esrrb','Fgf4','Fgfr2','Gata6','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Esrrb', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Esrrb',
                        genes_modeled = genes.main.components.Esrrb, genes_displayed = genes.main.components.Esrrb, layout_matrix = layout_5genes_BFgf4,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Esrrb)
      
    } # 5 Datasets MC + Esrrb
    
  {
      genes.main.components.Gata3 <- c('Code','Fgf4','Fgfr2','Gata3','Gata6','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Gata3', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Gata3',
                        genes_modeled = genes.main.components.Gata3, genes_displayed = genes.main.components.Gata3, layout_matrix = layout_5genes_AFgfr2BGata6,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Gata3)
      
    } # 5 Datasets MC + Gata3
    
  {
      genes.main.components.Gata4 <- c('Code','Fgf4','Fgfr2','Gata4','Gata6','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Gata4', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Gata4',
                        genes_modeled = genes.main.components.Gata4, genes_displayed = genes.main.components.Gata4, layout_matrix = layout_5genes_AFgfr2BGata6,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Gata4)
      
    } # 5 Datasets MC + Gata4
    
  {
      genes.main.components.Hnf4a <- c('Code','Fgf4','Fgfr2','Gata6','Hnf4a','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Hnf4a', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Hnf4a',
                        genes_modeled = genes.main.components.Hnf4a, genes_displayed = genes.main.components.Hnf4a, layout_matrix = layout_5genes_AGata6BNanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Hnf4a)
      
    } # 5 Datasets MC + Hnf4a
    
  {
      genes.main.components.Id2 <- c('Code','Fgf4','Fgfr2','Gata6','Id2','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Id2', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Id2',
                        genes_modeled = genes.main.components.Id2, genes_displayed = genes.main.components.Id2, layout_matrix = layout_5genes_AGata6BNanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Id2)
      
    } # 5 Datasets MC + Id2
    
  {
      genes.main.components.Klf2 <- c('Code','Fgf4','Fgfr2','Gata6','Klf2','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Klf2', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Klf2',
                        genes_modeled = genes.main.components.Klf2, genes_displayed = genes.main.components.Klf2, layout_matrix = layout_5genes_AGata6BNanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Klf2)
      
    } # 5 Datasets MC + Klf2
    
  {
      genes.main.components.Klf4 <- c('Code','Fgf4','Fgfr2','Gata6','Klf4','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Klf4', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Klf4',
                        genes_modeled = genes.main.components.Klf4, genes_displayed = genes.main.components.Klf4, layout_matrix = layout_5genes_AGata6BNanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Klf4)
      
    } # 5 Datasets MC + Klf4
    
  {
      genes.main.components.Klf5 <- c('Code','Fgf4','Fgfr2','Gata6','Klf5','Nanog')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Klf5', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Klf5',
                        genes_modeled = genes.main.components.Klf5, genes_displayed = genes.main.components.Klf5, layout_matrix = layout_5genes_AGata6BNanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Klf5)
      
    } # 5 Datasets MC + Klf5
  
  {
    genes.main.components.Pdgfa <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pdgfa')
    Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Pdgfa', threshold_load_path = 'threshold_0.0_kanto_0.4',
                      saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Pdgfa',
                      genes_modeled = genes.main.components.Pdgfa, genes_displayed = genes.main.components.Pdgfa, layout_matrix = layout_5genes_ANanog,
                      graph_title = '', threshold_vec = c(0.00001),
                      disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
    rm(genes.main.components.Pdgfa)
    
  } # 5 Datasets MC + Pdgfa
    
  {
      genes.main.components.Pdgfra <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pdgfra')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Pdgfra', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Pdgfra',
                        genes_modeled = genes.main.components.Pdgfra, genes_displayed = genes.main.components.Pdgfra, layout_matrix = layout_5genes_ANanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Pdgfra)
      
    } # 5 Datasets MC + Pdgfra
    
  {
      genes.main.components.Pecam1 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pecam1')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Pecam1', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Pecam1',
                        genes_modeled = genes.main.components.Pecam1, genes_displayed = genes.main.components.Pecam1, layout_matrix = layout_5genes_ANanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Pecam1)
      
    } # 5 Datasets MC + Pecam1
    
  {
      genes.main.components.Pou5f1 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pou5f1')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Pou5f1', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Pou5f1',
                        genes_modeled = genes.main.components.Pou5f1, genes_displayed = genes.main.components.Pou5f1, layout_matrix = layout_5genes_ANanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Pou5f1)
      
    } # 5 Datasets MC + Pou5f1
    
  {
      genes.main.components.Sox17 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Sox17')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Sox17', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Sox17',
                        genes_modeled = genes.main.components.Sox17, genes_displayed = genes.main.components.Sox17, layout_matrix = layout_5genes_ANanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Sox17)
      
    } # 5 Datasets MC + Sox17
    
  {
      genes.main.components.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Sox2')
      Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Sox2', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Sox2',
                        genes_modeled = genes.main.components.Sox2, genes_displayed = genes.main.components.Sox2, layout_matrix = layout_5genes_ANanog,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Sox2)
      
    } # 5 Datasets MC + Sox2

  
  {
    genes.main.components.Bmp4.Pecam1 <- c('Code','Bmp4','Fgf4','Fgfr2','Gata6','Nanog','Pecam1')
    Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Bmp4_Pecam1', threshold_load_path = 'threshold_0.0_kanto_0.4',
                      saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Bmp4_Pecam1',
                      genes_modeled = genes.main.components.Bmp4.Pecam1, genes_displayed = genes.main.components.Bmp4.Pecam1, layout_matrix = layout_6genes_BGata6ANanog,
                      graph_title = '', threshold_vec = c(0.00001),
                      disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
    rm(genes.main.components.Bmp4.Pecam1)
    
  } # 5 Datasets MC + Bmp4 + Pecam1 

  {
    genes.main.components.Pecam1.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Nanog','Pecam1','Sox2')
    Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Pecam1_Sox2', threshold_load_path = 'threshold_0.0_kanto_0.4',
                      saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Pecam1_Sox2',
                      genes_modeled = genes.main.components.Pecam1.Sox2, genes_displayed = genes.main.components.Pecam1.Sox2, layout_matrix = layout_6genes_ANanog,
                      graph_title = '', threshold_vec = c(0.00001),
                      disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
    rm(genes.main.components.Pecam1.Sox2)
    
  } # 5 Datasets MC + Pecam1 + Sox2 
  
  
  {
    genes.main.components.Hnf4a.Pecam1.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Hnf4a','Nanog','Pecam1','Sox2')
    Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Hnf4a_Pecam1_Sox2', threshold_load_path = 'threshold_0.0_kanto_0.4',
                      saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Hnf4a_Pecam1_Sox2',
                      genes_modeled = genes.main.components.Hnf4a.Pecam1.Sox2, genes_displayed = genes.main.components.Hnf4a.Pecam1.Sox2, layout_matrix = layout_7genes_AGata6BNanog,
                      graph_title = '', threshold_vec = c(0.00001),
                      disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
    rm(genes.main.components.Hnf4a.Pecam1.Sox2)
    
  } # 5 Datasets MC + Hnf4a + Pecam1 + Sox2 
  
  {
    genes.main.components.Klf2.Pecam1.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Klf2','Nanog','Pecam1','Sox2')
    Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Klf2_Pecam1_Sox2', threshold_load_path = 'threshold_0.0_kanto_0.4',
                      saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Klf2_Pecam1_Sox2',
                      genes_modeled = genes.main.components.Klf2.Pecam1.Sox2, genes_displayed = genes.main.components.Klf2.Pecam1.Sox2, layout_matrix = layout_7genes_AGata6BNanog,
                      graph_title = '', threshold_vec = c(0.00001),
                      disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
    rm(genes.main.components.Klf2.Pecam1.Sox2)
    
  } # 5 Datasets MC + Klf2 + Pecam1 + Sox2 
  
  
  {
    genes.main.components.Hnf4a.Klf2.Pecam1.Sox2 <- c('Code','Fgf4','Fgfr2','Gata6','Hnf4a','Klf2','Nanog','Pecam1','Sox2')
    Qgraph.MergedData(load_path = '5_Datasets_maincomponents_with_Hnf4a_Klf2_Pecam1_Sox2', threshold_load_path = 'threshold_0.0_kanto_0.4',
                      saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/maincomponents_Hnf4a_Klf2_Pecam1_Sox2',
                      genes_modeled = genes.main.components.Hnf4a.Klf2.Pecam1.Sox2, genes_displayed = genes.main.components.Hnf4a.Klf2.Pecam1.Sox2, layout_matrix = layout_8genes_AGata6BNanog,
                      graph_title = '', threshold_vec = c(0.00001),
                      disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
    rm(genes.main.components.Hnf4a.Klf2.Pecam1.Sox2)
    
  } # 5 Datasets MC + Hnf4a + Klf2 + Pecam1 + Sox2 
  
  
  {
    Qgraph.MergedData(load_path = '5_Datasets_21_genes', threshold_load_path = 'threshold_0.0_kanto_0.4',
                         saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/5_Datasets_21_genes',
                         genes_modeled = genes.21, genes_displayed = genes.21, layout_matrix = layout_matrix_20_genes,
                         graph_title = '', threshold_vec = c(0.00001),
                         disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
  } # 5 Datasets All Genes
  
  {
    Qgraph.Only.One.Gene(load_path = '5_Datasets_21_genes', threshold_load_path = 'threshold_0.0_kanto_0.4',
                         saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/5_Datasets_21_genes_Fgf4',
                         genes_modeled = genes.21, genes_displayed = genes.21, gene_name = "Fgf4", layout_matrix = layout_matrix_20_genes,
                         graph_title = '', threshold_vec = c(0.00001),
                         disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
  } # 5 Datasets All Genes, Only showing Fgf4 interactions
  
  {
    Qgraph.Only.One.Gene(load_path = '5_Datasets_21_genes', threshold_load_path = 'threshold_0.0_kanto_0.4',
                         saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/5_Datasets_21_genes_Hnf4a',
                         genes_modeled = genes.21, genes_displayed = genes.21, gene_name = "Hnf4a", layout_matrix = layout_matrix_20_genes,
                         graph_title = '', threshold_vec = c(0.00001),
                         disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
  } # 5 Datasets All Genes, Only showing Hnf4a interactions
  
  {
    Qgraph.Only.One.Gene(load_path = '5_Datasets_21_genes', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/5_Datasets_21_genes_Nanog',
                        genes_modeled = genes.21, genes_displayed = genes.21, gene_name = "Nanog", layout_matrix = layout_matrix_20_genes,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
  } # 5 Datasets All Genes, Only showing Nanog interactions
  
  {
    Qgraph.Only.One.Gene(load_path = '5_Datasets_21_genes', threshold_load_path = 'threshold_0.0_kanto_0.4',
                         saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/5_Datasets_21_genes_Pecam1',
                         genes_modeled = genes.21, genes_displayed = genes.21, gene_name = "Pecam1", layout_matrix = layout_matrix_20_genes,
                         graph_title = '', threshold_vec = c(0.00001),
                         disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
  } # 5 Datasets All Genes, Only showing Pecam1 interactions
  
  {
    Qgraph.Only.One.Gene(load_path = '5_Datasets_21_genes', threshold_load_path = 'threshold_0.0_kanto_0.4',
                         saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/5_Datasets_21_genes_Sox2',
                         genes_modeled = genes.21, genes_displayed = genes.21, gene_name = "Sox2", layout_matrix = layout_matrix_20_genes,
                         graph_title = '', threshold_vec = c(0.00001),
                         disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
    
  } # 5 Datasets All Genes, Only showing Sox2 interactions
  
  
} # Qgraphs 
######### #########  ######### #########  ######### ######### ######### ######### 
######### #########  ######### #########  ######### ######### ######### ######### 


















load_path <- '5_Datasets_maincomponents_with_Hnf4a'     
threshold_load_path <- 'threshold_0.0_kanto_0.4'
saving_path <- 'threshold_0.0_kanto_0.4/5_Datasets/Article/Pres/maincomponents_Hnf4a'
path <- paste('DONE/',load_path,'/cardamom/',threshold_load_path,'/',sep = '')
print(path)

setwd("C:/Users/franc/CARDAMOM_OG/cardamom/")
genes_modeled <- c('Code','Fgf4','Fgfr2','Gata6','Hnf4a','Nanog')
dataset <- matrix.for.qgraphs_fun(path, 'inter.csv', genes_modeled)
qgraph1 <- matrix.for.qgraphs_fun(path, 'inter_1.csv', genes_modeled)
qgraph2 <- matrix.for.qgraphs_fun(path, 'inter_2.csv', genes_modeled)
qgraph3 <- matrix.for.qgraphs_fun(path, 'inter_3.csv', genes_modeled)
qgraph4 <- matrix.for.qgraphs_fun(path, 'inter_4.csv', genes_modeled)
genes_displayed <- genes_modeled
genes_displayed[1] <-  'Stimulus'
dataset <- dataset[genes_displayed,genes_displayed]
qgraph1 <- qgraph1[genes_displayed,genes_displayed]
qgraph2 <- qgraph2[genes_displayed,genes_displayed]
qgraph3 <- qgraph3[genes_displayed,genes_displayed]
qgraph4 <- qgraph4[genes_displayed,genes_displayed]


dataset['Stimulus',] <- 0
dataset['Fgf4','Fgf4'] <- 0
dataset['Fgf4','Fgfr2'] <- 0
dataset['Fgf4','Gata6'] <- 0
dataset['Fgfr2',] <- 0
dataset['Gata6','Fgf4'] <- 0
dataset['Hnf4a','Fgfr2'] <- 0
dataset

setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Gene Network Inference/")
print(paste("qgraphs/MergedData/",saving_path,"_GRN.pdf", sep = ''))
pdf(paste("qgraphs/MergedData/",saving_path,"_GRN.pdf", sep = ''))
qgraph(input = as.matrix(dataset), directed = TRUE, layout = layout_5genes_AGata6BNanog,
      title = paste('',' Starting at 8C - Thresh ', 0.000001 , sep = ''), threshold = 0.000001, cut = 0.000001*5)
dev.off()





{ # Merging multiple datasets with a new column indicating their source
  # Processing Allegre dataset
  Allegre.data <- cbind(
    Allegre.mRNA[!str_detect(Allegre.mRNA$Code, 'C16'), genes.20], # Filter rows excluding 'C16'
    Dataset = 'Allegre' # Add 'Dataset' column
  )
  
  # Processing Goolam dataset
  Goolam.data <- cbind(
    Goolam.mRNA[
      str_detect(Goolam.mRNA[, 'Code'], '8C') | 
        str_detect(Goolam.mRNA[, 'Code'], 'C16') | 
        str_detect(Goolam.mRNA[, 'Code'], 'C32'),
      genes.20
    ],
    Dataset = 'Goolam'
  )
  
  # Processing Guo dataset
  Guo.data <- cbind(
    Guo.mRNA.for.CARDAMOM[!(
      str_detect(Guo.mRNA.for.CARDAMOM[, 'Code'], '2C') | 
        str_detect(Guo.mRNA.for.CARDAMOM[, 'Code'], '4C')
    ), genes.20],
    Dataset = 'Guo'
  )
  
  # Processing Posfai dataset
  Posfai.data <- cbind(
    Posfai.mRNA[, genes.20], # Include all rows
    Dataset = 'Posfai'
  )
  
  # Processing Wang dataset
  Wang.data  <- cbind(
    Wang.mRNA[
      str_detect(Wang.mRNA[, 'Code'], '8C') | 
        str_detect(Wang.mRNA[, 'Code'], 'C16') | 
        str_detect(Wang.mRNA[, 'Code'], 'C32'),
      genes.20
    ],
    Dataset = 'Wang'
  )
  
  # Combine all processed datasets into one
  RtqpCR_Datasets_20_genes <- rbind(Goolam.data,Posfai.data,Wang.data)
  
  # Cleanup intermediate variables
  rm(Allegre.data,Goolam.data,Guo.data,Posfai.data,Wang.data)
  
} # Creation of Only scRNA-seq datasets

{
   RtqpCR_Datasets_20_genes
  map <- umap2(RtqpCR_Datasets_20_genes[,genes.20[2:length(genes.20)]])
  
  # `umap2` might return a matrix, so convert it to a data frame for ggplot compatibility.
  map_df <- as.data.frame(map)
  
  # Rename columns for easy plotting
  colnames(map_df) <- c("UMAP1", "UMAP2")
  
  # Add the desired color column from `Five_Datasets_21_genes`, replace "your_column_name" with the actual column name.
  map_df$Code <- RtqpCR_Datasets_20_genes$Code
  map_df$Dataset <- RtqpCR_Datasets_20_genes$Dataset
  
  # Plotting
  p1 <- ggplot(map_df, aes(x = UMAP1, y = UMAP2, color = Code)) +
    geom_point() +
    scale_color_manual(values = c('blue',"purple", "red","orange",'green')) + 
    theme_minimal() +
    theme(legend.position = "bottom", legend.title=element_text(size=20),legend.text=element_text(size=20),
          plot.title = element_text(hjust = 0.5, size = 35),
          axis.title = element_text(size = 25), axis.text = element_text(size = 25)) +
    guides(color = guide_legend(override.aes = list(size = 6))) + # Increase legend point size
    labs(title = "UMAP Plot Colored by Cellular Stage", color = "Cellular Stage")
  
  p2 <- ggplot(map_df, aes(x = UMAP1, y = UMAP2, color = Dataset)) +
    geom_point() +
    scale_color_manual(values = c('darkred','blue','darkgreen','orange','purple')) + 
    
    theme_minimal() +
    theme(legend.position = "bottom", legend.title=element_text(size=20),legend.text=element_text(size=20),
          plot.title = element_text(hjust = 0.5, size = 35),
          axis.title = element_text(size = 25), axis.text = element_text(size = 25)) +
    guides(color = guide_legend(override.aes = list(size = 6))) + # Increase legend point size
    labs(title = "UMAP Plot Colored by Dataset", color = "Dataset")
  
  p <- grid.arrange(p1,p2,nrow = 1)
  
  
  
  genes.main.components.with.Bmp4 <- c('Code','Bmp4','Fgf4','Fgfr2','Gata6','Nanog')
  deg_df_main_components_with_Bmp4 <- degradation_rates_df_21_genes[c(1,2,6,7,10,16),]
  deg_df_main_components_with_Bmp4
  
  automate_analysis(mRNA_data = RtqpCR_Datasets_20_genes[,genes.main.components.with.Bmp4],
                    levels = c('8C','C16','C32','C64','C90'), hours = c(0,12,18,24,30),
                    degradation_rates = deg_df_main_components_with_Bmp4,
                    saving_name = 'scRNA_seq_with_Bmp4',
                    bootstrap = FALSE, n_bootstraps = 0)
  rm(genes.main.components.with.Bmp4,deg_df_main_components_with_Bmp4)
  
  {
    
    Qgraph.scRNAseq <- function(load_path,threshold_load_path,saving_path,
                                   genes_modeled,genes_displayed,layout_matrix,
                                   graph_title,threshold_vec,disposition_full = c(2,2),
                                   nbr_of_stage_transitions_or_columns = 4){
      
      path <- paste('DONE/',load_path,'/cardamom/',threshold_load_path,'/',sep = '')
      print(path)
      
      setwd("C:/Users/franc/CARDAMOM_OG/cardamom/")
      dataset <- matrix.for.qgraphs_fun(path, 'inter.csv', genes_modeled)
      qgraph1 <- matrix.for.qgraphs_fun(path, 'inter_1.csv', genes_modeled)
      qgraph2 <- matrix.for.qgraphs_fun(path, 'inter_2.csv', genes_modeled)
      qgraph3 <- matrix.for.qgraphs_fun(path, 'inter_3.csv', genes_modeled)
      genes_displayed[1] <-  'Stimulus'
      dataset <- dataset[genes_displayed,genes_displayed]
      qgraph1 <- qgraph1[genes_displayed,genes_displayed]
      qgraph2 <- qgraph2[genes_displayed,genes_displayed]
      qgraph3 <- qgraph3[genes_displayed,genes_displayed]
      setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Gene Network Inference/")
      print(paste("qgraphs/MergedData/",saving_path,"_GRN.pdf", sep = ''))
      pdf(paste("qgraphs/MergedData/",saving_path,"_GRN.pdf", sep = ''))
      par(mfrow = disposition_full) 
      for (i in 1:length(threshold_vec)){
        qgraph(input = as.matrix(dataset), directed = TRUE, layout = layout_matrix,
               #title = paste(graph_title,' Starting at 8C - Thresh ', threshold_vec[i] , sep = ''),
               threshold = threshold_vec[i], cut = threshold_vec[i]*5)
      }
      dev.off()
      
      
      print(paste("qgraphs/MergedData/",saving_path,"_GRN_By_Stage.pdf",sep = ''))
      pdf(paste("qgraphs/MergedData/",saving_path,"_GRN_By_Stage.pdf",sep = ''))
      par(mfrow = c(length(threshold_vec), nbr_of_stage_transitions_or_columns))
      for (i in 1:length(threshold_vec)){
        qgraph(input = as.matrix(qgraph1), directed = TRUE, layout = layout_matrix,
               title = paste('8C -> 16C'),# - Thresh ', threshold_vec[i] , sep = ''), 
               threshold = threshold_vec[i], cut = threshold_vec[i]*5)
        qgraph(input = as.matrix(qgraph2), directed = TRUE, layout = layout_matrix,
               title = paste('16C -> 32C'),# - Thresh ', threshold_vec[i] , sep = ''), 
               threshold = threshold_vec[i], cut = threshold_vec[i]*5)
        qgraph(input = as.matrix(qgraph3), directed = TRUE, layout = layout_matrix,
               title = paste('32C -> 64C'),# - Thresh ', threshold_vec[i] , sep = ''), 
               threshold = threshold_vec[i], cut = threshold_vec[i]*5)
      }
      dev.off()
    }
    
    
    {
      genes.main.components.Bmp4 <- c('Code','Bmp4','Fgf4','Fgfr2','Gata6','Nanog')
      Qgraph.scRNAseq(load_path = 'scRNA_seq_with_Bmp4', threshold_load_path = 'threshold_0.0_kanto_0.4',
                        saving_path = 'threshold_0.0_kanto_0.4/5_Datasets/Article/OnlyscRNAseq/maincomponents_Bmp4',
                        genes_modeled = genes.main.components.Bmp4, genes_displayed = genes.main.components.Bmp4, layout_matrix = layout_5genes_BFgf4,
                        graph_title = '', threshold_vec = c(0.00001),
                        disposition_full = c(1,1), nbr_of_stage_transitions_or_columns = 1)
      
      rm(genes.main.components.Bmp4)
      
    } # 5 Datasets MC + Bmp4
    
  } # 5 Datasets MC
  
  
} # Stuff with only scRNA seq datasets, not well done




