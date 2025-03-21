######################## Entropy #########################
## In this file, we generate the inter-cellular entropy for the 21 genes considered and do the heatmaps, barplots.

######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### Setting the directory for loading data ######### ######### 
setwd("C:/Users/franc/Dropbox/Université/Embryo/Francisco/Entropy")
load('FinalData/InterCellularEntropy/EntropyGamma/PreProcessedData.RData')
######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### 




######### ######### ######### ######### ######### ######### ######### #########  
######### ######### Load necessary libraries ######### ######### 
library(dplyr)
library(stringr)
library(tidyr)
library(fitdistrplus)
library(data.table)
library(tibble)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(viridis)
######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### 

######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### Genes chosen and degradation rates ######### ######### 
genes.20 <- c('Code','Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6',
              'Hnf4a','Id2','Klf2', 'Klf4','Klf5','Nanog','Pdgfra','Pecam1','Pou5f1','Sox17','Sox2')
genes.21 <- c('Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6',
                                'Hnf4a','Id2','Klf2', 'Klf4','Klf5','Nanog','Pdgfa','Pdgfra','Pecam1','Pou5f1','Sox17','Sox2')
genes.21.with.code.and.pop <- c('Code','Population','Bmp4','Cdx2','Dab2','Esrrb','Fgf4','Fgfr2','Gata3','Gata4','Gata6',
                                'Hnf4a','Id2','Klf2', 'Klf4','Klf5','Nanog','Pdgfa','Pdgfra','Pecam1','Pou5f1','Sox17','Sox2')
######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### 




######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### Apply Scaling Factors to datasets and keep only stages we want ######### ######### 
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
  Five_Datasets_20_genes <- rbind(Allegre.data,Goolam.data,Guo.data,Posfai.data,Wang.data)
  
  # Cleanup intermediate variables
  rm(Allegre.data,Goolam.data,Guo.data,Posfai.data,Wang.data)
  
  All.Data <- Five_Datasets_20_genes[,genes.20]
  #All.Data[,2:length(genes.20)] <- All.Data[,2:length(genes.20)]/10 #Scaling Factor necessary
  All.Data <- All.Data[order(All.Data$Code), ]
  
} # Creation of Five_Datasets_20_genes and scaling

{
  Allegre.data <- Allegre.mRNA[, genes.21.with.code.and.pop]
  # Allegre.data[,3:ncol(Allegre.data)] <- Allegre.data[,3:ncol(Allegre.data)]/10
  head(Allegre.data)
  
  Guo.data <- Guo.mRNA.for.entropy[str_detect(Guo.mRNA.for.entropy$Code,'8C') | str_detect(Guo.mRNA.for.entropy$Code,'C16') | str_detect(Guo.mRNA.for.entropy$Code,'C32') | str_detect(Guo.mRNA.for.entropy$Code,'C64'), genes.21.with.code.and.pop]
  head(Guo.data)
  
  Goolam.data <- Goolam.mRNA[str_detect(Goolam.mRNA[,'Code'],'8C') | str_detect(Goolam.mRNA$Code,'C16') | str_detect(Goolam.mRNA$Code,'C32'), genes.21.with.code.and.pop]
  # Goolam.data[,3:ncol(Goolam.data)] <- Goolam.data[,3:ncol(Goolam.data)]/10
  head(Goolam.data)
  
  Posfai.data <- Posfai.mRNA[, genes.21.with.code.and.pop]
  # Posfai.data[,3:ncol(Posfai.data)] <- Posfai.data[,3:ncol(Posfai.data)]/10
  head(Posfai.data)
  
  Wang.data <- Wang.mRNA[str_detect(Wang.mRNA[,'Code'],'8C') | str_detect(Wang.mRNA$Code,'C16') | str_detect(Wang.mRNA$Code,'C32'), genes.21.with.code.and.pop]
  head(Wang.data)
} # Creation of Individual datasets and scaling
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 




######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### ######### ######### Define the functions ######### ######### ######### #########
gamma_entropy <- function(shape, rate) {
  return(shape + log(gamma(shape)) + (1 - shape) * digamma(shape) - log(rate))
} # function to compute the inter-cellular gamma entropy
compute_entropy_with_bootstrap <- function(data, n_bootstrap = 100) {
  # Identify gene columns (assume non-gene columns are "Code" and "Population")
  gene_columns <- setdiff(names(data), c("Code", "Population"))
  
  # Get unique stages from the "Code" column
  stages <- unique(data$Code)
  
  # Initialize lists to store bootstrap results for entropy, KS statistic, and KS p-value
  entropy_results <- list()
  entropy_sd_results <- list()
  ks_stat_results <- list()
  ks_pval_results <- list()
  
  # Loop through each stage
  for (stage in stages) {
    # Subset data for the current stage
    stage_data <- data %>% filter(Code == stage)
    
    # Initialize lists to store bootstrap results for this stage
    stage_entropy <- list()
    stage_ks_stats <- list()
    stage_ks_pvals <- list()
    
    # Perform bootstrap iterations
    for (b in 1:n_bootstrap) {
      # Resample data with replacement
      resampled_data <- stage_data[sample(1:nrow(stage_data), replace = TRUE), ]
      
      # Compute entropy, KS statistic, and p-value for all genes in the resampled dataset
      resampled_results <- sapply(gene_columns, function(gene) {
        gene_values <- resampled_data[[gene]]
        entropy <- NA  # Default to NA for entropy
        ks_stat <- NA  # Default to NA for KS statistic
        ks_pval <- NA  # Default to NA for KS p-value
        
        # Proceed if there are enough non-NA values
        if (length(gene_values[!is.na(gene_values)]) > 1) {
          # Try to fit a Gamma distribution and calculate entropy and KS test
          tryCatch({
            # Fit Gamma distribution and compute entropy
            fit <- fitdist(gene_values, "gamma", method = "mme")
            entropy <- gamma_entropy(shape = fit$estimate["shape"], rate = fit$estimate["rate"])
            if (is.infinite(entropy)) {
              entropy <- NA
            }
            
            cat("\r",b)
            # Perform KS test comparing the gene values with the fitted Gamma distribution
            ks_result <- ks.test(gene_values , "pgamma", shape = fit$estimate["shape"], rate = fit$estimate["rate"])
            ks_stat <- ks_result$statistic
            ks_pval <- ks_result$p.value
            
            #ks_stat <- fit$aic
            #ks_pval <- fit$bic
          }, error = function(e) {
            message(paste("Failed to compute entropy and KS for gene", gene, "at stage", stage, ":", e$message))
          })
        }
        
        return(c(entropy, ks_stat, ks_pval))  # Return all three values
      })
      
      # Store the results for this bootstrap iteration (transpose so each row is a gene)
      stage_entropy[[b]] <- resampled_results[1,]
      stage_ks_stats[[b]] <- resampled_results[2,]
      stage_ks_pvals[[b]] <- resampled_results[3,]
    }
    
    # Combine bootstrap results into data frames
    stage_entropy_df <- do.call(rbind, stage_entropy) %>% as.data.frame() %>% setNames(paste0(gene_columns))
    stage_ks_stats_df <- do.call(rbind, stage_ks_stats) %>% as.data.frame() %>% setNames(paste0(gene_columns))
    stage_ks_pvals_df <- do.call(rbind, stage_ks_pvals) %>% as.data.frame() %>% setNames(paste0(gene_columns))
    
    # Compute summary statistics (mean, standard deviation) for each gene
    stage_entropy_mean <- as.data.frame(lapply(stage_entropy_df, function(x) mean(x, na.rm = TRUE)))
    stage_ks_stat_mean <- as.data.frame(lapply(stage_ks_stats_df, function(x) mean(x, na.rm = TRUE)))
    stage_ks_pval_mean <- as.data.frame(lapply(stage_ks_pvals_df, function(x) mean(x, na.rm = TRUE)))
    
    stage_entropy_sd <- as.data.frame(lapply(stage_entropy_df, function(x) sd(x, na.rm = TRUE)))
    #stage_ks_stat_sd <- as.data.frame(lapply(stage_ks_stats_df, function(x) sd(x, na.rm = TRUE)))
    #stage_ks_pval_sd <- as.data.frame(lapply(stage_ks_pvals_df, function(x) sd(x, na.rm = TRUE)))
    
    # Add stage information to the summaries
    stage_entropy_mean <- cbind(Stage = stage, stage_entropy_mean)
    stage_ks_stat_mean <- cbind(Stage = stage, stage_ks_stat_mean)
    stage_ks_pval_mean <- cbind(Stage = stage, stage_ks_pval_mean)
    
    stage_entropy_sd <- cbind(Stage = stage, stage_entropy_sd)
    #stage_ks_stat_sd <- cbind(Stage = stage, stage_ks_stat_sd)
    #stage_ks_pval_sd <- cbind(Stage = stage, stage_ks_pval_sd)
    
    # Store results in respective lists
    entropy_results[[stage]] <- stage_entropy_mean
    entropy_sd_results[[stage]] <- stage_entropy_sd
    ks_stat_results[[stage]] <- stage_ks_stat_mean
    ks_pval_results[[stage]] <- stage_ks_pval_mean
    
    
  }
  
  # Combine results for all stages into separate data frames
  entropy_df <- do.call(rbind, entropy_results)
  entropy_sd_df <- do.call(rbind, entropy_sd_results)
  ks_stat_df <- do.call(rbind, ks_stat_results)
  ks_pval_df <- do.call(rbind, ks_pval_results)
  
  return(list(entropy = entropy_df, entropy_sd = entropy_sd_df, ks_stat = ks_stat_df, ks_pval = ks_pval_df))
} # function to compute the entropy with bootstrap
prepare_entropy_by_dataset_for_plot <- function(all_data_list, all_data_list_sd, all_data_list_ks, all_data_list_ks_pval) {
  
  # Transform and merge all datasets
  plot_data <- bind_rows(lapply(names(all_data_list), function(dataset) {
    
    df <- all_data_list[[dataset]]
    df_sd <- all_data_list_sd[[dataset]]
    df_ks <- all_data_list_ks[[dataset]]
    df_ks_pval<- all_data_list_ks_pval[[dataset]]
    
    # Transform entropy data
    df <- df %>%
      mutate(Dataset = dataset) %>%
      pivot_longer(cols = -c(Stage, Dataset), names_to = "Gene", values_to = "Entropy")
    
    # Transform standard deviation data
    df_sd <- df_sd %>%
      mutate(Dataset = dataset) %>%
      pivot_longer(cols = -c(Stage, Dataset), names_to = "Gene", values_to = "Entropy_SD")
    
    # Transform ks data
    df_ks <- df_ks %>%
      mutate(Dataset = dataset) %>%
      pivot_longer(cols = -c(Stage, Dataset), names_to = "Gene", values_to = "KS")
    
    # Transform ks data
    df_ks_pval <- df_ks_pval %>%
      mutate(Dataset = dataset) %>%
      pivot_longer(cols = -c(Stage, Dataset), names_to = "Gene", values_to = "KS_pval")
    
    # Merge entropy and SD data
    a <- left_join(df, df_sd, by = c("Stage", "Dataset", "Gene"))
    
    # Merge ks and ks pval data
    b <- left_join(df_ks, df_ks_pval, by = c("Stage", "Dataset", "Gene"))
    
    # Merge all data
    left_join(a, b, by = c("Stage", "Dataset", "Gene"))
    
  }))
  
  return(plot_data)
} # function to prepare the entropy by dataset 
prepare_heatmap_data <- function(data, value_col, na_fill = NA) {
  data %>%
    # Combine Stage and Dataset for column names
    mutate(Stage_Dataset = paste(Stage, Dataset, sep = "_")) %>%
    # Pivot to wide format
    pivot_wider(
      id_cols = Gene, 
      names_from = Stage_Dataset, 
      values_from = all_of(value_col)
    ) %>%
    # Replace missing values with NA or a specified value
    replace(is.na(.), na_fill) %>%
    # Convert to a matrix
    column_to_rownames("Gene") %>%
    as.matrix()
} # Prepare data for heatmaps
generate_heatmap <- function(data_matrix, title, na_color = "grey") {
  pheatmap(
    data_matrix,
    color = inferno(100), # Customize color palette
    cluster_rows = FALSE,
    cluster_cols = FALSE,  # Do not cluster columns, retain stage/dataset order
    na_col = na_color,     # Set color for NA values
    main = title
  )
} # Generate heatmaps
entropy_tendencies <- function(data, code_labels){
  
  # Select only the entropy columns (exclude 'Code' column)
  entropy_data <- data$entropy[, -1]
  

  # Compute row-wise differences
  row_differences <- as.data.frame(diff(as.matrix(entropy_data)))
  row_differences$Code <- code_labels
  
  row_differences <- row_differences[, c(ncol(row_differences), 1:(ncol(row_differences) - 1))]
  rownames(row_differences) <- c()
  return(row_differences)
} # function to compute the entropy tendencies
barplots_of_comparisons <- function(list_data,gene_of_interest){
  
  # List to store genes with the same tendency as Nanog in each dataset
  matching_genes_list <- list()
  
  # List of dataset names
  list_data_name <- list("Allègre","Guo","Goolam","Posfai","Wang")
  
  # Loop through all datasets
  for (i in seq_along(list_data)) {
    
    df <- list_data[[i]]  # Get the dataset
    
    # Find available stages (we focus only on C16, C32, C64)
    available_stages <- intersect(df$Stage, c("C16", "C32", "C64"))
    
    # If there are fewer than 2 stages, we can't compute tendencies
    if (length(available_stages) < 2) {
      matching_genes_list[list_data_name[[i]]] <- character(0)
      next
    }
    
    # Filter dataset to keep only available stages
    df_filtered <- df[df$Stage %in% available_stages, ]
    
    # Remove the Stage column
    df_filtered <- df_filtered[, -1]
    
    # Compute tendencies (differences between consecutive stages)
    tendencies <- as.data.frame(diff(as.matrix(df_filtered)))
    
    # Check if Nanog is present in the dataset
    if (!"Nanog" %in% colnames(tendencies)) {
      matching_genes_list[list_data_name[[i]]] <- character(0)
      next
    }
    
    # Extract Nanog's pattern for the available stages
    nanog_trend <- sign(tendencies$Nanog)
    
    # Identify genes that have the same tendency pattern as Nanog
    matching_genes <- colnames(tendencies)[apply(tendencies, 2, function(gene_trend) {
      # Compare only for the stages that exist
      all(sign(gene_trend) == nanog_trend, na.rm = TRUE)
    })]
    # Store results for the dataset
    matching_genes_list[[list_data_name[[i]]]] <- matching_genes
    
  }
  
  # Convert matching_genes_list into a data frame
  matching_genes_df <- do.call(rbind, lapply(names(matching_genes_list), function(dataset) {
    data.frame(Dataset = dataset, Gene = matching_genes_list[[dataset]], stringsAsFactors = FALSE)
  }))
  
  # Remove empty datasets (where no genes matched Nanog)
  matching_genes_df <- matching_genes_df[matching_genes_df$Gene != "", ]
  
  # Count occurrences of each gene across datasets
  gene_counts <- matching_genes_df %>%
    group_by(Gene, Dataset) %>%
    summarise(Count = n(), .groups = "drop")
  
  # Sort genes alphabetically (instead of by total occurrences)
  gene_counts$Gene <- factor(gene_counts$Gene, levels = rev(sort(unique(gene_counts$Gene))))
  
  # **Exclude Nanog** from the plot (optional)
  gene_counts <- gene_counts %>% filter(Gene != gene_of_interest)
  
  
  return(gene_counts)
} # This creates the data to be able to do the stacked barplots comparing tendencies between datasets.
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 




######### ######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### ######### Compute the entropy and KS ######### ######### ######### #########
Allegre.entropy <- compute_entropy_with_bootstrap(Allegre.data, n_bootstrap = 1000) 
Guo.entropy <- compute_entropy_with_bootstrap(Guo.data, n_bootstrap = 1000)
Guo.entropy$entropy[,2:ncol(Guo.entropy$entropy)] <- Guo.entropy$entropy[,2:ncol(Guo.entropy$entropy)] + 7
Goolam.entropy <- compute_entropy_with_bootstrap(Goolam.data, n_bootstrap = 1000)
Posfai.entropy <- compute_entropy_with_bootstrap(Posfai.data, n_bootstrap = 1000)
Wang.entropy <- compute_entropy_with_bootstrap(Wang.data, n_bootstrap = 1000)
Entropy.all.datasets <- compute_entropy_with_bootstrap(All.Data, n_bootstrap = 1000)
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
load('FinalData/InterCellularEntropy/EntropyGamma/DataEntropyGamma.RData')


######### ######### ######### ######### ######### ######### ######### ######### #########
######### ######### ######### ######### Plot Data ######### ######### ######### #########
{
  all_data <- list(
  Allegre = Allegre.entropy$entropy,
  Guo = Guo.entropy$entropy,
  Goolam = Goolam.entropy$entropy,
  Posfai = Posfai.entropy$entropy,
  Wang = Wang.entropy$entropy
  )

  all_data_sd <- list(
    Allegre = Allegre.entropy$entropy_sd,
    Guo = Guo.entropy$entropy_sd,
    Goolam = Goolam.entropy$entropy_sd,
    Posfai = Posfai.entropy$entropy_sd,
    Wang = Wang.entropy$entropy_sd
  )
  
  all_data_ks <- list(
    Allegre = Allegre.entropy$ks_stat,
    Guo = Guo.entropy$ks_stat,
    Goolam = Goolam.entropy$ks_stat,
    Posfai = Posfai.entropy$ks_stat,
    Wang = Wang.entropy$ks_stat
  )
  
  all_data_ks_pval <- list(
    Allegre = Allegre.entropy$ks_pval,
    Guo = Guo.entropy$ks_pval,
    Goolam = Goolam.entropy$ks_pval,
    Posfai = Posfai.entropy$ks_pval,
    Wang = Wang.entropy$ks_pval
  )
} # Combine the datasets into a list

{
  plot_data_by_dataset <- prepare_entropy_by_dataset_for_plot(all_data, all_data_sd, all_data_ks, all_data_ks_pval) 
  plot_data_all_together <- prepare_entropy_by_dataset_for_plot(list(Full = Entropy.all.datasets$entropy), list(Full = Entropy.all.datasets$entropy_sd),
                                                                list(Full = Entropy.all.datasets$ks_stat), list(Full = Entropy.all.datasets$ks_pval)) 
  # Define the colors for the datasets
  dataset_colors <- c('darkred', 'blue', 'darkgreen', 'orange', 'purple')
  
} # Create the datasets for ggplots and heatmaps

{
  # Prepare matrices for KS values and P-values
  ks_matrix <- prepare_heatmap_data(plot_data_by_dataset, "KS", na_fill = NA)
  ks_matrix <- ks_matrix[, order(colnames(ks_matrix))]
  # Replace "C" at the start of column names with an empty string and add "C" at the end
  colnames(ks_matrix) <- gsub("^C(\\d+)", "\\1C", colnames(ks_matrix))
  
  pval_matrix <- prepare_heatmap_data(plot_data_by_dataset, "KS_pval", na_fill = NA)
  pval_matrix <- pval_matrix[,order(colnames(pval_matrix))]
  # Replace "C" at the start of column names with an empty string and add "C" at the end
  colnames(pval_matrix) <- gsub("^C(\\d+)", "\\1C", colnames(pval_matrix))
  
  
  # Generate heatmaps
  ks_pheatmap <- generate_heatmap(ks_matrix, title = "KS Statistic Heatmap")
  ks_pval_pheatmap <- generate_heatmap(pval_matrix, title = "P-value Heatmap")
  
  p <- grid.arrange(ks_pheatmap$gtable,ks_pval_pheatmap$gtable, nrow = 2)
  
  ggsave('GGplotsFigs/EntropyGamma/Heatmaps/Correlations/KS_Stat_boot1000.pdf', p, width = 5, height = 10)
  
} # Heatmap of KS and PValues By Dataset

{
  # Prepare matrices for KS values and P-values
  ks_matrix <- prepare_heatmap_data(plot_data_all_together, "KS", na_fill = NA)
  ks_matrix <- ks_matrix[, order(colnames(ks_matrix))]
  # Replace "C" at the start of column names with an empty string and add "C" at the end
  colnames(ks_matrix) <- gsub("^C(\\d+)", "\\1C", colnames(ks_matrix))
  
  pval_matrix <- prepare_heatmap_data(plot_data_all_together, "KS_pval", na_fill = NA)
  pval_matrix <- pval_matrix[,order(colnames(pval_matrix))]
  # Replace "C" at the start of column names with an empty string and add "C" at the end
  colnames(pval_matrix) <- gsub("^C(\\d+)", "\\1C", colnames(pval_matrix))
  
  
  # Generate heatmaps
  ks_pheatmap <- generate_heatmap(ks_matrix, title = "KS Statistic Heatmap")
  ks_pval_pheatmap <- generate_heatmap(pval_matrix, title = "P-value Heatmap")
  
  p <- grid.arrange(ks_pheatmap$gtable,ks_pval_pheatmap$gtable, nrow = 1)
  
  ggsave('GGplotsFigs/EntropyGamma/Heatmaps/Correlations/KS_Stat_boot1000_All_Data.pdf', p, width = 5, height = 5)
  
} # Heatmap of KS and Pvalues for Full Dataset

{
  
  pdf("GGplotsFigs/EntropyGamma/HistogramFittings/Esrrb/Posfai_C64.pdf", width = 8, height = 8)
  plot(fitdist(Posfai.data[Posfai.data$Code == 'C64','Esrrb'], distr = 'gamma', method = 'mme'))
  dev.off()
} # Plots of Gamma Distributions generated

{
  data.for.plot <- plot_data_by_dataset[str_detect(plot_data_by_dataset$Gene, 'Nanog') |
                                          str_detect(plot_data_by_dataset$Gene, 'Gata6') |
                                          str_detect(plot_data_by_dataset$Gene, 'Fgf4') |
                                          str_detect(plot_data_by_dataset$Gene, 'Fgfr2'),]
  
  # Ensure the Gene column is a factor with the specified order
  data.for.plot$Gene <- factor(data.for.plot$Gene, levels = c("Nanog", "Gata6", "Fgf4", "Fgfr2"))
  
  
  p1 <- ggplot(data = data.for.plot,
              aes(x = factor(Stage, levels = c('8C', 'C16', 'C32', 'C64', 'C90')),
                  y = Entropy)) +
    geom_line(aes(group = Dataset, color = Dataset), linewidth = 1.5) +
    scale_color_manual(values = dataset_colors, name = 'Dataset Name', drop = FALSE) +
    geom_errorbar(aes(ymin = Entropy - Entropy_SD,
                      ymax = Entropy + Entropy_SD), linewidth = 0.6, width = 0.5) +
    facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none') +
    scale_x_discrete(labels = c('8C', '16C', '32C', '64C', '90C'), drop = FALSE) +
    scale_fill_discrete(drop = FALSE) +
    xlab('') +
    ylab("Inter-cellular entropy")
  p1
  
  ggsave('GGplotsFigs/EntropyGamma/GenePlots/MC.pdf', p1, width = 20, height = 5)
  
  p2 <- ggplot(data = plot_data_all_together,
              aes(x = factor(Stage, levels = c('8C', 'C16', 'C32', 'C64', 'C90')),
                  y = Entropy)) +
    geom_line(aes(group = Dataset), linewidth = 1.5) +
    scale_color_manual(values = dataset_colors, name = 'Dataset Name', drop = FALSE) +
    geom_errorbar(aes(ymin = Entropy - Entropy_SD,
                      ymax = Entropy + Entropy_SD), linewidth = 1, width = 0.5) +
    facet_wrap(~ Gene, scales = "free_y", nrow = 8) +
    theme(text = element_text(size = 25),
          axis.text.x = element_text(size = 25),
          plot.title = element_text(hjust = 0.5),
          legend.position = 'none') + 
    ggtitle("Entropy Variation by Gene") +
    scale_x_discrete(labels = c('8C', '16C', '32C', '64C', '90C'), drop = FALSE) +
    scale_fill_discrete(drop = FALSE) +
    xlab('') +
    ylab("Inter-cellular Entropy")
  p2
  #grid.arrange(p1,p2)
  
} # Entropy plots by Dataset and Full Dataset

{
  # Define row names for different datasets
  rownames_list <- list(
    c("16C -> 32C", "32C -> 64C", "64C -> 90C"),
    c("8C -> 16C", "16C -> 32C", "32C -> 64C"),
    c("8C -> 16C", "16C -> 32C"),
    c("16C -> 32C", "32C -> 64C"),
    c("8C -> 16C", "16C -> 32C")
  )
  p_list <- list()
  
  list_data_name <- list("Allègre","Guo","Goolam","Posfai","Wang")
  # Define custom breaks from -3 to 3
  breaksList <- seq(-3, 3, length.out = 50)
  
  # Loop through all datasets
  for (i in seq_along(all_data)) {
    
    df <- all_data[[i]]  # Get the dataset
    
    # Remove the Stage column
    stage_names <- df$Stage
    df <- df[, -1]
    
    # Compute the differences between consecutive stages
    tendencies <- as.data.frame(diff(as.matrix(df)))
    
    # Rename rows to represent stage transitions
    rownames(tendencies) <- rownames_list[[i]]  # Correct indexing
    
    # Create a heatmap with capped color scale
    p_list[[i]] <- pheatmap(
        t(tendencies),  # Transpose for correct orientation
        cluster_rows = FALSE, cluster_cols = FALSE,
        scale = "none",  
        color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(50), 
        breaks = breaksList,  # Apply correct breaks
        main = paste(list_data_name[[i]]),
        angle_col = 90, 
        legend = FALSE
      )
    }
  
  p <- grid.arrange(p_list[[1]]$gtable,p_list[[2]]$gtable,p_list[[3]]$gtable,p_list[[4]]$gtable,p_list[[5]]$gtable,
                    nrow = 1)
  
  ggsave("GGplotsFigs/EntropyGamma/Heatmaps/Tendencies/Entropy/All_Data.pdf",p , width = 13, height = 7)

} # Heatmaps of tendencies

{
  data.for.barplot <- barplots_of_comparisons(all_data,"Nanog")
  
  # Define custom colors for datasets
  dataset_colors <- c("Allègre" = "darkred", 
                      "Guo" = "darkgreen", 
                      "Goolam" = "darkblue", 
                      "Posfai" = "orange", 
                      "Wang" = "purple")  
  
  # Create stacked bar plot
  p <- ggplot(data.for.barplot, aes(x = Gene, y = Count, fill = Dataset)) +
        geom_bar(stat = "identity", position = "stack", width = 0.7, color = "black") +  # **Black outline for boxes**
        scale_fill_manual(values = dataset_colors) +
        labs(x = "Gene",
             y = "Number of Matching Datasets",
             fill = "Dataset") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
              axis.text.y = element_text(size = 12),
              legend.position = "right",
              legend.text = element_text(size = 12)) +
        coord_flip()  # Flip for better readability if needed
  
  ggsave("GGplotsFigs/EntropyGamma/Barplots/GeneProfileComparisons/Nanog/Entropy_boot1000.pdf", p, width = 4, height = 7)
  
} # Barplots of comparison of tendencies.
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 
######### ######### ######### ######### ######### ######### ######### ######### ######### ######### 

