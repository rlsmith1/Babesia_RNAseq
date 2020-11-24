


# libraries ---------------------------------------------------------------

      library(DESeq2)
      library(tidyverse)
      library(magrittr)
      library(pheatmap)
      library(RColorBrewer)
      library(plotly)
      library(stringi)
      library(purrr)
      library(ggrepel)



# data --------------------------------------------------------------------

      load("objects/analysis_data.Rdata")



# Create DESeq2 object ----------------------------------------------------


      # round counts to nearest integer
      df_rawcounts_l[3:ncol(df_rawcounts_l)] <- round(df_rawcounts_filtered[3:ncol(df_rawcounts_l)], digits = 0)
      df_rawcounts_h[3:ncol(df_rawcounts_h)] <- round(df_rawcounts_filtered[3:ncol(df_rawcounts_h)], digits = 0)
      
      #convert to matrix
      m_rawcounts_l <- df_rawcounts_l %>% dplyr::select(-c("geneID", "ensgene")) %>% as.matrix()
      rownames(m_rawcounts_l) <- df_rawcounts_l$ensgene
      
      m_rawcounts_h <- df_rawcounts_h %>% dplyr::select(-c("geneID", "ensgene")) %>% as.matrix()
      rownames(m_rawcounts_h) <- df_rawcounts_h$ensgene
      
      #create object - ***design indicates condition of interest for differential testing!!!!****
      dds_dogs_l <- DESeqDataSetFromMatrix(countData = m_rawcounts_l, colData = df_metadata %>% filter(dose == "low"), design = ~ day)
      dds_dogs_h <- DESeqDataSetFromMatrix(countData = m_rawcounts_h, colData = df_metadata  %>% filter(dose == "high"), design = ~ day)
      
      

# Normalize counts (median of ratios) --------------------------------------------------------

  
      #estimate size factors
      dds_dogs_l <- estimateSizeFactors(dds_dogs_l)
      dds_dogs_h <- estimateSizeFactors(dds_dogs_h)
      
      #extraction
      m_norm_counts_l <- counts(dds_dogs_l, normalized = TRUE)
      m_norm_counts_h <- counts(dds_dogs_h, normalized = TRUE)
      
      

# QC check ----------------------------------------------------------------

      
      # 1. Hierarchical heatmap

            # improve visualization
            vsd_l <- vst(dds_dogs_l, blind = TRUE)
            vsd_h <- vst(dds_dogs_h, blind = TRUE)
      
            # extract vst matrix from the object
            vsd_mat_l <- assay(vsd_l) 
            vsd_mat_h <- assay(vsd_h) 
            
            # relabel columns
            colnames(vsd_mat_l) <- df_metadata %>% filter(dose == "low") %>% .$sample
            colnames(vsd_mat_h) <- df_metadata %>% filter(dose == "high") %>% .$sample
            
            # compute pairwise correlation values
            vsd_cor_l <- cor(vsd_mat_l)
            vsd_cor_h <- cor(vsd_mat_h)
            
            # plot heatmap
            pheatmap(vsd_cor_l)
            pheatmap(vsd_cor_h)
            
      # 2. PCA
      
            pcaData_l <- plotPCA(vsd_l, intgroup = c("day", "dog"), returnData = TRUE)
            pcaData_h <- plotPCA(vsd_h, intgroup = c("day", "dog"), returnData = TRUE)
            
            ggplot(pcaData_l, aes(PC1, PC2, color = day)) +
                  geom_text_repel(aes(label = dog), color = "black", size = 5) +
                  geom_point(size = 3)
            
            ggplot(pcaData_h, aes(PC1, PC2, color = day), text = "dog") +
                  geom_text_repel(aes(label = dog), color = "black", size = 5) +
                  geom_point(size = 3)
            

            
# Differential expression analysis ----------------------------------------
            
            
      de_l <- DESeq(dds_dogs_l)  
      de_h <- DESeq(dds_dogs_h)  
      

      

# How well do the data fit the model --------------------------------------

      
      # Mean-variance relationship
      
            # low
            mean_counts_l <-  df_rawcounts_l[3:ncol(df_rawcounts_l)] %>% rowMeans() 
            var_counts_l <- df_rawcounts_l[3:ncol(df_rawcounts_l)] %>% as.matrix() %>% rowVars() 
            
            df_mean_var_l <- data.frame(mean_counts_l, var_counts_l)
            
            ggplot(df_mean_var_l) +
                  
                  geom_point(aes(mean_counts_l, var_counts_l)) +
                  scale_y_log10() +
                  scale_x_log10() +
                  xlab("Mean counts per gene") +
                  ylab("Variance per gene")
            
            # high
            mean_counts_h <-  df_rawcounts_h[3:ncol(df_rawcounts_h)] %>% rowMeans() 
            var_counts_h <- df_rawcounts_h[3:ncol(df_rawcounts_h)] %>% as.matrix() %>% rowVars() 
            
            df_mean_var_h <- data.frame(mean_counts_h, var_counts_h)
            
            ggplot(df_mean_var_h) +
                  
                  geom_point(aes(mean_counts_h, var_counts_h)) +
                  scale_y_log10() +
                  scale_x_log10() +
                  xlab("Mean counts per gene") +
                  ylab("Variance per gene")

      # Dispersion  
      plotDispEsts(de_l) 
      plotDispEsts(de_h) 
      

      
# Differential expression analysis ----------------------------------------


      # Set contrasts and return results
      df_degs_l <- results(de_l, 
                           alpha = 0.05, # FDR threshold
                           contrast = c("day", 4, 0), # compare day 4 to day 0
                           lfcThreshold = log2(1.5), # 1.5 absolute FC threshold
                           independentFiltering = FALSE) %>% 
         data.frame() %>% 
         rownames_to_column("ensgene") %>% 
         as_tibble() %>% 
         na.omit() %>% 
         left_join(df_cpm_l, by = "ensgene") %>% 
         
         filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))

      df_degs_h <- results(de_h, 
                           alpha = 0.05, # FDR threshold
                           contrast = c("day", 3, 0), # compare day 3 to day 0
                           lfcThreshold = log2(1.5)) %>%  # 1.5 absolute FC threshold
         data.frame() %>% 
         rownames_to_column("ensgene") %>% 
         as_tibble() %>% 
         na.omit() %>% 
         left_join(df_cpm_h, by = "ensgene") %>% 
         
         filter(padj < 0.05 & abs(log2FoldChange) > log2(1.5))
      

      
      
      
      
      
      
      
      
      
      
      
   
            
            
            
            
            
            
               
      
      
      

