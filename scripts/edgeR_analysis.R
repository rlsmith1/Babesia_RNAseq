


# libraries ---------------------------------------------------------------

      library(edgeR)
      library(tidyverse)
      library(magrittr)
      library(pheatmap)
      library(RColorBrewer)
      library(plotly)
      library(stringi)
      library(purrr)
      library(ggrepel)
      library(matrixStats)



# data --------------------------------------------------------------------

      load("objects/analysis_data.Rdata")




# Major sources of variation in dataset -----------------------------------

      # Mean-variance relationship
      mean_counts <-  df_rawcounts_filtered[3:30] %>% rowMeans() 
      var_counts <- df_rawcounts_filtered[3:30] %>% as.matrix() %>% rowVars() 
      
      df_mean_var <- data.frame(mean_counts, var_counts)
      
      ggplot(df_mean_var) +
         
         geom_point(aes(mean_counts, var_counts)) +
         scale_y_log10() +
         scale_x_log10() +
         xlab("Mean counts per gene") +
         ylab("Variance per gene") # we see an increase in variance associated with increase in mean counts

      # Variance stabilizing transformation to visualize counts
      m_rawcounts_filtered_trans <- df_rawcounts_filtered[3:30] %>% as.matrix()
      rownames(m_rawcounts_filtered_trans) <- df_rawcounts_filtered$ensgene
      m_rawcounts_filtered_trans <- log(m_rawcounts_filtered_trans + 1) # natural log, add constant to deal with 0s
      
      # Relabel columns
      df_metadata1 <- df_metadata %>% mutate(day = paste0("D",day))
      
      df_metadata2 <- df_metadata1 %>% mutate(dose_spl = paste0(dose, 
                                                                sep = ifelse(dose == "low",
                                                                             "__",
                                                                             "_"), 
                                                                case_when(
                                                                   
                                                                   dog == 0 ~ "0",
                                                                   dog == 1 ~ "1",
                                                                   dog == 2 ~ "2", 
                                                                   dog == 3 ~ "1",
                                                                   dog == 4 ~ "2", 
                                                                   dog == 5 ~ "3"
                                                                   
                                                                )))
      
      df_metadata3 <- df_metadata2 %>% mutate(label = paste0(dose_spl, sep = "_", day))
                                                
      colnames(m_rawcounts_filtered_trans) <- df_metadata3$label
      
      # Visualize transformed counts
      
            # Heatmap
      
                  # Correlate samples based on transformed counts
                  m_cor_samples <- cor(m_rawcounts_filtered_trans)
                  
                  # Plot
                  pheatmap(m_cor_samples)
                  
                  # Mean correlation of highly correlated samples
                  
                        # Day 0 + low day 1
                        m_cor_samples %>% data.frame() %>% rownames_to_column("label") %>% as_tibble() %>% 
                           pivot_longer(2:29, names_to = "label1", values_to = "corr") %>% 
                           filter(label != label1) %>% 
                           
                           # first filter
                           left_join(df_metadata3, by = "label") %>% 
                           filter(day == "D0" | dose_day == "low_1") %>% 
                           
                           # second filter
                           dplyr::select(1:3) %>% 
                           dplyr::rename("label2" = "label", "label" = "label1") %>% 
                           left_join(df_metadata3, by = "label") %>% 
                           filter(day == "D0" | dose_day == "low_1") %>% 
                           
                           # take mean
                           summarise(mean(corr))
                           
                        # High day 1 + low day 3
                        m_cor_samples %>% data.frame() %>% rownames_to_column("label") %>% as_tibble() %>% 
                           pivot_longer(2:29, names_to = "label1", values_to = "corr") %>% 
                           filter(label != label1) %>% 
                           
                           # first filter
                           left_join(df_metadata3, by = "label") %>% 
                           filter(dose_day == "high_1" | dose_day == "low_3") %>% 
                           
                           # second filter
                           dplyr::select(1:3) %>% 
                           dplyr::rename("label2" = "label", "label" = "label1") %>% 
                           left_join(df_metadata3, by = "label") %>% 
                           filter(dose_day == "high_1" | dose_day == "low_3") %>% 
                           
                           # take mean
                           summarise(mean(corr))

            # PCA
                  
                  # Calculate PCA data
                  l_pca_data <- prcomp(m_rawcounts_filtered_trans)
                  df_pca_data <- l_pca_data$rotation %>% data.frame() %>% rownames_to_column("label") %>% as_tibble()
                  df_pca_meta_data <- df_metadata3 %>% left_join(df_pca_data, by = "label")
                  
                  round(summary(l_pca_data)$importance[2,1]*100, 2)
               
                  # plot
                  df_pca_meta_data %>% ggplot(aes(x = PC1, y = PC2)) +
                     
                     geom_point(aes(color = dose), size = 3) +
                     geom_text_repel(aes(label = substr(day, start = 2, stop = 2)), color = "black") +
                     
                     labs(color = "inoculum") +
                     
                     xlab(paste0("PC1: ", round(summary(l_pca_data)$importance[2,1]*100, 2), "% of variance")) +
                     ylab(paste0("PC2: ", round(summary(l_pca_data)$importance[2,2]*100, 2), "% of variance")) +
                     
                     theme_bw()
                     

            
# format data for edgeR -------------------------------------------------------------


      #convert to matrix
      m_rawcounts_l <- df_rawcounts_l %>% dplyr::select(-c("geneID", "ensgene")) %>% as.matrix()
      rownames(m_rawcounts_l) <- df_rawcounts_l$ensgene
      
      m_rawcounts_h <- df_rawcounts_h %>% dplyr::select(-c("geneID", "ensgene")) %>% as.matrix()
      rownames(m_rawcounts_h) <- df_rawcounts_h$ensgene

      # specify group
      
            # group 1 = day 0
            # group 2 = day 1
            # group 3 = day 3
            # group 4 = day 4
            # group 5 = day 6
            # group 6 = day 8
      
      group_l <- factor(c(1, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6))
      
      group_h <- factor(c(1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5, 6, 1, 2, 3))
      
      # Create DGEList
      l_dge_l <- DGEList(counts = m_rawcounts_l, group = group_l)

      l_dge_h <- DGEList(counts = m_rawcounts_h, group = group_h)

      # Already filtered out lowly expressed genes, don't need to use edgeR function
      
      # Create design matrix
      m_design_l <- model.matrix(~group_l)
      
      m_design_h <- model.matrix(~group_h)
      
      # Effective library size
      l_dge_l <- calcNormFactors(l_dge_l) # scales lib size to minimize LFC between samples for most genes
      
      l_dge_h <- calcNormFactors(l_dge_h)
      
      # Estimate dispersions ***is this a single factor experiment????***
      l_dge_l <- estimateDisp(l_dge_l, m_design_l)
      
      l_dge_h <- estimateDisp(l_dge_h, m_design_h)



# Testing for DEGs --------------------------------------------------------

      
      # General test
      glm_l <- glmQLFit(l_dge_l, m_design_l)

      glm_h <- glmQLFit(l_dge_h, m_design_h)
      
      # Write function to format DEG results
      format_deg_results <- function(DGELRT) {
      
            
            DGELRT$table %>% 
                  rownames_to_column("ensgene") %>% 
                  as_tibble() %>% 
                  mutate(FDR = p.adjust(.$PValue, method = "BH")) %>% 
                  arrange(FDR) %>% 
                  #filter(FDR <= 0.05) %>% 
                  left_join(df_geneID, by = "ensgene") %>% 
                  dplyr::select(c(ensgene, geneID), everything())
            
      }
      
      # all results

         # With a L2FC threshold (can't use all coefficients in one)
         l_pairwise_l <- c(2:6) %>% purrr::map(~glmQLFTest(glm_l, coef = .x) %>% format_deg_results())
         names(l_pairwise_l) <- c("Day 1 vs Day 0", "Day 3 vs Day 0", "Day 4 vs Day 0", "Day 6 vs Day 0", "Day 8 vs Day 0")
         
         l_pairwise_h <- c(2:6) %>% purrr::map(~glmQLFTest(glm_h, coef = .x) %>% format_deg_results())
         names(l_pairwise_h) <- c("Day 1 vs Day 0", "Day 3 vs Day 0", "Day 4 vs Day 0", "Day 6 vs Day 0", "Day 8 vs Day 0")
         
      # DEGs (FDR <= 0.05)
         
         # With a L2FC threshold (can't use all coefficients in one)
         l_pairwise_degs_l <- c(1:length(l_pairwise_l)) %>% purrr::map(~l_pairwise_l[[.x]] %>% 
                                                                          filter(FDR <= 0.05 & abs(logFC) >= 1))
         names(l_pairwise_degs_l) <- names(l_pairwise_l)
         
         l_pairwise_degs_h <- c(1:length(l_pairwise_h)) %>% purrr::map(~l_pairwise_h[[.x]] %>% 
                                                                          filter(FDR <= 0.05 & abs(logFC) >= 1))
         names(l_pairwise_degs_h) <- names(l_pairwise_h)
         

         
# Find DEGs on all days ---------------------------------------------------
         
         # low
         df_all_degs_l <- l_pairwise_degs_l[[2]] %>% 
            inner_join(l_pairwise_degs_l[[3]], by = c("ensgene", "geneID")) %>% dplyr::select(c(ensgene, geneID)) %>% 
            inner_join(l_pairwise_degs_l[[4]], by = c("ensgene", "geneID")) %>% dplyr::select(c(ensgene, geneID)) %>% 
            inner_join(l_pairwise_degs_l[[5]], by = c("ensgene", "geneID")) %>% dplyr::select(c(ensgene, geneID))
         
         # high
         df_all_degs_h <- l_pairwise_degs_h[[1]] %>% 
            inner_join(l_pairwise_degs_h[[2]], by = c("ensgene", "geneID")) %>% dplyr::select(c(ensgene, geneID)) %>% 
            inner_join(l_pairwise_degs_h[[3]], by = c("ensgene", "geneID")) %>% dplyr::select(c(ensgene, geneID)) %>% 
            inner_join(l_pairwise_degs_h[[4]], by = c("ensgene", "geneID")) %>% dplyr::select(c(ensgene, geneID)) %>% 
            inner_join(l_pairwise_degs_h[[5]], by = c("ensgene", "geneID")) %>% dplyr::select(c(ensgene, geneID))
         
         # overlap
         df_all_degs <- df_all_degs_h %>% inner_join(df_all_degs_l, by = c("ensgene", "geneID"))
      

# Find clusters -----------------------------------------------------------

         df_all_degs_l <- df_all_degs_l %>% left_join(df_clust_l, by = "ensgene")
         
         df_all_degs_h <- df_all_degs_h %>% left_join(df_clust_h, by = "ensgene") 
         
         df_all_degs <- df_all_degs %>% 
            left_join(df_clust_l, by = "ensgene") %>% left_join(df_clust_h, by = "ensgene") %>% 
            dplyr::rename("low cluster" = "cluster.x", "high cluster" = "cluster.y")
         

# save DEG lists ----------------------------------------------------------


      save(l_pairwise_l, l_pairwise_h,
           l_pairwise_degs_l, l_pairwise_degs_h, 
           df_all_degs_l, df_all_degs_h, df_all_degs, file = "objects/degs.Rdata")

      load("objects/degs.Rdata")

      

# volcano plots -----------------------------------------------------------

      source("functions/volcano_plot.R")

      c(1:length(l_pairwise_l)) %>% map(~(volcano_plot(l_pairwise_l[.x])))
      
      c(1:length(l_pairwise_h)) %>% map(~(volcano_plot(l_pairwise_h[.x])))
      
      c(1:length(l_pairwise_l)) %>% map(~volcano_plot(l_pairwise_l[.x]), ~volcano_plot(l_pairwise_h[.x]))
      
      
            
      
      

