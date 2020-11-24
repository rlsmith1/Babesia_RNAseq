


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
      


# data --------------------------------------------------------------------

      load("objects/analysis_data.Rdata")



# format data -------------------------------------------------------------


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

      # Already filtered out lowly expressed genes, don't need to use their function
      
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
      
      
            
      
      

