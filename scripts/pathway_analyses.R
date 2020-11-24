



# libraries ---------------------------------------------------------------

      library(tidyverse)
      library(purrr)
      library(ggrepel)


# data objects ------------------------------------------------------------

      load("objects/analysis_data.Rdata")
      load("objects/degs.Rdata")
      

# Pathway functions ----------------------------------------------------------

      
      source("functions/topGO.R")
      source("functions/reactome.R")
      


# functions to create pathway tibble ---------------------------------------

      # input tibble/data.frame of DEGs, where ensgene is an available column
      
      # topGO
      topGO_analysis <- function(df) {
            
            df$ensgene %>% 
                  
                  topGO() %>% 
                  as_tibble() %>% 
                  filter(weightFisher <= 0.05) %>% 
                  dplyr::rename("Pathway" = "Term", "Number of genes in set" = "Annotated", "Number of DEGs" = "Significant", "p-value" = "weightFisher")
            

      }
      
      # Reactome
      
      reactome_analysis <- function(df) {
            
            df %>% 
                  
                  reactome() %>% 
                  dplyr::rename("Pathway" = "path_name", "Number of genes in set" = "universe_n", "Number of DEGs" = "deg_n", "FDR" = "padj")
                  

      }
                  
# topGO pathway analyses --------------------------------------------------------

      
      # Pairwise
      l_pairwise_gopathways_l <- c(1:length(l_pairwise_degs_l)) %>% map(~topGO_analysis(l_pairwise_degs_l[[.x]]))
      names(l_pairwise_gopathways_l) <- names(l_pairwise_degs_l)
      
      l_pairwise_gopathways_h <- c(1:length(l_pairwise_degs_h)) %>% map(~topGO_analysis(l_pairwise_degs_h[[.x]]))
      names(l_pairwise_gopathways_h) <- names(l_pairwise_degs_h[1:length(l_pairwise_degs_h)])
      
      # All
      df_all_gopathways_l <- topGO_analysis(df_all_degs_l)
      
      df_all_gopathways_h <- topGO_analysis(df_all_degs_h)
      
      df_all_gopathways <- topGO_analysis(df_all_degs)


# reactome pathway analyses -----------------------------------------------


      # Pairwise
      l_pairwise_rpathways_l <- c(1:length(l_pairwise_degs_l)) %>% map(~reactome_analysis(l_pairwise_degs_l[[.x]]))
      names(l_pairwise_rpathways_l) <- names(l_pairwise_degs_l)
      
      l_pairwise_rpathways_h <- c(1:length(l_pairwise_degs_h)) %>% map(~reactome_analysis(l_pairwise_degs_h[[.x]]))
      names(l_pairwise_rpathways_h) <- names(l_pairwise_degs_h[1:length(l_pairwise_degs_h)])
      
      # All
      df_all_rpathways_l <- reactome_analysis(df_all_degs_l)
      
      df_all_rpathways_h <- reactome_analysis(df_all_degs_h)
      
      df_all_rpathways <- reactome_analysis(df_all_degs)
      

# pathway plots -----------------------------------------------------------

      
      # combine paths
      l_gorpaths_l <- c(1:length(l_pairwise_gopathways_l)) %>% purrr::map(~
                                                                             rbind(l_pairwise_gopathways_l[[.x]] %>% filter(!is.na(Pathway)) %>% .[1:10,] %>% 
                                                                                      select(-c(GO.ID, `p-value`)) %>% mutate(database = "GO"), 
                                                                                   l_pairwise_rpathways_l[[.x]] %>% filter(!is.na(Pathway)) %>% .[1:10,] %>% 
                                                                                      select(-c(path_id, FDR)) %>% mutate(database = "Reactome")) 
      )
      names(l_gorpaths_l) <- names(l_pairwise_gopathways_l)
      
      
      l_gorpaths_h <- c(1:length(l_pairwise_gopathways_h)) %>% purrr::map(~
                                                             rbind(l_pairwise_gopathways_h[[.x]] %>% filter(!is.na(Pathway)) %>% .[1:10,] %>% 
                                                                      select(-c(GO.ID, `p-value`)) %>% mutate(database = "GO"), 
                                                                   l_pairwise_rpathways_h[[.x]] %>% filter(!is.na(Pathway)) %>% .[1:10,] %>% 
                                                                      select(-c(path_id, FDR)) %>% mutate(database = "Reactome")) 
      )
      names(l_gorpaths_h) <- names(l_pairwise_gopathways_h)
     
      
      source("functions/pathways_plot.R")
      
      pathways_plot(l_gorpaths_h[1])

            
# save object -------------------------------------------------------------
      
      save(l_pairwise_gopathways_l, l_pairwise_gopathways_h, 
           l_pairwise_rpathways_l, l_pairwise_rpathways_h,
           l_gorpaths_l, l_gorpaths_h, 
           file = "objects/pathways.Rdata")
      
      load("objects/pathways.Rdata")
      
      
      
      
      
      
      
      
      
      
