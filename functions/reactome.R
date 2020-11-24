
# Function that returns results table after performing Reactome pathway analysis on a given list of genes.
# input is tibble or data.frame of DEGs


# libraries ---------------------------------------------------------------


   library(biomaRt)
   library(reactome.db)
   library(tidyverse)



# load data ---------------------------------------------------------------

   load("objects/analysis_data.Rdata")
   load("objects/degs.Rdata")



# define function ---------------------------------------------------------


reactome <- function(df) {
   
   # connect to ensembl genome
   ensembl <- useMart("ensembl")
   dog_genome <- useDataset("clfamiliaris_gene_ensembl", mart = ensembl)
   
   # establish gene universe
   l_gene_universe <- df_rawcounts_filtered[[1]]
   
   # connect ensembl id to reactome pathway id
   convert_id_universe <- getBM(attributes = c("ensembl_gene_id", "reactome"), 
                                filters = "ensembl_gene_id",
                                values = l_gene_universe,
                                mart = dog_genome)
   
   convert_id_degs <- getBM(attributes = c("ensembl_gene_id", "reactome"), 
                       filters = "ensembl_gene_id",
                       values = df$ensgene,
                       mart = dog_genome)
   
   # count number of genes in each pathway
   df_reactome_id_universe <- convert_id_universe %>% as_tibble() %>% dplyr::rename("path_id" = "reactome", "ensgene" = "ensembl_gene_id")
   df_reactome_id_degs <- convert_id_degs %>% as_tibble() %>% dplyr::rename("path_id" = "reactome", "ensgene" = "ensembl_gene_id")

   df_universe_counts <- df_reactome_id_universe %>% dplyr::count(path_id) %>% filter(path_id != "") %>% dplyr::rename("universe_n" = "n")
   df_deg_counts <- df_reactome_id_degs %>% dplyr::count(path_id) %>% filter(path_id != "") %>% dplyr::rename("deg_n" = "n")

   df_all_counts <- df_universe_counts %>% left_join(df_deg_counts)

   # hypergeometric distribution
   df_enrich_path <- df_all_counts %>%
      mutate(p = dhyper(deg_n, universe_n, length(l_gene_universe), nrow(df))) %>%
      arrange(p) %>%
      mutate(padj = p.adjust(p, method = "BH")) %>%
      filter(padj < 0.05)

   # align IDs with pathway names
   df_path_id2name <- as.data.frame(reactomePATHID2NAME) %>% 
      as_tibble() %>%
      filter(grepl("R-CFA", DB_ID)) %>% 
      mutate(path_name = substr(path_name, 19, nchar(path_name))) %>% 
      dplyr::rename("path_id" = "DB_ID")
   
   df_enrich_path <- df_enrich_path %>% 
      left_join(df_path_id2name) %>% 
      dplyr::select("path_id", "path_name", "universe_n", "deg_n", "padj")
   
   # return results
   return(df_enrich_path)
   

}







# # add names of genes in pathway -------------------------------------------
# 
# 
# 
#    #high
#    df_gene_in_path_high <- df_high_path_names %>% left_join(df_reactome_id_high) %>% left_join(df_geneID)
#    
#    
#    l_genes_in_path_high <- unique(df_gene_in_path_high$path_name) %>% 
#          purrr::map(~df_gene_in_path_high %>% filter(path_name == .x) %>% pull(geneID) %>% paste0(collapse = ","))
#    
#    names(l_genes_in_path_high) <- unique(df_gene_in_path_high$path_name)
#    
#    df_genes_by_path_high <- bind_cols(l_genes_in_path_high) %>% pivot_longer(cols = 1:length(l_genes_in_path_high), names_to = "path_name", values_to = "geneID")
#    
#    
#    df_rpaths_high <- df_high_path_names %>% 
#       left_join(df_genes_by_path_high) %>% 
#       dplyr::rename("Reactome.ID" = "path_id", "Biological pathway" = "path_name", "Total genes in set" = "universe_n", "Number of DEGs" = "high_n")%>% 
#       mutate(padj = round(padj, 5))
# 
# 
#    #low
#    df_gene_in_path_low <- df_low_path_names %>% left_join(df_reactome_id_low) %>% left_join(df_geneID)
#    
#    
#    l_genes_in_path_low <- unique(df_gene_in_path_low$path_name) %>% 
#          purrr::map(~df_gene_in_path_low %>% filter(path_name == .x) %>% pull(geneID) %>% paste0(collapse = ","))
#    
#    names(l_genes_in_path_low) <- unique(df_gene_in_path_low$path_name)
#    
#    df_genes_by_path_low <- bind_cols(l_genes_in_path_low) %>% pivot_longer(cols = 1:length(l_genes_in_path_low), names_to = "path_name", values_to = "geneID")
#    
#    df_rpaths_low <- df_low_path_names %>% 
#       left_join(df_genes_by_path_low) %>% 
#       dplyr::rename("Reactome.ID" = "path_id", "Biological pathway" = "path_name", "Total genes in set" = "universe_n", "Number of DEGs" = "low_n") %>% 
#       mutate(padj = round(padj, 5))
# 
# 




  




