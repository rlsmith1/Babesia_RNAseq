
# Function that returns results table after performing topGO pathway analysis on a given list of genes.
# input is list of ensembl gene IDs



# load libraries and objects ----------------------------------------------

load("objects/analysis_data.Rdata")
load("objects/degs.Rdata")

library(biomaRt)
library(topGO)
library(GO.db)
library(Rgraphviz)
library(tidyverse)
library(purrr)
library(magrittr)


# define function ---------------------------------------------------------



topGO <- function(gene_list) {
      
      #gene universe
      l_gene_universe <- df_rawcounts_filtered[[1]]
      
      #create GO annotation
      
            #connect to ensembl
            ensembl <- useMart("ensembl")
            dog_genome <- useDataset("clfamiliaris_gene_ensembl", mart = ensembl)
            
            #map each gene to GO terms
            m_go_ids <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "go_id", "namespace_1003"), 
                              filters = "ensembl_gene_id",
                              values = l_gene_universe,
                              mart = dog_genome)
            
            #build annotation list
            l_gene_2_GO <- unstack(m_go_ids[,c(1,3)])
            
      #remove DEGs without GO annotation
      keep <- gene_list %in% m_go_ids[,1]
      keep <- which(keep == TRUE)
      gene_list <- gene_list[keep]
      
      # make named factor showing DEGs
      
      l_geneList <- factor(as.integer(l_gene_universe %in% gene_list))
      names(l_geneList) <-  l_gene_universe
      
      #create topGO data object
      GOdata <- new("topGOdata",
                    ontology = "BP",
                    allGenes = l_geneList,
                    annotationFun = annFUN.GO2genes,
                    GO2genes = l_gene_2_GO,
                    description = "GO analysis of cluster")
      
      #test for significance
      weight_fisher_result <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
      
      #generate results table
      allGO <- usedGO(GOdata)
      all_res <- GenTable(GOdata, 
                          weightFisher = weight_fisher_result,
                          orderBy = "weightFisher",
                          topNodes = length(allGO))
      
      #fix names so they don't get cut off
         
         #pull list of all GO terms from database
         l_go_terms <- as.list(GOTERM)
         
         #extract GO ID and term name
         l_go_ids <- 1:length(l_go_terms) %>% map(~GOID(l_go_terms[[.x]])) %>% unlist()
         l_go_terms <- 1:length(l_go_terms) %>% map(~Term(l_go_terms[[.x]])) %>% unlist()
         
         df_go_terms <- tibble(GO.ID = l_go_ids, Term = l_go_terms)
         
         #merge with results df
         all_res %<>% dplyr::select(-Term) %>% left_join(df_go_terms, by = "GO.ID")
      
         #sort columns
         all_res %<>% dplyr::select(GO.ID, Term, Annotated, Significant, weightFisher)
      
      #return results table
      return(all_res)
      
      
      
}











