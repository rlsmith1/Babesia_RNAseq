


# Libraries ---------------------------------------------------------------

   library(fastcluster)
   library(purrr)
   library(cluster)
   library(corrr)
   library(tidyverse)
   library(data.table)
   library(plotly)



# data --------------------------------------------------------------------

   load("objects/analysis_data.Rdata")   



# Calculate L2FC on all genes ---------------------------------------------

   # Structure
   
      # pivot and add metadata
      df_cpm_pivot_metadata <- df_cpm_filtered[3:30] %>% 
         cbind(c(df_geneID)) %>% 
         as_tibble() %>% 
         pivot_longer(cols = 1:28, values_to = "CPM", names_to = "sample") %>% 
         left_join(df_metadata %>% rownames_to_column(var = "sample"), by = "sample") %>% 
         ungroup()
      
      # create nested df
      df_nest <- df_cpm_pivot_metadata %>% 
         filter(dog != 0) %>% 
         mutate(ensgene_dog = paste0(ensgene, sep = "_", dog)) %>% 
         nest_by(ensgene_dog)
   
   # calculate L2FC for each gene in high and low groups each
   l_l2fc <- 1:length(df_nest$data) %>% map(~df_nest$data[[.x]] %>% 
                                               mutate(L2FC = log2(CPM/.$CPM[.$day==0] + .0000000000000001)))
   
   # create df of L2FC from day 0 in each dog
   df_l2fc_filtered <- rbindlist(l_l2fc) %>% 
      as_tibble() %>% 
      pivot_wider(id_cols = c(ensgene, geneID), values_from = L2FC, names_from = sample)


   

# Separate high and low inoculum ------------------------------------------


   # pivot and add metadata
   df_l2fc_pivot_metadata <- df_l2fc_filtered %>% 
      
      pivot_longer(cols = 3:29, names_to = "sample", values_to = "L2FC") %>% 
      dplyr::select(-"geneID") %>% 
      left_join(df_metadata %>% rownames_to_column("sample"))
   
   # filter for low vs high
   df_l2fc_pivot_l <- df_l2fc_pivot_metadata %>% filter(dose == "low")
   df_l2fc_pivot_h <- df_l2fc_pivot_metadata %>% filter(dose == "high")
   
   # pivot wider to structure for correlation
   
   # low
   df_l2fc_by_sample_l <- df_l2fc_pivot_l %>% 
      pivot_wider(id_cols = sample, names_from = ensgene, values_from = L2FC) %>% 
      dplyr::select(-sample)
   
   # high
   df_l2fc_by_sample_h <- df_l2fc_pivot_h %>% 
      pivot_wider(id_cols = sample, names_from = ensgene, values_from = L2FC) 
   
   # find columns (genes) with 0 stdev
   
   # calculate standard deviation of each column
   l_row_sds <- df_l2fc_by_sample_h %>% select(-sample) %>% as.matrix() %>% colSds()
   l_row_sds <- c("sd", l_row_sds)
   df_l2fc_by_sample_h_sd <- df_l2fc_by_sample_h %>% rbind(l_row_sds)
   df_l2fc_by_sample_h_sd <- df_l2fc_by_sample_h_sd %>%
      mutate_at(vars(c(2:8002)), as.numeric) %>%
      as_tibble()
   
   # where does sd equal zero?
   df_l2fc_by_sample_h_sd[16,]
   
   df_l2fc_by_sample_h <- df_l2fc_by_sample_h %>% select(-sample)


# Structure L2FC data (for plotting) -----------------------------------------------------


   # pivot and add metadata
   df_l2fc_pivot <- df_l2fc_filtered %>% 
      pivot_longer(cols = 3:29, names_to = "sample", values_to = "L2FC") %>% 
      dplyr::select(-"geneID") %>% 
      left_join(df_metadata %>% rownames_to_column("sample"))
   
   # take mean CPM of each gene for dose_day and pivot wider
   df_avg_l2fc <- df_l2fc_pivot %>% 
      group_by(ensgene, dose_day) %>% 
      summarise(mean_L2FC = mean(L2FC, na.rm = TRUE)) %>% 
      
      pivot_wider(id_cols = ensgene, names_from = dose_day, values_from = mean_L2FC) 




# Cluster all genes! ------------------------------------------------------


   # LOW
   
      # correlate genes based on expression values on each day
      df_corr_l <- correlate(df_l2fc_by_sample_l, diagonal = 1)
      
      # convert to a matrix
      m_corr_l <- as.matrix(df_corr_l %>% dplyr::select(-rowname))
      
      # calculate distances based on correlation values
      dist_l <- dist(m_corr_l)
      
      # cluster
      hclust_l <- hclust(dist_l, method = "ward.D2") 
      
      # optimize k based on average width
      sil_width_l <- map_dbl(2:10, function(k) {
         model <- pam(dist_l, k = k)
         model$silinfo$avg.width
      })
      
      df_sil_l <- tibble(k = 2:10, sil_width = sil_width_l)
      
      df_sil_l %>% ggplot(aes(x = k, y = sil_width)) +
         geom_point() +
         geom_line()
      
      # Determine which observations belong to which cluster    
      my_clusters_l <- cutree(tree = hclust_l, k = 3) 
      
      df_corr_l$cluster <- my_clusters_l
      df_clust_l <- df_corr_l %>% dplyr::select(rowname, cluster)
      names(df_clust_l)[1] = "ensgene"
      
      # plot
      
         # calculate number of genes in each cluster
         df_n_l <- df_clust_l %>% dplyr::count(cluster)
         df_clust_l <- df_clust_l %>% left_join(df_n_l, by = "cluster")
         
         # add cluster to low genes
         df_avg_l2fc_l <- df_avg_l2fc[, c(1,grep("low_*", colnames(df_avg_l2fc)))] %>% left_join(df_clust_l, by = "ensgene")
         
         # pivot
         df_avg_l2fc_pivot_l <- df_avg_l2fc_l %>% pivot_longer(cols = 2:7, names_to = "day", values_to = "L2FC") %>% 
            mutate(day = substr(day, start = 5, stop = 5))
         
         df_clust_day_avgl2fc_l <- df_avg_l2fc_pivot_l %>%
            
            mutate(n = as.numeric(n)) %>% 
            
            group_by(cluster, n, day) %>%
            
            summarise(mean_L2FC = mean(L2FC),
                      stdev_L2FC = sd(L2FC)) %>% 
            
            mutate(error = qnorm(0.975)*stdev_L2FC/sqrt(n)) %>% 
            mutate(left = mean_L2FC - error, right = mean_L2FC + error) %>% 
            
            mutate(day = as.numeric(day)) %>% 
            mutate(cluster = as.factor(cluster)) %>% 
            
            mutate(cluster_n = paste0(cluster, sep = ", n = ", n))
         
         # plot all clusters
         df_clust_day_avgl2fc_l %>% ggplot() + 
            geom_line(aes(x = day, y = mean_L2FC, col = cluster)) + 
            geom_errorbar(aes(x = day, ymin = left, ymax = right), width = 0.2) +
            facet_wrap(~cluster_n) +
            theme_minimal() +
            ggtitle("low inoculum clusters")


   # HIGH
   
      # correlate genes based on expression values on each day
      df_corr_h <- correlate(df_l2fc_by_sample_h, diagonal = 1)
      
      # convert to a matrix
      m_corr_h <- as.matrix(df_corr_h %>% dplyr::select(-rowname))
      
      # calculate distances based on correlation values
      dist_h <- dist(m_corr_h)
      
      # cluster
      hclust_h <- hclust(dist_h, method = "ward.D2") 
      
      # optimize k based on average width
      sil_width_h <- map_dbl(2:10, function(k) {
         model <- pam(dist_h, k = k)
         model$silinfo$avg.width
      })
      
      df_sil_h <- tibble(k = 2:10, sil_width = sil_width_h)
      
      df_sil_h %>% ggplot(aes(x = k, y = sil_width)) +
         geom_point() +
         geom_line()
      
      # Determine which observations belong to which cluster    
      my_clusters_h <- cutree(tree = hclust_h, k = 3) 
      
      df_corr_h$cluster <- my_clusters_h
      df_clust_h <- df_corr_h %>% dplyr::select(rowname, cluster)
      names(df_clust_h)[1] = "ensgene"
      
      # plot
      
         # calculate number of genes in each cluster
         df_n_h <- df_clust_h %>% dplyr::count(cluster)
         df_clust_h <- df_clust_h %>% left_join(df_n_h, by = "cluster")
         
         # add cluster to low genes
         df_avg_l2fc_h <- df_avg_l2fc[, c(1,grep("high_*", colnames(df_avg_l2fc)))] %>% left_join(df_clust_h, by = "ensgene")
         
         # pivot
         df_avg_l2fc_pivot_h <- df_avg_l2fc_h %>% pivot_longer(cols = 2:7, names_to = "day", values_to = "L2FC") %>% 
            mutate(day = substr(day, start = 6, stop = 6))
         
         df_clust_day_avgl2fc_h <- df_avg_l2fc_pivot_h %>%
            
            mutate(n = as.numeric(n)) %>% 
            
            group_by(cluster, n, day) %>%
            
            summarise(mean_L2FC = mean(L2FC),
                      stdev_L2FC = sd(L2FC)) %>% 
            
            mutate(error = qnorm(0.975)*stdev_L2FC/sqrt(n)) %>% 
            mutate(left = mean_L2FC - error, right = mean_L2FC + error) %>% 
            
            mutate(day = as.numeric(day)) %>% 
            mutate(cluster = as.factor(cluster)) %>% 
            
            mutate(cluster_n = paste0(cluster, sep = ", n = ", n))
         
         
         # plot all clusters
         df_clust_day_avgl2fc_h %>% ggplot() + 
            geom_line(aes(x = day, y = mean_L2FC, col = cluster)) + 
            geom_errorbar(aes(x = day, ymin = left, ymax = right), width = 0.2) +
            facet_wrap(~cluster_n) +
            theme_minimal() +
            ggtitle("high inoculum clusters")




# save objects ------------------------------------------------------------


   save(df_l2fc_filtered, df_avg_l2fc,
        df_corr_l, m_corr_l, dist_l, hclust_l, df_sil_l, 
        df_corr_h, m_corr_h, dist_h, hclust_h, df_sil_h, file = "objects/cluster_all_objects.Rdata")
         
   save(df_l2fc_filtered, df_avg_l2fc, df_corr_l, hclust_l, df_corr_h, hclust_h, df_clust_h, df_clust_l,
        file = "objects/cluster_plots.Rdata")
   
   load("objects/cluster_all_objects.Rdata")
   load("objects/cluster_plots.Rdata")




















