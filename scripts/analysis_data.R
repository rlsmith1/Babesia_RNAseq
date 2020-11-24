

# libraries ---------------------------------------------------------------

      library(tidyverse)



# Data --------------------------------------------------------------------

      # count data
      df_rawcounts <- read_csv("data/RawCountFiles_RSEM_genes.csv")
      
      df_rawcounts <- df_rawcounts %>% separate(symbol, into = c("ensgene", "geneID"), sep = "_") %>% as_tibble()
      
      # metadata
      df_metadata <- read_csv("data/metadata.csv") %>% mutate(dose = as.factor(dose)) %>% mutate(day = as.factor(day))



# Convert rawcounts to CPM ----------------------------------------------------------


      df_counts <- df_rawcounts[3:30]
      df_ids <- df_rawcounts[1:2]
      
      # FORMULA: CPM = readsMappedToGene * 1/totalNumReads * 10^6
      
      df_cpm <- df_counts %>% 
            apply(2, function(x){x/sum(x) * 10^6}) %>% 
            as_tibble() 
      
      df_cpm <- cbind(df_ids, df_cpm) %>% as_tibble()



# filter by CPM -----------------------------------------------------------


      # pivot expression matrix
      df_cpm_pivot <- df_cpm %>% 
            pivot_longer(cols = 3:30, names_to = "sample", values_to = "CPM") %>%
            left_join(df_metadata, by = "sample") %>% 
            dplyr::select(-c("dog"))
      
      # Filter at 1 CPM in at least two dose-day treatment groups
      df_cpm_pivot2 <- df_cpm_pivot %>% 
            arrange(ensgene, dose_day) %>% 
            mutate(thresh = ifelse(CPM >= 1, 1, 0)) %>% 
            group_by(ensgene, dose_day) %>% 
            mutate(thresh2 = sum(thresh)) %>% 
            filter(thresh2 >= 2) %>% 
            dplyr::select(-c("thresh", "thresh2"))
      
      df_cpm_filtered <- df_cpm_pivot2 %>% 
            pivot_wider(id_cols = c("ensgene", "geneID"), names_from = "sample", values_from = "CPM") %>% 
            na.omit()
      
      # also subset the rawcounts data frame using genes that meet CPM reqs for all genes 
      df_cpm_filtered_pivot <- df_cpm_pivot2 %>% ungroup() %>% na.omit()
      
      df_rc_pivot <- df_rawcounts %>% pivot_longer(cols = 3:30, names_to = "sample", values_to = "raw_counts")
      
      df_rcf_pivot <- df_cpm_filtered_pivot %>% left_join(df_rc_pivot) 
      
      df_rawcounts_filtered <- df_rcf_pivot %>% pivot_wider(id_cols = c("ensgene", "geneID"), names_from = sample, values_from = raw_counts) %>% na.omit()

      
      
# Confirm filtering by visualizing count distributions  ------------------------------------------------------
      
      
      # before filtering
      
      ggplot(df_rc_pivot) +
            geom_histogram(aes(x = log(raw_counts)), stat = "bin", bins = 200) +
            facet_wrap(~sample) +
            xlab("log(Raw counts)") +
            ylab("Number of genes") +
            title("Distribution of raw counts") 
      
      
      # after filtering
      
      ggplot(df_rcf_pivot) +
            geom_histogram(aes(x = log(raw_counts)), stat = "bin", bins = 200) +
            facet_wrap(~sample) +
            xlab("log(Raw counts)") +
            ylab("Number of genes") +
            title("Distribution of filtered raw counts") 
      
      
      
      
      # before filtering
      
      ggplot(df_cpm_pivot) +
            geom_histogram(aes(x = log(CPM)), stat = "bin", bins = 200) +
            facet_wrap(~sample) +
            geom_vline(xintercept = log(1), color = "red", lty = 2) +
            geom_vline(xintercept = log(5), color = "blue", lty = 2) +
            xlab("log(CPM)") +
            ylab("Number of genes") +
            ggtitle("Distribution of CPM") 
      
      
      # after filtering at 1 CPM
      
      ggplot(df_cpm_filtered_pivot) +
            geom_histogram(aes(x = log(CPM)), stat = "bin", bins = 200) +
            facet_wrap(~sample) +
            xlab("log(CPM)") +
            ylab("Number of genes") +
            ggtitle("Distribution of filtered CPM (1 CPM threshold)") 
      
      
      

# Order samples in data -----------------------------------------------------------


      # for rawcounts
      df_rc_counts <- df_rawcounts_filtered[3:30]
      df_rc_labels <- df_rawcounts_filtered[1:2]
      
      idx <- match(rownames(df_metadata %>% column_to_rownames("sample")), colnames(df_rc_counts))
      df_rc_counts <- df_rc_counts[ , idx]
      
      df_rawcounts_filtered <- bind_cols(df_rc_labels, df_rc_counts)
      
      
      # and for CPM
      df_cpm_counts <- df_cpm_filtered[3:30]
      df_cpm_labels <- df_cpm_filtered[1:2]
      
      idx <- match(rownames(df_metadata %>% column_to_rownames("sample")), colnames(df_cpm_counts))
      df_cpm_counts <- df_cpm_counts[ , idx]
      
      df_cpm_filtered <- bind_cols(df_cpm_labels, df_cpm_counts)
      
      
      
# split into low vs high --------------------------------------------------
      
      
      df_rawcounts_filtered_pivot <- df_rawcounts_filtered %>% 
            
            pivot_longer(cols = 3:ncol(df_rawcounts_filtered), names_to = "sample", values_to = "rawcounts") %>% 
            left_join(df_metadata, by = "sample")
      
      df_cpm_filtered_pivot <- df_cpm_filtered %>% 
            
            pivot_longer(cols = 3:ncol(df_cpm_filtered), names_to = "sample", values_to = "CPM") %>% 
            left_join(df_metadata, by = "sample")
      
      
      # Low
      df_rawcounts_l <- df_rawcounts_filtered_pivot %>% 
            
            filter(dose == "low") %>% 
            pivot_wider(id_cols = c("ensgene", "geneID"), names_from = "sample", values_from = "rawcounts")
      
      df_cpm_l <- df_cpm_filtered_pivot %>% 
            
            filter(dose == "low") %>% 
            pivot_wider(id_cols = c("ensgene", "geneID"), names_from = "sample", values_from = "CPM")
      
      
      # High
      df_rawcounts_h <- df_rawcounts_filtered_pivot %>% 
            
            filter(dose == "high") %>% 
            pivot_wider(id_cols = c("ensgene", "geneID"), names_from = "sample", values_from = "rawcounts")
      
      
      df_cpm_h <- df_cpm_filtered_pivot %>% 
            
            filter(dose == "high") %>% 
            pivot_wider(id_cols = c("ensgene", "geneID"), names_from = "sample", values_from = "CPM")
      
      

# create df_geneID --------------------------------------------------------


      df_geneID <- df_cpm_filtered %>% select(c(ensgene, geneID))
      
      # fill in names of unknown ENSCAFG (not in database)
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000000500"] <- "DLA-64"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000000819"] <- "DLA-DOB"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000000896"] <- "DLA-DOA"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000032222"] <- "DLA-12"
      
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000010511"] <- "OASL1"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000024447"] <- "OASL2"
      
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000015901"] <- "FTH1a"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000030465"] <- "FTH1b"
      
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000011638"] <- "SPTA1"
      
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000014723"] <- "UBA52a"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000018652"] <- "UBA52b"
      
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000020220"] <- "HP"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000028980"] <- "SIGLEC14"
      
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000028569"] <- "HBZ"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000030286"] <- "HBB"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000032615"] <- "HBA"
      
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000000415"] <- "RAP1B"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000000641"] <- "HSPA1A"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000000340"] <- "CMC2"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000008154"] <- "FER-1 LIKE PROTEIN 4"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000000061"] <- "BLOC1S1"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000005637"] <- "MMADHC"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000006946"] <- "CSNK2A1"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000000121"] <- "STAT2"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000029624"] <- "IGHZ"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000015217"] <- "TUBA1"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000016263"] <- "TUBA1A"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000031588"] <- "SEC62"
      
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000031706"] <- "FABP5"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000012657"] <- "IRGM (IFI1)"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000014312"] <- "ITGAX"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000029509"] <- "MCEMP1"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000013769"] <- "SQOR"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000025128"] <- "TRAC"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000029830"] <- "APOL2"
      df_geneID$geneID[df_geneID$ensgene == "ENSCAFG00000014478"] <- "TRBV28"
      
      
      
      
      df_geneID %>% filter(grepl("ENSCAFG", geneID))
      
      
      
      

# save data structures for future analyses --------------------------------


      save(df_rawcounts_filtered, df_cpm_filtered, df_metadata, df_geneID, 
           df_rawcounts_l, df_rawcounts_h, df_cpm_l, df_cpm_h, file = "objects/analysis_data.Rdata")
      
    
      
      
      
      
      
      
      
      
      
      
      




