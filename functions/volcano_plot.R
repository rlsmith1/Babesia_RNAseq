
# Function: returns volcano plot to visualize padj and LFC of genes between specified days 
#     (use all_results.R to generate results table)


# libraries ---------------------------------------------------------------

   library(tidyverse)
   library(ggrepel)


# define function ---------------------------------------------------------


   volcano_plot <- function(element) {
         
         df <- element[[1]]
      
         #add threshold column for padj < 0.05, color by padj and LFC
         df_res_volc <- df %>% 
               mutate(
                     
                     my_col =
                           
                           case_when(
                                 
                                 #-log10(FDR) < -log10(0.05) ~ "Gray",
                                 logFC > -1 & logFC < 1 ~ "Gray",
                                 logFC < 0  ~ "#4575B4",
                                 logFC > 0  ~ "#D73027",
                                 
                                 
                           )
               )
         
         # list of labeled genes
         l_top_genes <- df %>% filter(FDR <= 0.05 & !grepl("ENSCAFG", geneID)) %>% arrange(desc(abs(logFC))) %>% head(20) %>% .$geneID
         
         # plot
         df_res_volc %>% ggplot(aes(x = logFC, y = -log(FDR))) +
            geom_point(aes(color = my_col), shape = 1, alpha = 0.5) +
            scale_color_manual(values = c("#4575B4", "#D73027",  "Grey")) +
            geom_vline(aes(xintercept = log2(2)), lty = 2, color = "black") +
            geom_vline(aes(xintercept = -log2(2)), lty = 2, color = "black") +
            geom_hline(aes(yintercept = -log10(0.05)), lty = 2, color = "black") +
            geom_text_repel(data = df_res_volc %>% filter(geneID %in% l_top_genes),
                            aes(label = geneID), color = "black", size = 6) +
            xlim(-5, 10) +
            ylim(0, 20) +
            ylab("-log10(FDR)") +
            xlab("log2(Fold Change)") +
            ggtitle(element %>% names()) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.title = element_text(size = 20),
                  axis.text = element_text(size = 17),
                  plot.title = element_text(size = 23))

   }






