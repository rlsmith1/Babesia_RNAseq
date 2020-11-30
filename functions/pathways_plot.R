

# input combined pathways list (includes both GO ane reactome pathways)

# libraries ---------------------------------------------------------------

      library(tidyverse)
      library(RColorBrewer)

# write function ----------------------------------------------------------


      pathways_plot <- function(object) {
            
            object[[1]] %>% 
                  
                  ggplot(aes(x = `Number of genes in set`, y = (`Number of DEGs`/`Number of genes in set`)*100)) +
                  
                  geom_point(aes(color = `Number of DEGs`, shape = database), size = 4) +
                  # geom_text(aes(x = `Number of genes in set` + 2, y = (`Number of DEGs`/`Number of genes in set`)*100 + 5,
                  #               label = Pathway), hjust = 0) +
                  geom_text_repel(aes(label = str_wrap(Pathway, width = 35)), color = "black", size = 4, box.padding = 0.15) +
                  facet_wrap(vars(database), ncol = 1, scales = "free") +
                  scale_color_gradient(low = "#4575B4", high = "#D73027") +
                  xlab("Number of genes in pathway") +
                  ylab("Percent of DEGs in pathway") +
                  ggtitle(names(object)) +
                  guides( shape = FALSE) +
                  #ylim(0,105) +
                  theme_bw() +
                  theme(legend.position = "top",
                        legend.justification = "right",
                        legend.direction = "horizontal",
                        axis.title = element_text(size = 20),
                        axis.text = element_text(size = 17),
                        plot.title = element_text(size = 23))
            
      }

      
      

      
      