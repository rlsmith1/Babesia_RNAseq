



# libraries ---------------------------------------------------------------

      library(tidyverse)
      library(readxl)
      library(stringi)


# read data --------------------------------------------------------------------


      read_excel_allsheets <- function(filename, tibble = FALSE) {
            
            sheets <- readxl::excel_sheets(filename)
            x <- lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
            if(!tibble) x <- lapply(x, as.data.frame)
            names(x) <- sheets
            x
            
      }
      
      l_clinicaldata <- read_excel_allsheets(filename = "clinical data/Case dataFINAL.xlsx", tibble = TRUE)
      
      l_clinicaldata_2 <- read_excel_allsheets(filename = "clinical data/Case dataFINAL_2.xlsx", tibble = TRUE)




# split data -------------------------------------------------------------

      df_habitus <- l_clinicaldata[[2]] %>% select(c(`Dog number`, `Exp group`, contains("H")))
      df_resp_rate <- l_clinicaldata[[2]] %>% select(c(`Dog number`, `Exp group`, contains("RR")))
      df_sys_bp <- l_clinicaldata[[2]] %>% select(c(`Dog number`, `Exp group`, contains("BP(SYS")))
      df_dias_bp <- l_clinicaldata[[2]] %>% select(c(`Dog number`, `Exp group`, contains("BP(DIA")))
      df_temp <- l_clinicaldata[[2]] %>% select(c(`Dog number`, `Exp group`, starts_with("T")))
      df_pulse <- l_clinicaldata[[2]] %>% select(c(`Dog number`, `Exp group`, starts_with("P")))
      df_appetite <- l_clinicaldata[[2]] %>% select(c(`Dog number`, `Exp group`, starts_with("A")))
      df_Hb <- l_clinicaldata[[3]] %>% select(c(`Dog number`, `Exp group`, starts_with("Hb")))
      df_par <- l_clinicaldata[[9]] 
      df_oxy_sat <- l_clinicaldata_2[[5]] %>% select(c(`Dog number`, `Exp group`, starts_with("OS"))) %>% .[1:6,]
      df_retic <- l_clinicaldata_2[[3]] %>% select(c(179:190)) %>% .[1:5,]
      df_white_cell <- l_clinicaldata[[3]] %>% select(c(`Dog number`, `Exp group`, starts_with("WCC")))
      df_neutrophil <- l_clinicaldata[[3]] %>% select(c(`Dog number`, `Exp group`, starts_with("SNeut")))
      df_lymophocyte <- l_clinicaldata[[3]] %>% select(c(`Dog number`, `Exp group`, starts_with("Lymph")))
      df_monocyte <- l_clinicaldata[[3]] %>% select(c(`Dog number`, `Exp group`, starts_with("Mono")))
      df_band_cell <- l_clinicaldata[[3]] %>% select(c(`Dog number`, `Exp group`, starts_with("BNeut")))
      df_platelet <- l_clinicaldata[[3]] %>% select(c(`Dog number`, `Exp group`, starts_with("Plt")))
      df_cfHb <- l_clinicaldata[[10]] 
      df_urea <- l_clinicaldata_2[[4]] %>% select(c(70:80)) %>% .[1:5,]
      df_creatinine <- l_clinicaldata[[4]] %>% select(c(`Dog number`, `Exp group`, starts_with("CREA")))
      df_lactate <- l_clinicaldata[[2]] %>% select(c(`Dog number`, `Exp group`, starts_with("LAC")))
      df_hco3 <- l_clinicaldata[[5]] %>% select(c(`Dog number`, `Exp group`, starts_with("HCO3")))
      df_pco2 <- l_clinicaldata[[5]] %>% select(c(`Dog number`, `Exp group`, starts_with("CO2")))
      df_po2 <- l_clinicaldata[[5]] %>% select(c(`Dog number`, `Exp group`, starts_with("O2")))
      df_pH <- l_clinicaldata[[5]] %>% select(c(`Dog number`, `Exp group`, starts_with("ph")))
      df_glucose <- l_clinicaldata[[2]] %>% select(c(`Dog number`, `Exp group`, starts_with("GLC")))
      df_crp <- l_clinicaldata[[4]] %>% select(c(`Dog number`, `Exp group`, starts_with("CRP")))
      df_alt <- l_clinicaldata[[4]] %>% select(c(`Dog number`, `Exp group`, starts_with("ALT")))
      df_cortisol <- l_clinicaldata_2[[6]] %>% select(c(1:12)) %>% .[1:5,]
      df_t4 <- l_clinicaldata_2[[6]] %>% select(c(15:24)) %>% .[1:5,]
      
      df_gm_csf <- l_clinicaldata_2[[8]] %>% select(c(`Dog number`, `Exp group`, 4:13)) %>% .[1:5,]
      df_keratinocyte <- l_clinicaldata_2[[8]] %>% select(c(`Dog number`, `Exp group`, 30:39)) %>% .[1:5,]
      df_il6 <- l_clinicaldata_2[[8]] %>% select(c(`Dog number`, `Exp group`, 65:74)) %>% .[1:5,]
      df_il8 <- l_clinicaldata_2[[8]] %>% select(c(`Dog number`, `Exp group`, 91:100)) %>% .[1:5,]
      df_il10 <- l_clinicaldata_2[[8]] %>% select(c(`Dog number`, `Exp group`, 104:113)) %>% .[1:5,]
      df_macrophage <- l_clinicaldata_2[[8]] %>% select(c(`Dog number`, `Exp group`, 139:149)) %>% .[1:5,]
      df_tnfa <- l_clinicaldata_2[[8]] %>% select(c(`Dog number`, `Exp group`, 154:163)) %>% .[1:5,]
      
   
save(df_habitus, df_resp_rate, df_sys_bp, df_dias_bp, df_temp, df_pulse, df_appetite, df_Hb, df_white_cell, df_band_cell,
     df_platelet, df_cfHb, df_creatinine, df_lactate, df_hco3, df_pco2, df_po2, df_pH, df_glucose, df_alt, df_par,
     df_retic, df_urea, df_oxy_sat, df_cortisol, df_t4, df_crp, df_neutrophil,
     df_gm_csf, df_keratinocyte, df_il6, df_il8, df_il10, df_macrophage, df_tnfa,
     file = "clinicalfigures.Rdata"
     )


# format data function -------------------------------------------------------------


      format_clin_data <- function(df){
         
         df1 <- df %>% 
            pivot_longer(cols = 3:ncol(.), names_to = "day", values_to = "value" ) 
         
         
         df1$day <- gsub("[^0-9.-]", "", df1$day) %>% as.factor()
         df1$`Exp group`[is.na(df1$`Exp group`)] <- 0
         df2 <- df1 %>% mutate(inoculum = ifelse(`Exp group` == 0, "low", "high"))
         
         df2
         
         
      }   
      

# plot function -------------------------------------------------------------------
      
      
      
      plot_clin_data <- function(df) {
         
         
         df %>% ggplot(aes(x = day, y = value, color = inoculum)) +
         geom_point(size = 3, alpha = 0.7) +
         stat_summary(aes(group = `Exp group`), fun = mean, geom = "line", size = 2) +
            xlim(c("-1","2", "3", "4", "5", "6", "7", "8")) +
         xlab("") +
         ylab("") +
         theme_classic() +
         theme(plot.title = element_text(hjust = 0.5))
      
      }
      

# plots -------------------------------------------------------------------

      
      # Vital signs
      
         # Temp
      
            df_temp %>% format_clin_data() %>% plot_clin_data()
      
         # Heart rate
      
            df_pulse %>% format_clin_data() %>% plot_clin_data()
      
         # Blood pressure (systolic/diastolic)
            
            df_sys_bp1 <- df_sys_bp %>% format_clin_data()
            df_dias_bp1 <- df_dias_bp %>% format_clin_data()
            
            df_bp <- df_sys_bp1 %>% 
               inner_join(df_dias_bp1, by = c("day", "Dog number", "Exp group", "inoculum")) %>% 
               mutate(value = value.x/value.y)
            
            df_bp %>% plot_clin_data()
            
         # Respiratory rate
            
            df_resp_rate %>% format_clin_data %>% plot_clin_data()
      
         # Oxygen saturation
      
         # Habitus
      
            df_habitus %>% format_clin_data() %>% plot_clin_data()
      
      # Hematology
      
         # Parasitemia
            
            df_par %>% format_clin_data() %>% plot_clin_data()
      
         # Hemoglobin
            
            df_Hb %>% format_clin_data() %>% plot_clin_data()
      
         # Reticulocyte
      
         # White blood cell
            
            df_white_cell %>% format_clin_data() %>% plot_clin_data()
      
         # Band forms
            
            df_band_cell %>% format_clin_data() %>% plot_clin_data()
      
         # Platelets
            
            df_platelet %>% format_clin_data() %>% plot_clin_data()
      
      # Biochemistry
      
         # plasma hemoglobin
            
            df_cfHb %>% format_clin_data() %>% plot_clin_data()
      
         # serum urea
      
         # serum creatinine
            
            df_creatinine %>% format_clin_data() %>% plot_clin_data()
      
         # blood HCO3
            
            df_hco3_1 <- df_hco3 %>% .[,1:10] %>% pivot_longer(cols = 3:ncol(.), names_to = "day", values_to = "value" )
            
            df_hco3_1$day <- gsub("HCO3", "", df_hco3_1$day)
            df_hco3_1$day <- gsub("[^0-9.-]", "", df_hco3_1$day) %>% as.factor()
            df_hco3_1$`Exp group`[is.na(df_hco3_1$`Exp group`)] <- 0
            df_hco3_2 <- df_hco3_1 %>% mutate(inoculum = ifelse(`Exp group` == 0, "low", "high"))
            
            df_hco3_2 %>% plot_clin_data()
            
         # blood gluose
            
            df_glucose %>% format_clin_data() %>% plot_clin_data()
      
         # ALT/AST
            
            df_alt %>% format_clin_data() %>% plot_clin_data()
      
      # Acid/base and endocrine
      
         # pH
            
            df_pH %>% format_clin_data() %>% plot_clin_data()
      
         # pCO2
            
            df_pco2_1 <- df_pco2 %>% .[,1:10] %>% pivot_longer(cols = 3:ncol(.), names_to = "day", values_to = "value" )
            
            df_pco2_1$day <- gsub("CO2", "", df_pco2_1$day)
            df_pco2_1$day <- gsub("[^0-9.-]", "", df_pco2_1$day) %>% as.factor()
            df_pco2_1$`Exp group`[is.na(df_pco2_1$`Exp group`)] <- 0
            df_pco2_2 <- df_pco2_1 %>% mutate(inoculum = ifelse(`Exp group` == 0, "low", "high"))
            
            df_pco2_2 %>% plot_clin_data()
            
      
         # pO2
            
            df_po2_1 <- df_po2 %>% .[,1:10] %>% pivot_longer(cols = 3:ncol(.), names_to = "day", values_to = "value" )
            
            df_po2_1$day <- gsub("O2", "", df_po2_1$day)
            df_po2_1$day <- gsub("[^0-9.-]", "", df_po2_1$day) %>% as.factor()
            df_po2_1$`Exp group`[is.na(df_po2_1$`Exp group`)] <- 0
            df_po2_2 <- df_po2_1 %>% mutate(inoculum = ifelse(`Exp group` == 0, "low", "high"))
            
            df_po2_2 %>% plot_clin_data()
            
         # lactate
            
            df_lactate %>% format_clin_data() %>% plot_clin_data()
      
         # cortisol
         
         # T4
      

     
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      







