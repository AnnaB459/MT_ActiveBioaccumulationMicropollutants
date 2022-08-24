# mass balance
# R-script which takes the QuEChERS data and the aquatic data (both only sorption
# control) and makes a mass balance to see whether actually about 100% of the 
# spiked substances can be found again

# Anna Boesel, 25.7.2022

# load libraries which are used
library(readr) # to read csv-file
library(tidyverse) # contains data science packages
library(naniar) # used for handling missing values
library(fs) # used for handling files 
library(dplyr) # contains many data manipulation functions, useful for dataframes
library(tidyr) # even more functions for data analysis
#install.packages("ggplot2")          # Install & load ggplot2 package
library("ggplot2") # for more plotting possibilities

# store filepaths where csv-files are located
filepaths_quechers <- fs::dir_ls("Rdata") # can directly use output file of QuEChERS data which is located here
filepaths_aq <- fs::dir_ls("Rdata\\biotransformation_aq")

# constants
V_init <- 0.5 # L, initial volume in batch
tt_weights <- c(0.689,0.5185,0.488,0.448)
mean_tt_weights <- c(0.60375,0.468)
# tt_weights <- c(0.1,0.1,0.1,0.1) # this is just hypothetic data trying to find the error
# import quechers data
quechers_data <- read_csv(file = "Rdata\\bt_exp_output.csv")


for (i in seq_along(filepaths_aq)) { #c(1:5))
  
  # Importing data and making some basic manipulation ------------------------------------------------------------------------------------------
  # read data from csv-file and store it in dataframe
  # structure: rows are different samples, columns are associated values to samples
  aq_data <- read_csv(file = filepaths_aq[[i]], show_col_types = FALSE)
  
  # replace all missing values (N/F) with NA
  aq_data[aq_data == "N/F"] <- NA
  aq_data[aq_data == "N/A"] <- NA
  
  # convert the columns which are later needed for calculations to numeric. If not done, the calculations cannot be executed properly. All NA are replaced with 0.
  # this contains the columns: 5 (Area), 12 (Calculated Amount), 23 (theoretical amount), 7 (response ratio), 6 (ISTD response)
  # in addition, wherever the calculated amount, area etc. is smaller than 0, it is replaced by 0
  co_to_conv <- c(5,6,7,12,21,23)
  for (j in co_to_conv) {
    aq_data[,j] <- as.numeric(unlist(aq_data[,j]))
    aq_data[,j][aq_data[,j] <= 0] <- 0
    aq_data[,j][is.na(aq_data[,j])] <- as.double(0)
  }
  
 
  # get name of current compound
  current_compound <- aq_data$Compound[1]
  # get molar mass of the compound
  if (aq_data$Adduct[10] == "M+H") { mm <- aq_data$`m/z (Expected)`[10] - 1.007825} # molar mass of current compound in g/mol
  if (aq_data$Adduct[10] == "M-H") { mm <- aq_data$`m/z (Expected)`[10] + 1.007825} # molar mass of current compound in g/mol
  
  # get the data you need
  sc_aq_bf <- aq_data %>% filter(str_detect(Filename,'biof_ab_scq')) %>% select('Calculated Amt') # concentration in ng/L of compound in aqueos biofilm samples for timepoints 0, 48h and 96h
  sc_aq_as <- aq_data %>% filter(str_detect(Filename,'asb_ab_scq')) %>% select('Calculated Amt') # concentration in ng/L of compound in aqueos activated sludge samples for timepoints 0, 48h and 96h
  bt_aq_bf <- aq_data %>% filter(str_detect(Filename,'biof_ab_btq')) %>% select('Calculated Amt') # concentration in ng/L of compound in aqueos biofilm samples for timepoints 0, 48h and 96h
  bt_aq_as <- aq_data %>% filter(str_detect(Filename,'asb_ab_btq')) %>% select('Calculated Amt') # concentration in ng/L of compound in aqueos activated sludge samples for timepoints 0, 48h and 96h
  
  
  # now calculate to how many nmol this corresponds, taking into account the initial volume and converting it to nmol 
  nmolsc_aq_bf <- t(sc_aq_bf*V_init/mm)
  nmolsc_aq_as <- t(sc_aq_as*V_init/mm)
  nmolbt_aq_bf <- t(bt_aq_bf*V_init/mm)
  nmolbt_aq_as <- t(bt_aq_as*V_init/mm)
  colnames(nmolsc_aq_bf) <- c('t0','t48','t96')
  colnames(nmolsc_aq_as) <- c('t0','t48','t96')
  colnames(nmolbt_aq_bf) <- c('t0','t48','t96')
  colnames(nmolbt_aq_as) <- c('t0','t48','t96')
  
  # now find corresponding quechers-data
  sc_qq_bf <- quechers_data %>% filter(str_detect(type,paste('biof_scq', sub("_.*", "", current_compound)))) %>% select(`C0 [nmol/g biomass]`,`mean C48 [nmol/g biomass]`,`mean C96 [nmol/g biomass]`)
  sc_qq_as <- quechers_data %>% filter(str_detect(type,paste('asb_scq', sub("_.*", "", current_compound)))) %>% select(`C0 [nmol/g biomass]`,`mean C48 [nmol/g biomass]`,`mean C96 [nmol/g biomass]`)
  bt_qq_bf <- quechers_data %>% filter(str_detect(type,paste('biof_btq', sub("_.*", "", current_compound)))) %>% select(`C0 [nmol/g biomass]`,`mean C48 [nmol/g biomass]`,`mean C96 [nmol/g biomass]`)
  bt_qq_as <- quechers_data %>% filter(str_detect(type,paste('asb_btq', sub("_.*", "", current_compound)))) %>% select(`C0 [nmol/g biomass]`,`mean C48 [nmol/g biomass]`,`mean C96 [nmol/g biomass]`)
  
  # calculate how many nmol are thus in the solid phase by multiplying it with the total biomass
  nmolsc_qq_bf <- sc_qq_bf*mean_tt_weights[1]
  nmolsc_qq_as <- sc_qq_as*mean_tt_weights[2]
  nmolbt_qq_bf <- bt_qq_bf*mean_tt_weights[1]
  nmolbt_qq_as <- bt_qq_as*mean_tt_weights[2]
  
  # substract the amount which was initially present (before spiking)
  nmolsc_qq_bf[1,1:3] <- nmolsc_qq_bf[1,1:3] - nmolsc_qq_bf[1,1]
  nmolsc_qq_as[1,1:3] <- nmolsc_qq_as[1,1:3] - nmolsc_qq_as[1,1]
  colnames(nmolsc_qq_bf) <- c('t0','t48','t96')
  colnames(nmolsc_qq_as) <- c('t0','t48','t96')
  nmolbt_qq_bf[1,1:3] <- nmolbt_qq_bf[1,1:3] - nmolbt_qq_bf[1,1]
  nmolbt_qq_as[1,1:3] <- nmolbt_qq_as[1,1:3] - nmolbt_qq_as[1,1]
  colnames(nmolbt_qq_bf) <- c('t0','t48','t96')
  colnames(nmolbt_qq_as) <- c('t0','t48','t96')
  
  # combine liquid and solid measurements
  liq_sol_bf <- as.matrix(rbind(nmolsc_aq_bf,nmolsc_qq_bf))
  liq_sol_as <- as.matrix(rbind(nmolsc_aq_as,nmolsc_qq_as))
  
  liq_sol_bf_bt <- as.matrix(rbind(nmolbt_aq_bf,nmolbt_qq_bf))
  liq_sol_as_bt <- as.matrix(rbind(nmolbt_aq_as,nmolbt_qq_as))
  
  # calculate the total nmol (which should be 5nmol ideally)
  nmol_total_bf <-nmolsc_qq_bf + nmolsc_aq_bf
  nmol_total_as <- nmolsc_qq_as + nmolsc_aq_as
  
  nmol_total_bf_bt <-nmolbt_qq_bf + nmolbt_aq_bf
  nmol_total_as_bt <- nmolbt_qq_as + nmolbt_aq_as
  
  # sorption control plotting
  data <- data.frame(time=c('t0','t48','t96','t0','t48','t96','t0','t48','t96','t0','t48','t96'),biocommunity=c('BF','BF','BF','AS','AS','AS','BF','BF','BF','AS','AS','AS'),phase=c('aqueous','aqueous','aqueous','aqueous','aqueous','aqueous','adsorbed','adsorbed','adsorbed','adsorbed','adsorbed','adsorbed'),amount=c(unlist(nmolsc_aq_bf),unlist(nmolsc_aq_as),unlist(as.matrix(nmolsc_qq_bf)),unlist(as.matrix(nmolsc_qq_as))))
  
  print(ggplot(data,                         # Draw barplot with grouping & stacking col=c('aliceblue','green3')
         aes(x = biocommunity,
             y = amount,
             fill = phase)) +
    geom_bar(stat = "identity",
             position = "stack") +
    facet_grid(~ time) +
    scale_fill_manual(values=c('grey','cadetblue2'))+ 
    geom_abline(slope=0,intercept=5) +
    theme_bw(18) +
    theme(legend.position = 'right') +
    theme(legend.key.size = unit(1.5, 'cm')) +
    labs(title = (paste(sub("_.*", "", current_compound), 'mass balance sorption control')),y = "amount [nmol]"))
  
  print(paste(current_compound, 'biofilm',nmolsc_aq_bf,t(nmolsc_qq_bf),'activated sludge',nmolsc_aq_as,t(nmolsc_qq_as)))
  
  
  # for biotic experiment
  data <- data.frame(time=c('t0','t48','t96','t0','t48','t96','t0','t48','t96','t0','t48','t96'),biocommunity=c('BF','BF','BF','AS','AS','AS','BF','BF','BF','AS','AS','AS'),phase=c('aqueous','aqueous','aqueous','aqueous','aqueous','aqueous','adsorbed','adsorbed','adsorbed','adsorbed','adsorbed','adsorbed'),amount=c(unlist(nmolbt_aq_bf),unlist(nmolbt_aq_as),unlist(as.matrix(nmolbt_qq_bf)),unlist(as.matrix(nmolbt_qq_as))))
  
  print(ggplot(data,                         # Draw barplot with grouping & stacking col=c('aliceblue','green3')
               aes(x = biocommunity,
                   y = amount,
                   fill = phase)) + 
          geom_bar(stat = "identity",
                   position = "stack") +
          facet_grid(~ time) +
          scale_fill_manual(values=c('grey','cadetblue2'))+ 
          geom_abline(slope=0,intercept=5) +
          theme_bw(18) +
          theme(legend.position = 'right') +
          theme(legend.key.size = unit(1.5, 'cm')) +
          labs(title = (paste(sub("_.*", "", current_compound), 'mass balance biotic experiment')),y = "amount [nmol]"))
  
  print(paste(current_compound, 'biofilm',nmolbt_aq_bf,t(nmolbt_qq_bf),'activated sludge',nmolbt_aq_as,t(nmolbt_qq_as)))
  
}