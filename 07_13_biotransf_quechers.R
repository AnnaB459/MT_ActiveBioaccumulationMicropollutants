# File for data evaluation of the biotransformation experiment with quechers
# Anna Boesel
# created on 8th of July, 2022

# load libraries which are used
library(readr) # to read csv-file
library(tidyverse) # contains data science packages
library(naniar) # used for handling missing values
library(fs) # used for handling files 
library(dplyr) # contains many data manipulation functions, useful for dataframes
library(tidyr) # even more functions for data analysis
library(ggplot2) # for plotting
library(ggpattern) # for plotting patterns

# store filepaths where csv-files are located
filepaths_quechers <- fs::dir_ls("Rdata\\biotransformation_quechers")

# import weights of biomasses
bweights <- read_csv(file = "Rdata\\Biomass_weights.csv")
bweights[,2] <- as.numeric(unlist(bweights[,2]))

# location of all the functions which will be needed in this code
source("functions_for_quechers.R") # where the functions are located

# initialise some variables
paecklis <- c("biof_btq","biof_scq", "asb_btq", "asb_scq")
output <- list()
all_bioacc <- list()

# run through all substances
for (i in seq_along(filepaths_quechers)) {
  
  # Importing data and making some basic manipulation --------------------------
  # read data from csv-file and store it in dataframe called ISTDx_data
  # structure: rows are different samples, columns are associated values to samples
  quechers_data <- read_csv(file = filepaths_quechers[[i]], show_col_types = FALSE)
  
  # replace all missing values (N/F) with NA
  quechers_data[quechers_data == "N/F"] <- NA
  quechers_data[quechers_data == "N/A"] <- NA
  
  # convert the columns which are later needed for calculations to numeric. If not done, the calculations cannot be executed properly. All NA are replaced with 0.
  # this contains the columns: 5 (Area), 12 (Calculated Amount), 23 (theoretical amount), 7 (response ratio), 6 (ISTD response)
  # in addition, wherever the calculated amount, area etc. is smaller than 0, it is replaced by 0
  co_to_conv <- c(5,6,7,12,21,23)
  for (j in co_to_conv) {
    quechers_data[,j] <- as.numeric(unlist(quechers_data[,j]))
    quechers_data[,j][quechers_data[,j] <= 0] <- 0
    quechers_data[,j][is.na(quechers_data[,j])] <- as.double(0)
  }
  
  # add the column with the weights
  quechers_data$bweights <- bweights$Weight
  
  # remove rows where not sufficient biomass was spiked
  quechers_data <- quechers_data[-c(27,39,98),]
  
  # get name of current compound
  current_compound <- quechers_data$Compound[1]
  # get molar mass of the compound
  if (quechers_data$Adduct[10] == "M+H") { mm <- quechers_data$`m/z (Expected)`[10] - 1.007825} # molar mass of current compound in g/mol
  if (quechers_data$Adduct[10] == "M-H") { mm <- quechers_data$`m/z (Expected)`[10] + 1.007825} # molar mass of current compound in g/mol
  
  # check whether it is a compound which has an IS or not
  if (quechers_data$`ISTD Response`[1]==0) {IS=FALSE} else {IS=TRUE}
  
  # calculations for substances with IS ----------------------------------------
  # if it has an IS, call function which calculates the mean concentration and standard deviations
  if (IS) {
    # select only the needed compound data and store it in a dataframe
    compound_data <- as.data.frame(cbind(quechers_data$Filename,quechers_data$`Calculated Amt`,quechers_data$bweights))
    colnames(compound_data) <- c("Filename","Concentration in ng/L","biomass weight in g")
    compound_data[,2:3] <- as.numeric(unlist(compound_data[,2:3])) # store as numeric in order to be able to perform calculations
    # initialize list
    final_conc <- list()
    
    # run for-loop in which for one "bottle experiment" after the other. The concentration changes in ng substance / g biomass are calculated
    for (j in 1:length(paecklis)) {
      # choose current experiment type
      compound_data_inside <- compound_data %>% filter(str_detect(Filename,paecklis[j]))
      
      # call the function which calculates the concentrations for all time-points and their respective standard deviation which is called IS_concentrations
      concentrations <- as.data.frame(IS_concentrations(compound_data_inside)) # C0, C48, C48 stdv, C96, C96 stdv. In ng substance / g biomass
      
      # do some formatting
      conc_for_which <- cbind(paste(paecklis[j],sub("_.*", "", current_compound)),concentrations,0)
      colnames(conc_for_which) <- c('type', 'C0 [nmol/g biomass]','mean C48 [nmol/g biomass]','stdv C48 [nmol/g biomass]', 'mean C96 [nmol/g biomass]','stdv C96 [nmol/g biomass]','stdv C0 if existent')
      final_conc <- rbind(final_conc,conc_for_which)
      
    }
    final_conc[1:4,2:6] <- final_conc[1:4,2:6]/mm # now in units nmol substance / g biomass
    
  }
  
  # calculations for substances with standard addition -------------------------
  else {
    # group data needed to calculate the concentration
    compound_data <- as.data.frame(cbind(quechers_data$Filename,quechers_data$`Area`,quechers_data$bweights))
    colnames(compound_data) <- c("Filename","Area","biomass weight in g")
    compound_data[,2:3] <- as.numeric(unlist(compound_data[,2:3])) # store as numeric in order to be able to perform calculations
    final_conc <- list() # initialize
    # run through all experiment types
    for (j in 1:length(paecklis)) {
      # choose current experiment type
      compound_data_inside <- compound_data %>% filter(str_detect(Filename,paecklis[j]))
      
      # now call the function which is stda_function which calculates the concentrations
      stda_output <- as.data.frame(stda_function(compound_data_inside)) # C0, C48, C48 stdv, C96, C96 stdv, C0 stdv, R2 of the linear regression. Already in nmol/g
      # do some formatting
      concentrations <- stda_output[1:6]
      conc_for_which <- cbind(paste(paecklis[j],sub("_.*", "", current_compound)),concentrations)
      colnames(conc_for_which) <- c('type', 'C0 [nmol/g biomass]','mean C48 [nmol/g biomass]','stdv C48 [nmol/g biomass]', 'mean C96 [nmol/g biomass]','stdv C96 [nmol/g biomass]','stdv C0 if existent')
      final_conc <- rbind(final_conc,conc_for_which)
    }
    
  }
  
  # calculate the bioaccumulation ----------------------------------------------
  #call bioaccumulation function
  bioacc_data <- bioaccumulation(final_conc) # gives back the bioaccumulation in nmol/g compared to timepoint 0h
  # do some formatting
  bio_names <- c('bioaccumulation_biof','bioaccumulation_asb')
  bioacc_data <- cbind(paste(bio_names,sub("_.*", "", current_compound)),bioacc_data)
  colnames(bioacc_data) <- c('type', 'C0 [nmol/g biomass]','mean C48 [nmol/g biomass]','stdv C48 [nmol/g biomass]', 'mean C96 [nmol/g biomass]','stdv C96 [nmol/g biomass]','stdv C0')
  
  # save in output-frame
  all_bioacc <- rbind(all_bioacc,bioacc_data)
  output <- rbind(output,final_conc[,1:6],bioacc_data[,1:6])
  
  # and plot with bar plots
  bar_plotting(final_conc,current_compound)
}

# write results to an csv-file
write.csv(output,"Rdata\\bt_exp_output.csv", row.names = FALSE)