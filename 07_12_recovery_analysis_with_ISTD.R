# File to calculate recoveries for substances with an ISTD
# Anna Boesel, 06/17/2022

# load libraries which are used
library(readr) # to read csv-file
library(tidyverse) # contains data science packages
library(naniar) # used for handling missing values
library(fs) # used for handling files 
library(dplyr) # contains many data manipulation functions, useful for dataframes
library(tidyr) # even more functions for data analysis


# constants
conv_fac <- 1.2/1.75 #conversion factor for schemes 3 and 4 (there was a mistake in the calculation), can be used to correct the concentration
theo_conc= 4 #theoretical concentration = 4 nmol/L
AS_3_5_conv_fac = 4.2/7.7 #conversion factor needed due to a manipulation error in sample activated sludge, scheme 3, replicate 5
IS_conv <- 4/5.25
# store filepaths where csv-files are located
# all ISTD3-files are located in one folder, all ISTD4-files are located in another folder
filepaths_ISTD3 <- fs::dir_ls("Rdata\\recovery_ISTD3_data") 
filepaths_ISTD4 <- fs::dir_ls("Rdata\\recovery_ISTD4_data")

# initialise variables
ISTD3_data <- list()
ISTD4_data <- list()
finaloutput <- list()
final_cc_output <- list()


# start a for loop which runs through all the csv-files located in the given filepaths, one after the other and analyses it
for (i in seq_along(filepaths_ISTD3)) {
  
  # Importing data and making some basic manipulation ------------------------------------------------------------------------------------------
  # read data from csv-file and store it in dataframe called ISTDx_data
  # structure: rows are different samples, columns are associated values to samples
  ISTD3_data <- read_csv(file = filepaths_ISTD3[[i]], show_col_types = FALSE)
  ISTD4_data <- read_csv(file = filepaths_ISTD4[[i]], show_col_types = FALSE)
  
  # replace all missing values (N/F) with NA
  ISTD3_data[ISTD3_data == "N/F"] <- NA
  ISTD4_data[ISTD4_data == "N/F"] <- NA # replacing all N/F entries with NA
  
  # convert the columns which are later needed for calculations to numeric. If not done, the calculations cannot be executed properly. All NA are replaced with 0.
  # this contains the columns: 5 (Area), 12 to (Calculated Amount), 23 (theoretical amount), 7 (response ratio), 6 (ISTD response)
  # in addition, wherever the calculated amount, area etc. is smaller than 0, it is replaced by 0
  co_to_conv <- c(5,6,7,12,21,23)
  for (j in co_to_conv) {
    ISTD3_data[,j] <- as.numeric(unlist(ISTD3_data[,j]))
    ISTD3_data[,j][ISTD3_data[,j] <= 0] <- 0
    ISTD3_data[,j][is.na(ISTD3_data[,j])] <- as.double(0)
    ISTD4_data[,j] <- as.numeric(unlist(ISTD4_data[,j]))
    ISTD4_data[,j][is.na(ISTD4_data[,j])] <- as.double(0)
  }
  
  # get the name and the molar mass [g/mol] of the current compound
  current_compound <- ISTD3_data$Compound[1]
  if (ISTD3_data$Adduct[1] == "M+H") { mm <- ISTD3_data$`m/z (Expected)`[1] - 1.007825} # molar mass of current compound in g/mol
  if (ISTD3_data$Adduct[1] == "M-H") { mm <- ISTD3_data$`m/z (Expected)`[1] + 1.007825} # molar mass of current compound in g/mol
  
  #check if correct two csv were chosen, if it is really the same compound in both csv
  if (current_compound != ISTD4_data$Compound[1]) {
    print(cat('there is a problem, files are not corresponding. One file is for ', current_compound, 'the other is for', ISTD4_data$Compound[1]))
    the_problem_is_present = TRUE
  }
  
  
  # Take out all specific information needed from data and make dataframes -----------------------------------------------------------------------  
  
  # make a dataframe containing only the concentrations of all the samples
  # correct the concentrations in scheme 3 and scheme 4 as there was too much spiked
  # each column corresponds to one scheme, each row corresponds to the same standard addition spiking amount
  AS_conc <- cbind(ISTD3_data %>% filter(str_detect(Filename, 'AS_S1')) %>% select('Calculated Amt'), 
                   ISTD3_data %>% filter(str_detect(Filename, 'AS_S2')) %>% select('Calculated Amt')*conv_fac,
                   ISTD4_data %>% filter(str_detect(Filename, 'AS_S3')) %>% select('Calculated Amt')*conv_fac,
                   ISTD4_data %>% filter(str_detect(Filename, 'AS_S4')) %>% select('Calculated Amt')*conv_fac,
                   #ISTD4_data %>% filter(str_detect(Filename, 'AS_S3')) %>% select('Calculated Amt'),
                   #ISTD4_data %>% filter(str_detect(Filename, 'AS_S4')) %>% select('Calculated Amt'),
                   ISTD3_data %>% filter(str_detect(Filename, 'AS_S5')) %>% select('Calculated Amt'))
  #correct concentration in scheme 3, replicate 5
  AS_conc[5,3] <- AS_conc[5,3]*AS_3_5_conv_fac
  # name the columns
  colnames(AS_conc) <- c("scheme 1", "scheme 2", "scheme 3", "scheme 4", "scheme 5")
  
  # the same for biofilm
  BF_conc <- cbind(ISTD3_data %>% filter(str_detect(Filename, 'BF_S1')) %>% select('Calculated Amt'), 
                   ISTD3_data %>% filter(str_detect(Filename, 'BF_S2')) %>% select('Calculated Amt')*conv_fac,
                   ISTD4_data %>% filter(str_detect(Filename, 'BF_S3')) %>% select('Calculated Amt')*conv_fac,
                   ISTD4_data %>% filter(str_detect(Filename, 'BF_S4')) %>% select('Calculated Amt')*conv_fac,
                   #ISTD4_data %>% filter(str_detect(Filename, 'BF_S3')) %>% select('Calculated Amt'),
                   #ISTD4_data %>% filter(str_detect(Filename, 'BF_S4')) %>% select('Calculated Amt'),
                   ISTD3_data %>% filter(str_detect(Filename, 'BF_S5')) %>% select('Calculated Amt'))
  # name the columns
  colnames(BF_conc) <- c("scheme 1", "scheme 2", "scheme 3", "scheme 4", "scheme 5")
  
  
  # do the same for the areas
  # make a dataframe containing only the areas of all the samples
  # each column corresponds to one scheme, each row corresponds to the same standard addition spiking amount
  # ATTENTION: areas are not corrected yet
  AS_area <- cbind(ISTD3_data %>% filter(str_detect(Filename, 'AS_S1')) %>% select('Area'), 
                   ISTD3_data %>% filter(str_detect(Filename, 'AS_S2')) %>% select('Area'),
                   ISTD4_data %>% filter(str_detect(Filename, 'AS_S3')) %>% select('Area'),
                   ISTD4_data %>% filter(str_detect(Filename, 'AS_S4')) %>% select('Area'),
                   ISTD3_data %>% filter(str_detect(Filename, 'AS_S5')) %>% select('Area'))
  # name the columns
  colnames(AS_area) <- c("scheme 1", "scheme 2", "scheme 3", "scheme 4", "scheme 5")
  
  # the same for biofilm
  BF_area <- cbind(ISTD3_data %>% filter(str_detect(Filename, 'BF_S1')) %>% select('Area'), 
                   ISTD3_data %>% filter(str_detect(Filename, 'BF_S2')) %>% select('Area'),
                   ISTD4_data %>% filter(str_detect(Filename, 'BF_S3')) %>% select('Area'),
                   ISTD4_data %>% filter(str_detect(Filename, 'BF_S4')) %>% select('Area'),
                   ISTD3_data %>% filter(str_detect(Filename, 'BF_S5')) %>% select('Area'))
  # name the columns
  colnames(BF_area) <- c("scheme 1", "scheme 2", "scheme 3", "scheme 4", "scheme 5")
  
  
  # do the same for ISTD-response
  # make a dataframe containing only the ISTD Responses of all the samples
  # each column corresponds to one scheme, each row corresponds to the same standard addition spiking amount
  AS_ISTD_Response <- cbind(ISTD3_data %>% filter(str_detect(Filename, 'AS_S1')) %>% select('ISTD Response'), 
                            ISTD3_data %>% filter(str_detect(Filename, 'AS_S2')) %>% select('ISTD Response'),
                            ISTD4_data %>% filter(str_detect(Filename, 'AS_S3')) %>% select('ISTD Response'),
                            ISTD4_data %>% filter(str_detect(Filename, 'AS_S4')) %>% select('ISTD Response'),
                            ISTD3_data %>% filter(str_detect(Filename, 'AS_S5')) %>% select('ISTD Response'))
  
  # name the columns
  colnames(AS_ISTD_Response) <- c("scheme 1", "scheme 2", "scheme 3", "scheme 4", "scheme 5")
  
  # the same for biofilm
  BF_ISTD_Response <- cbind(ISTD3_data %>% filter(str_detect(Filename, 'BF_S1')) %>% select('ISTD Response'), 
                            ISTD3_data %>% filter(str_detect(Filename, 'BF_S2')) %>% select('ISTD Response'),
                            ISTD4_data %>% filter(str_detect(Filename, 'BF_S3')) %>% select('ISTD Response'),
                            ISTD4_data %>% filter(str_detect(Filename, 'BF_S4')) %>% select('ISTD Response'),
                            ISTD3_data %>% filter(str_detect(Filename, 'BF_S5')) %>% select('ISTD Response'))
  # name the columns
  colnames(BF_ISTD_Response) <- c("scheme 1", "scheme 2", "scheme 3", "scheme 4", "scheme 5")
  
  # get areas of ISTD in scheme 1 (needed for matrix factor calculation)
  AS_S1_ISTD_area <- (ISTD3_data %>% filter(str_detect(Filename, 'AS_S1')) %>% select('ISTD Response'))
  BF_S1_ISTD_area <- (ISTD3_data %>% filter(str_detect(Filename, 'BF_S1')) %>% select('ISTD Response'))

  # get the calibration data and make a linear regression 
  ISTD3_cal <- ISTD3_data %>% filter(str_detect(Filename, 'Cal'))
  ISTD3_cal <- ISTD3_cal[ISTD3_cal$`Response Ratio` !=0,] # remove the entries where nothing is detected
  ISTD4_cal <- ISTD4_data %>% filter(str_detect(Filename, 'Cal'))
  ISTD4_cal <- ISTD4_cal[ISTD4_cal$`Response Ratio` !=0,]
  
  # linear regression
  lin_regr3 <- coefficients(lm(unlist(ISTD3_cal%>%select('Response Ratio'))~unlist(ISTD3_cal%>%select('Theoretical Amt')), weights = 1/unlist(ISTD3_cal%>%select('Theoretical Amt'))))
  lin_regr4  <- coefficients(lm(unlist(ISTD4_cal%>%select('Response Ratio'))~unlist(ISTD4_cal%>%select('Theoretical Amt')), weights = 1/unlist(ISTD4_cal%>%select('Theoretical Amt'))))
  m <- data.frame((matrix(ncol = 5, nrow = 6))) # reserve space for a dataframe containing all m (slopes) of the regression. They are different depending on the scheme (as there is a different calibration curve / there was an error)
  b <- data.frame((matrix(ncol = 5, nrow = 6))) # reserve space for a dataframe containing all b (intersections with y-axis) of the regression
  for (i in 1:6){m[i,1:5] <- c(lin_regr3[2],lin_regr3[2],lin_regr4[2],lin_regr4[2],lin_regr3[2])} # slopes
  for (i in 1:6){b[i,1:5] <- c(lin_regr3[1],lin_regr3[1],lin_regr4[1],lin_regr4[1],lin_regr3[1])} # intersections
  
  # calculate the corrected, "true" area
  true_area_AS <- AS_ISTD_Response*m*AS_conc+AS_ISTD_Response*b
  true_area_BF <- BF_ISTD_Response*m*BF_conc+BF_ISTD_Response*b
  
  # calculate the area at 4nM in the calibration curve
  ISTD3_area <- mean(unlist(ISTD3_cal %>% select('ISTD Response')), na.rm=TRUE)
  ISTD4_area <- mean(unlist(ISTD4_cal %>% select('ISTD Response')), na.rm=TRUE)
  
  cali_level2 <- ISTD4_data %>% filter(str_detect(Filename, 'Cali')) %>% select('Level')#*0.97/0.96
  cali_area2 <- ISTD4_data %>% filter(str_detect(Filename, 'Cali')) %>% select('Area')
  
  area_interp2 <- approx(unlist(cali_level2), unlist(cali_area2),theo_conc)
  area_at_4nM4 <- area_interp2[2]
  
  # calculate recoveries -----------------------------------------------------------------------------------------------------------------------
  AS_abs_recovery <- (true_area_AS[,3]-true_area_AS[,4])/unlist(area_at_4nM4)#[3]) # 3rd column is scheme 3 minus 4rd column (schmeme 4), divide by area in calibration
  AS_rel_recovery <- (AS_conc[,3]-AS_conc[,4])/mm/(theo_conc*IS_conv)
  AS_extr_efficiency <- (AS_conc[,2]-AS_conc[,4]/IS_conv)/(AS_conc[,1]-AS_conc[,4]/IS_conv)
  AS_matrix_factor <- AS_S1_ISTD_area/ISTD3_area
  
  BF_abs_recovery <-(true_area_BF[,3]-true_area_BF[,4])/unlist(area_at_4nM4)#[3])
  BF_rel_recovery <- (BF_conc[,3]-BF_conc[,4])/mm/(theo_conc*IS_conv)
  BF_extr_efficiency <- (BF_conc[,2]-BF_conc[,4]/IS_conv)/(BF_conc[,1]-BF_conc[,4]/IS_conv)
  BF_matrix_factor <- BF_S1_ISTD_area/ISTD3_area
  
  # store the output ----------------------------------------------------------------------------------------------------------------------------
  the_output <- tibble()
  for (k in 1:6) {
    the_output <- rbind(the_output, tibble('Compound name' = current_compound, 'AS_abs_recovery' = AS_abs_recovery[k], 'AS_rel_recovery' = AS_rel_recovery[k], 'AS_extr_efficiency' = AS_extr_efficiency[k], 'AS_matrix_factor' = AS_matrix_factor[k,1], 'BF_abs_recovery' = BF_abs_recovery[k], 'BF_rel_recovery' = BF_rel_recovery[k], 'BF_extr_efficiency' = BF_extr_efficiency[k], 'BF_matrix_factor' = BF_matrix_factor[k,1]))
   }
  the_output <- rbind(the_output, tibble('Compound name' = 'mean', 'AS_abs_recovery' = mean(unlist(AS_abs_recovery[c(1:4,6)])), 'AS_rel_recovery' = mean(unlist(AS_rel_recovery[c(1:4,6)])), 'AS_extr_efficiency' = mean(unlist(AS_extr_efficiency[c(1:4,6)])), 'AS_matrix_factor' = mean(unlist(AS_matrix_factor[c(1:4,6),1])), 'BF_abs_recovery' = mean(unlist(BF_abs_recovery)), 'BF_rel_recovery' = mean(unlist(BF_rel_recovery)), 'BF_extr_efficiency' = mean(unlist(BF_extr_efficiency)), 'BF_matrix_factor' = mean(unlist(BF_matrix_factor))))
  the_output <- rbind(the_output, tibble('Compound name' = 'standard deviation', 'AS_abs_recovery' = sd(unlist(AS_abs_recovery[c(1:4,6)])), 'AS_rel_recovery' = sd(unlist(AS_rel_recovery[c(1:4,6)])), 'AS_extr_efficiency' = sd(unlist(AS_extr_efficiency[c(1:4,6)])), 'AS_matrix_factor' = sd(unlist(AS_matrix_factor[c(1:4,6),1])), 'BF_abs_recovery' = sd(unlist(BF_abs_recovery)), 'BF_rel_recovery' = sd(unlist(BF_rel_recovery)), 'BF_extr_efficiency' = sd(unlist(BF_extr_efficiency)), 'BF_matrix_factor' = sd(unlist(BF_matrix_factor[c(1:4,6),1]))))
  
  the_output[ nrow(the_output) + 1 , ] <- NA
  finaloutput <- rbind(finaloutput,the_output)
  
  cc_output <- paste(sub("_.*", "", current_compound), 'activated sludge,', format(round(mean(unlist(AS_abs_recovery[c(1:4,6)])),3),nsmall=3), '+/-', format(round(sd(unlist(AS_abs_recovery[c(1:4,6)])),3),nsmall=3), ',', format(round(mean(unlist(AS_rel_recovery[c(1:4,6)])),3),nsmall=3), '+/-', format(round(sd(unlist(AS_rel_recovery[c(1:4,6)])),3),nsmall=3), ',', format(round(mean(unlist(AS_extr_efficiency[1])),3),nsmall=3), ',', format(round(mean(unlist(AS_matrix_factor[c(1:4,6),1])),3),nsmall=3), '+/-', format(round(sd(unlist(AS_matrix_factor[c(1:4,6),1])),3),nsmall=3))
  cc_output <- rbind(cc_output, paste(sub("_.*", "", current_compound), 'biofilm,', format(round(mean(unlist(BF_abs_recovery[c(1:4,6)])),3),nsmall=3), '+/-', format(round(sd(unlist(BF_abs_recovery[c(1:4,6)])),3),nsmall=3), ',', format(round(mean(unlist(BF_rel_recovery[c(1:4,6)])),3),nsmall=3), '+/-', format(round(sd(unlist(BF_rel_recovery[c(1:4,6)])),3),nsmall=3), ',', format(round(mean(unlist(BF_extr_efficiency[1])),3),nsmall=3), ',', format(round(mean(unlist(BF_matrix_factor[c(1:4,6),1])),3),nsmall=3), '+/-', format(round(sd(unlist(BF_matrix_factor[c(1:4,6),1])),3),nsmall=3)))
  
  final_cc_output <- rbind(final_cc_output,cc_output)
}

# write output to csv-file-----------------------------------------------------------------------------------------------------------------------
write.csv(finaloutput,"Rdata\\output_recovery_ISTD.csv", row.names = FALSE)
write.csv(final_cc_output,"Rdata\\ISTD_output_condensed.csv",row.names = FALSE)