# load libraries which are used
library(readr)
library(tidyverse)
library(naniar)
library(fs)
library(dplyr)
library(tidyr)
library("ggplot2")

options(digits = 4) # WARNING: changes number of digits in whole code, globally. To have more digit, comment this line and restart Rstudio

# constants
conv_fac <- 1.2/1.75 #conversion factor for schemes 3 and 4 (there was a mistake in the calculation)
AS_3_5_conv_fac = 4.2/7.7 #conversion factor needed due to a manipulation error
theo_conc <- 4 #theoretical concentration = 4 nmol/L
stda_conc <- c(0, 2, 2, 10, 20, 20) # nmol/L
stda_conc_AS <- c(0, 2, 2, 10, 20) # nmol/L
stda_conc_AS3 <- c(0, 2, 2, 10,40,20) # nmol/L
x <- seq(from=-5,to=25,length.out=6)
AS_conc <- list()
BF_conc <- list()
r2_AS <- list()
r2_BF <- list()
AS_conc_stdv <- list()
BF_conc_stdv <- list()
the_output <- tibble()
condensed_output <- tibble()


# import data from csv
filepaths_stda <- fs::dir_ls("Rdata\\recovery_stda_data")

stda_data <- list()
finaloutput <- list()
final_cc_output <- list()

for (i in seq_along(filepaths_stda)) {
  stda_data <- read_csv(file = filepaths_stda[[i]], show_col_types = TRUE)
  
  stda_data[stda_data == "N/F"] <- NA # replacing all N/F entries with NA
  
  # convert column 5 to numeric (Area)
  stda_data[,5] <- as.numeric(unlist(stda_data[,5]))
  stda_data[,5][stda_data[,5] <= 0] <- 0
  stda_data[,5][is.na(stda_data[,5])] <- as.double(0)
  
  #molar mass [g/mol]
  current_compound <- stda_data$Compound[1]
  mm <- stda_data$`m/z (Expected)`[1] # molar mass of current compound in g/mol
  if (stda_data$Adduct[1] == "M+H") { mm <- stda_data$`m/z (Expected)`[1] - 1.007825} # molar mass of current compound in g/mol
  if (stda_data$Adduct[1] == "M-H") { mm <- stda_data$`m/z (Expected)`[1] + 1.007825} # molar mass of current compound in g/mol
  
  
  # get calibration information
  cali_level <- stda_data %>% filter(str_detect(Filename, 'Cali')) %>% select('Level')#*0.97/0.96
  cali_area <- stda_data %>% filter(str_detect(Filename, 'Cali')) %>% select('Area')
  
  area_interp <- approx(unlist(cali_level), unlist(cali_area),theo_conc)
  area_at_4nM <- area_interp[2]
  
  
  # get areas of samples
  AS_area <- cbind(stda_data %>% filter(str_detect(Filename, 'AS_S1')) %>% select('Area'), 
                   stda_data %>% filter(str_detect(Filename, 'AS_S2')) %>% select('Area')*conv_fac,
                   stda_data %>% filter(str_detect(Filename, 'AS_S3')) %>% select('Area')*conv_fac,
                   stda_data %>% filter(str_detect(Filename, 'AS_S4')) %>% select('Area')*conv_fac,
                   stda_data %>% filter(str_detect(Filename, 'AS_S5')) %>% select('Area'))
  
  #correct concentration in scheme 3, replicate 5
  #AS_area[5,3] <- AS_area[5,3]*AS_3_5_conv_fac
  # I will remove this sample, no correction needed
  
  BF_area <- cbind(stda_data %>% filter(str_detect(Filename, 'BF_S1')) %>% select('Area'), 
                   stda_data %>% filter(str_detect(Filename, 'BF_S2')) %>% select('Area')*conv_fac,
                   stda_data %>% filter(str_detect(Filename, 'BF_S3')) %>% select('Area')*conv_fac,
                   stda_data %>% filter(str_detect(Filename, 'BF_S4')) %>% select('Area')*conv_fac,
                   stda_data %>% filter(str_detect(Filename, 'BF_S5')) %>% select('Area'))
  
  
  # make a least-square linear regression to find the concentration thanks to standard additions
  for (j in 1:dim(AS_area)[2]) {
    if (j==3) {   # if scheme 3, the regression will be done without sample 5 for the activated sludge samples
      AS_regr <- lm(unlist(AS_area[c(1:4,6),j]) ~stda_conc_AS)
      AS_i_coeffs <- coefficients(AS_regr) # the coefficients of the regression 
      AS_i_coeffs_stdv <- unlist(summary(AS_regr)[4])[3:4] # the standard deviation of the coefficients
      r2_AS[j] = format(summary(AS_regr)$r.squared,digits=3) # r2 to evaluate goodness of fit
      
      # make a plot with the found regression
      plot(stda_conc_AS,AS_area[c(1:4,6),j],
           main=paste(current_compound, 'in activated sludge. Scheme', j),
           xlab='concentration of standard addition [nM]', ylab='area', xlim=c(-5,22))
      grid()
      abline(AS_i_coeffs)
      text(15,AS_area[1,j],paste0("R^2 = ", r2_AS[j]),1)
      lines(x,(AS_i_coeffs[2]-AS_i_coeffs_stdv[2])*x+AS_i_coeffs_stdv[1]+AS_i_coeffs[1],col="red")
      lines(x,(AS_i_coeffs[2]+AS_i_coeffs_stdv[2])*x-AS_i_coeffs_stdv[1]+AS_i_coeffs[1],col="orange")
      legend("topleft",
             c("measured points","linear regression","left stdv","right stdv"),
             col=c("black","black","red","orange"), pch = c(1,15,15,15))
      
      # do a regression for biofilm
      BF_regr <- lm(unlist(BF_area[,j])~stda_conc)
      BF_i_coeffs <- coefficients(BF_regr)
      BF_i_coeffs_stdv <- unlist(summary(BF_regr)[4])[3:4]
      r2_BF[j] = format(summary(BF_regr)$r.squared,digits=3)
      
      # make a plot with the found regression
      plot(stda_conc,BF_area[,j],
           main=paste(current_compound, 'in biofilm. Scheme', j),
           xlab='concentration of standard addition [nM]', ylab='area', xlim=c(-5,22), ylim=c(0,max(BF_area[,j])+10^6))
      grid()
      abline(BF_i_coeffs)
      text(15,BF_area[1,j],paste0("R^2 = ", r2_BF[j]),1)
      lines(x,(BF_i_coeffs[2]-BF_i_coeffs_stdv[2])*x+BF_i_coeffs_stdv[1]+BF_i_coeffs[1],col="red")
      lines(x,(BF_i_coeffs[2]+BF_i_coeffs_stdv[2])*x-BF_i_coeffs_stdv[1]+BF_i_coeffs[1],col="orange")
      legend("topleft",
             c("measured points","linear regression","left stdv","right stdv"),
             col=c("black","black","red","orange"), pch = c(1,15,15,15)
      )
      
      this_data <- as.data.frame(t(rbind(stda_conc, BF_area[,j])))
      colnames(this_data) <- c('concentration','area')
      print(ggplot(this_data, aes(x=concentration, y=area)) +
              geom_point(color='#2980B9', size=4) +
              xlim(-10,22) +
              geom_smooth(method=lm, color='#2C3E50', fullrange = TRUE) +
              labs(title = (paste(sub("_.*", "", current_compound), 'in biofilm. Scheme',j)),x = "added standard solution [nmol]") +
              geom_text(x=15,y=1,label=paste0('R^2 = ', r2_BF[j]) ) +
              theme_bw(18)
      )
    }
    
    # in all other schemes, the regression can be done with all samples
    else {
      AS_regr <- lm(unlist(AS_area[,j]) ~stda_conc)
      AS_i_coeffs <- coefficients(AS_regr)
      AS_i_coeffs_stdv <- unlist(summary(AS_regr)[4])[3:4]
      
      BF_regr <- lm(unlist(BF_area[,j])~stda_conc)
      BF_i_coeffs <- coefficients(BF_regr)
      BF_i_coeffs_stdv <- unlist(summary(BF_regr)[4])[3:4]
      
      r2_AS[j] = format(summary(AS_regr)$r.squared,digits=3)
      r2_BF[j] = format(summary(BF_regr)$r.squared,digits=3)
      
      # plot activated sludge 
      plot(stda_conc,AS_area[,j],
           main=paste(current_compound, 'in activated sludge. Scheme', j),
           xlab='concentration of standard addition [nM]', ylab='area', xlim=c(-5,22))#, ylim=c(0,10^8))
      grid()
      abline(AS_i_coeffs)
      text(15,AS_area[1,j],paste0("R^2 = ", r2_AS[j]),1)
      lines(x,(AS_i_coeffs[2]-AS_i_coeffs_stdv[2])*x+AS_i_coeffs_stdv[1]+AS_i_coeffs[1],col="red")
      lines(x,(AS_i_coeffs[2]+AS_i_coeffs_stdv[2])*x-AS_i_coeffs_stdv[1]+AS_i_coeffs[1],col="orange")
      legend("topleft",
             c("measured points","linear regression","left stdv","right stdv"),
             col=c("black","black","red","orange"), pch = c(1,15,15,15)
      )
      
      # plot biofilm
      plot(stda_conc,BF_area[,j],
           main=paste(current_compound, 'in biofilm. Scheme', j),
           xlab='concentration of standard addition [nM]', ylab='area', xlim=c(-5,22))#, ylim=c(0,10^8))
      grid()
      abline(BF_i_coeffs)
      text(15,BF_area[1,j],paste0("R^2 = ", r2_BF[j]),1)
      lines(x,(BF_i_coeffs[2]-BF_i_coeffs_stdv[2])*x+BF_i_coeffs_stdv[1]+BF_i_coeffs[1],col="red")
      lines(x,(BF_i_coeffs[2]+BF_i_coeffs_stdv[2])*x-BF_i_coeffs_stdv[1]+BF_i_coeffs[1],col="orange")
      legend("topleft",
             c("measured points","linear regression","left stdv","right stdv"),
             col=c("black","black","red","orange"), pch = c(1,15,15,15)
      )
      this_data <- as.data.frame(t(rbind(stda_conc, BF_area[,j])))
      colnames(this_data) <- c('concentration','area')
      print(ggplot(this_data, aes(x=concentration, y=area)) +
               geom_point(color='#2980B9', size=4) +
              xlim(-10,22) +
              geom_smooth(method=lm, color='#2C3E50', fullrange = TRUE) +
              labs(title = (paste(sub("_.*", "", current_compound), 'in biofilm. Scheme',j)),x = "added standard solution [nmol]") +
              geom_text(x=15,y=1,label=paste0('R^2 = ', r2_BF[j]) ) +
              theme_bw(18)
      )
    }
    
    # calculate the concentrations based on the regression and the standard deviation of the concentrations
    AS_conc[j] <- AS_i_coeffs[1]/AS_i_coeffs[2]
    AS_conc[j][AS_conc[j] <= 0] <- 0
    AS_conc_stdv[j] <- as.numeric(AS_conc[j])*sqrt((AS_i_coeffs_stdv[1]/AS_i_coeffs[1])^2+(AS_i_coeffs_stdv[2]/AS_i_coeffs[2])^2)
    
    BF_conc[j] <- BF_i_coeffs[1]/BF_i_coeffs[2]
    BF_conc[j][BF_conc[j] <= 0] <- 0
    BF_conc_stdv[j] <- as.numeric(BF_conc[j])*sqrt((BF_i_coeffs_stdv[1]/BF_i_coeffs[1])^2+(BF_i_coeffs_stdv[2]/BF_i_coeffs[2])^2)
    
  }
  r2_AS <- as.numeric(r2_AS)
  r2_BF <- as.numeric(r2_BF)
  
  # calculate recovery
  AS_abs_recovery <- (unlist(AS_area[,3])-unlist(AS_area[,4]))/unlist(area_at_4nM)
  AS_rel_recovery <- (unlist(AS_conc[3])-unlist(AS_conc[4]))/theo_conc
  AS_rr_stdv <- sqrt(as.numeric(AS_conc_stdv[3])^2+as.numeric(AS_conc_stdv[4])^2)/theo_conc
  AS_extr_efficiency <- (unlist(AS_conc[2])-unlist(AS_conc[4]))/(unlist(AS_conc[1])-unlist(AS_conc[4]))
  AS_matrix_factor <- (unlist(AS_area[,1])-unlist(AS_area[,5]))/unlist(area_at_4nM)
  # 
  BF_abs_recovery <- (unlist(BF_area[,3])-unlist(BF_area[,4]))/unlist(area_at_4nM)
  BF_rel_recovery <- (unlist(BF_conc[3])-unlist(BF_conc[4]))/theo_conc
  BF_rr_stdv <- sqrt(as.numeric(BF_conc_stdv[3])^2+as.numeric(BF_conc_stdv[4])^2)/theo_conc
  BF_extr_efficiency <- (unlist(BF_conc[2])-unlist(BF_conc[4]))/(unlist(BF_conc[1])-unlist(BF_conc[4]))
  BF_matrix_factor <- (unlist(BF_area[,1])-unlist(BF_area[,5]))/unlist(area_at_4nM)
  

  for (k in 1:6) {
    the_output <- rbind(the_output, tibble('Compound name' = current_compound, 'AS_abs_recovery' = AS_abs_recovery[k], 'AS_rel_recovery' = AS_rel_recovery, 'AS_extr_efficiency' = AS_extr_efficiency, 'AS_matrix_factor' = AS_matrix_factor[k], 'AS_R^2' = r2_AS[k], 'BF_abs_recovery' = BF_abs_recovery[k], 'BF_rel_recovery' = BF_rel_recovery, 'BF_extr_efficiency' = BF_extr_efficiency, 'BF_matrix_factor' = BF_matrix_factor[k], 'BF_R^2' = r2_BF[k]))
  }
  # Attention, R2 are for each scheme not for each standard addition, therefore actually not correct to have them in the same table -> be careful when interpreting!
  the_output <- rbind(the_output, tibble('Compound name' = 'mean', 'AS_abs_recovery' = mean(unlist(AS_abs_recovery[c(1:4,6)])), 'AS_rel_recovery' = AS_rel_recovery, 'AS_extr_efficiency' = AS_extr_efficiency, 'AS_matrix_factor' = mean(unlist(AS_matrix_factor[c(1:4,6)])), 'AS_R^2' = mean(r2_AS), 'BF_abs_recovery' = mean(unlist(BF_abs_recovery)), 'BF_rel_recovery' = BF_rel_recovery, 'BF_extr_efficiency' = BF_extr_efficiency, 'BF_matrix_factor' = mean(unlist(BF_matrix_factor)), 'BF_R^2' = mean(r2_BF)))
  the_output <- rbind(the_output, tibble('Compound name' = 'standard deviation', 'AS_abs_recovery' = sd(unlist(AS_abs_recovery[c(1:4,6)])), 'AS_rel_recovery' = 0, 'AS_extr_efficiency' = 0, 'AS_matrix_factor' = sd(unlist(AS_matrix_factor[c(1:4,6)])), 'AS_R^2' = sd(r2_AS), 'BF_abs_recovery' = sd(unlist(BF_abs_recovery)), 'BF_rel_recovery' = 0, 'BF_extr_efficiency' = 0, 'BF_matrix_factor' = sd(unlist(BF_matrix_factor)), 'BF_R^2' = sd(r2_BF)))
  the_output[ nrow(the_output) + 1 , ] <- NA
  
  cc_output <- paste(sub("_.*", "", current_compound), 'activated sludge,', format(round(mean(unlist(AS_abs_recovery[c(1:4,6)])),3),nsmall=3), '+/-', format(round(sd(unlist(AS_abs_recovery[c(1:4,6)])),3),nsmall=3), ',', format(round(AS_rel_recovery,3),nsmall=3), '+/-', format(round(AS_rr_stdv,3),nsmall=3), ',', format(round(AS_extr_efficiency,3),nsmall=3), '+/-', format(round(sd(unlist(AS_extr_efficiency[c(1:4,6)])),3),nsmall=3), ',', format(round(mean(unlist(AS_matrix_factor[c(1:4,6)])),3),nsmall=3), '+/-', format(round(sd(unlist(AS_matrix_factor[c(1:4,6)])),3),nsmall=3),',', format(round(mean(r2_AS),3),nsmall=3), '+/-', format(round(sd(r2_AS),3),nsmall=3))
  cc_output <- rbind(cc_output, paste(sub("_.*", "", current_compound), 'biofilm,', format(round(mean(unlist(BF_abs_recovery)),3), nsmall=3), '+/-', format(round(sd(unlist(BF_abs_recovery)),3),nsmall=3), ',', format(round(BF_rel_recovery,3),nsmall=3), '+/-', format(round(BF_rr_stdv,3),nsmall=3),',', format(round(BF_extr_efficiency,3),nsmall=3), '+/-', format(round(BF_extr_efficiency,3),nsmall=3), ',', format(round(mean(unlist(BF_matrix_factor)),3),nsmall=3), '+/-', format(round(sd(unlist(BF_matrix_factor)),3),nsmall=3),',', format(round(mean(r2_AS),3),nsmall=3), '+/-', format(round(sd(r2_BF),3),nsmall=3)))
  finaloutput <- rbind(finaloutput,the_output)
  final_cc_output <- rbind(final_cc_output,cc_output)
}

# write the results to a csv
write.csv(finaloutput,"Rdata\\stdaoutput.csv", row.names = FALSE)
write.csv(final_cc_output,"Rdata\\stda_output_condensed.csv",row.names = FALSE)