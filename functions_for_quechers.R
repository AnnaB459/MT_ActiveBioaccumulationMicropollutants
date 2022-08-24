# function IS_concentrations: concentration calculation for substances with IS
# input: all filenames, calculated amount and weight of a compound in one of the two matrices (biofilm or
# activated sludge) and either sorption or accumulation batch -> the function is called four times
# calculation: the function calculates the actual accumulation / sorption amount by normalizing it to the
# respective biomass weight and making a mean throughout the replicates
IS_concentrations <- function(compound_data) {
  # initialize lists
  C48 <- list()
  C96 <- list()
  # calculate for each measurement the weight substance per weight biomass
  compound_data$c_per_m <- compound_data$`Concentration in ng/L`/compound_data$`biomass weight in g` *10^-3 # concentration per mass in ng substance / g biomass. Factor *10^-3 because it initially was in 1mL (and substance conc in ng/L)
  
  # calculate C0
  C0 <- compound_data %>% filter(str_detect(Filename, '0_r1')) %>% select('c_per_m')
  
  # calculate other concentrations (C48 and C96)
  
  replicates <- c('_r1','_r2','_r3','_r4','_r5','_r6')
  for (i in 1:length(replicates)) {
    detect_replicates <- str_detect(compound_data$Filename,replicates[i]) # gives TRUE / FALSE vector
    pos_r <- which(detect_replicates) # gives index of true in detect_replicates, position_replicate
    if (length(pos_r)==0) {next}
    if (pos_r[1]>6) {next} # if there is no measurement at the first index, i.e. at time 0
    if (pos_r[2]>12 || is.na(pos_r[2])) {C48 <- cbind(C48,NA)
                      C96 <- cbind(C96, compound_data$c_per_m[pos_r[2]] - compound_data$c_per_m[pos_(1)])
                      } else {C48 <- cbind(C48,compound_data$c_per_m[pos_r[2]] - compound_data$c_per_m[pos_r[1]])
                              C96 <- cbind(C96,compound_data$c_per_m[pos_r[3]] - compound_data$c_per_m[pos_r[1]]) }
  }
  pos_r <- list() # empty pos_r
  # calculate means and standard deviations
  C48_mean <- mean(unlist(C48), na.rm=TRUE)
  if (C48_mean <=0) {C48_mean <- 0}
  C96_mean <- mean(unlist(C96), na.rm=TRUE)
  if (C96_mean <=0) {C96_mean <- 0}
  C48_stdv <- sd(C48, na.rm=TRUE)
  C96_stdv <- sd(C96, na.rm=TRUE)
  # return mean in ng substance / g biomass and standard deviation for the current substance
  return(c(C0,C48_mean,C48_stdv,C96_mean,C96_stdv))
  
}

# function stda_function: concentration calculation for substances with stda
# input: all filenames, area and weight of a compound in one of the two matrices (biofilm or
# activated sludge) and either sorption or accumulation batch -> the function is called four times
# calculation: the function calculates the actual accumulation / sorption amount by first calculating the concentrations from the standard addition,
# then normalizing it to the respective biomass weight and calculate the standard deviation due to the linear regression
stda_function <- function(compound_data) {
  
  # for each sample, add the corresponding standard addition to the dataframe and divide by mass to have conc/mass
  stda_conc <- c(0, 2, 2, 10, 20, 20) # nmol/L
  compound_data$add_p_m <- NA
  
  # associate the added standard addition amount to each sample
  for (i in 1:dim(compound_data)[1]) {
    if (grepl('_r1', compound_data$Filename[i])) {
      compound_data$add_p_m[i] <- stda_conc[1]/compound_data$`biomass weight in g`[i]
    }
    else if (grepl('_r2', compound_data$Filename[i]) || grepl('_r3', compound_data$Filename[i])) {
      compound_data$add_p_m[i] <- stda_conc[2]/compound_data$`biomass weight in g`[i]
    }
    else if (grepl('_r4', compound_data$Filename[i])) {
      compound_data$add_p_m[i] <- stda_conc[4]/compound_data$`biomass weight in g`[i]
    }
    else if (grepl('_r5', compound_data$Filename[i]) || grepl('_r6', compound_data$Filename[i])) {
      compound_data$add_p_m[i] <- stda_conc[5]/compound_data$`biomass weight in g`[i]
    }
  }
  
  # now make a linear regression for each timepoint
  timepoints <- c('_0_', '_48_','_96_')
  conc <-list()
  conc_stdv <- list()
  regr_R2 <- list()
  # loop through all timepoints (separate linear regression for each timepoint)
  for (i in 1:length(timepoints)) {
    # linear regression
    regr <- lm(unlist(compound_data %>% filter(str_detect(Filename, timepoints[i])) %>% select('Area')) ~ unlist(compound_data %>% filter(str_detect(Filename, timepoints[i])) %>% select('add_p_m')))
    regr_coeffs <- coefficients(regr)
    regr_stdv <- unlist(summary(regr)[4])[3:4]
    regr_R2[i] = format(summary(regr)$r.squared,digits=3)
    
    # uncomment following lines to see how well the fit of the linear regression is
    #plot(unlist(compound_data %>% filter(str_detect(Filename, timepoints[i])) %>% select('add_p_m')),unlist(compound_data %>% filter(str_detect(Filename, timepoints[i])) %>% select('Area')))
    #abline(regr_coeffs)
    
    # calculate the concentrations
    conc[i] <- regr_coeffs[1]/regr_coeffs[2]
    conc[i][conc[i] <= 0] <- 0
    conc_stdv[i] <- as.numeric(conc[i])*as.numeric(sqrt((regr_stdv[1]/regr_coeffs[1])^2+(regr_stdv[2]/regr_coeffs[2])^2))
  }

  # return the concentrations for each timepoint in nmol/g
  return(c(as.numeric(conc[1])*10^-3,as.numeric(conc[2])*10^-3,as.numeric(conc_stdv[2])*10^-3,as.numeric(conc[3])*10^-3,as.numeric(conc_stdv[3])*10^-3,as.numeric(conc_stdv[1])*10^-3,regr_R2)) # converted to unit nmol/g
  
}

# function bioaccumulation calculates the bioaccumulated amount in nmol/g
# input: data-frame with all the concentrations in the different batches (should already be in nmol/g)
# done (as concentration calculation before) for one substance after the other
bioaccumulation <- function(conc_frame) {
  # calculate the bioaccumulation with its respective standard-deviations for biofilm
  accum0_bf <- conc_frame[1,2] - abs(conc_frame[2,2])
  accum0_bf_sd <- sqrt(conc_frame[1,7]^2+conc_frame[2,7]^2)
  accum48_bf <- conc_frame[1,3] - abs(conc_frame[2,3])
  accum48_bf_sd <- sqrt(conc_frame[1,4]^2+conc_frame[2,4]^2)
  accum96_bf <- conc_frame[1,5] - abs(conc_frame[2,5])
  accum96_bf_sd <- sqrt(conc_frame[1,6]^2+conc_frame[2,6]^2)
  
  # do the same for activated sludge
  accum0_as <- conc_frame[3,2] - abs(conc_frame[4,2])
  accum0_as_sd <- sqrt(conc_frame[3,7]^2+conc_frame[4,7]^2)
  accum48_as <- conc_frame[3,3] - abs(conc_frame[4,3])
  accum48_as_sd <- sqrt(conc_frame[3,4]^2+conc_frame[4,4]^2)
  accum96_as <- conc_frame[3,5] - abs(conc_frame[4,5])
  accum96_as_sd <- sqrt(conc_frame[3,6]^2+conc_frame[4,6]^2)
  
  # make dataframe to return
  all_accumulation <- as.data.frame(rbind(c(accum0_bf,accum48_bf,accum48_bf_sd,accum96_bf,accum96_bf_sd,accum0_bf_sd),c(accum0_as,accum48_as,accum48_as_sd,accum96_as,accum96_as_sd,accum0_as_sd)))
  return(all_accumulation)
}

# this function makes bar plots for t48 and t96, showing the amount in each batch
# input: the final concentrations and the name of the current compound
# output: no output, produces plot
bar_plotting <- function(final_conc,current_compound) {
  library(ggplot2)
  biocommunity <- c("BF","AS")
  mean_weight <- c(0.60375,0.468)
  bc_colors <- c("green4","tan4")
  
  for (i in 1:length(biocommunity)) {
    active_accum <- c(0,(final_conc[((2*i)-1),3]-final_conc[(2*i),3]),0,(final_conc[((2*i)-1),5]-final_conc[(2*i),5]))
    active_accum[active_accum <= 0] <- 0
    stdv <- c(final_conc[((2*i)-1):(2*i),4],final_conc[((2*i)-1):(2*i),6])
    value=c(final_conc[((2*i)-1):(2*i),3],final_conc[((2*i)-1):(2*i),5],active_accum)
    miny = c(((value[1:4]-stdv)[1:4]),0,0,0,0)
    maxy = c(((value[1:4]+stdv)[1:4]),0,0,0,0)
    #print(paste("this is miny", miny,"this is maxy",maxy))
    data <- data.frame(timepoint=c('t48','t48','t96','t96','t48','t48','t96','t96'),
               reactor=c('BE','SC', 'BE','SC','BE','SC', 'BE','SC'),
               accumulation=c('total (measured)','total (measured)','total (measured)','total (measured)','min active (calculated)','min active (calculated)','min active (calculated)','min active (calculated)'),
               value=c(final_conc[((2*i)-1):(2*i),3],final_conc[((2*i)-1):(2*i),5],active_accum))
    # now make the plots
    # where to automatically save all plots "eawag","userdata","boeselan","My Documents","Rdata","PlotOutputs","Accumulation", 
    mypath <- file.path("Rdata","PlotOutputs","Accumulation",paste("accumulation_", sub("_.*", "", current_compound), "_", biocommunity[i], ".jpeg", sep=''))
    jpeg(file=mypath)
    #pdfpath <- file.path("Rdata","PlotOutputs","Accumulation","justapdf.pdf")
    #pdf(file=pdfpath)
    print(ggplot(data, aes(x=reactor, y=value, pattern=accumulation)) +
      facet_grid(~timepoint) +
      geom_bar(stat = 'identity') +
      geom_col_pattern(
        position="stack",
        aes(reactor, value, pattern_fill = accumulation), 
        pattern = c('none','none','none','stripe','none','none','none','stripe'),
        fill    = c(bc_colors[i],bc_colors[i],'white','white',bc_colors[i],bc_colors[i],'white','white'),
        colour  = 'black'
      ) +
      theme_bw(18) +
      theme(legend.position = 'right') +
      theme(legend.key.size = unit(1.5, 'cm')) +
      scale_y_continuous(
        name='amount [nmol/g]', #limits=c(0,4.2),  #if wanted, can set limits for graphs
        sec.axis = sec_axis(trans=~.*mean_weight[i]/5*100, name="accumulation [%]")) +
        geom_errorbar(aes(ymin=miny, ymax=maxy), position=position_dodge(.0)) +
      labs(title = paste(sub("_.*", "", current_compound), 'in solid phase of', biocommunity[i])) +
      guides(pattern = guide_legend(override.aes = list(fill = c("white",bc_colors[i])))) +
      guides(pattern_fill=guide_legend(override.aes=list(pattern = c("stripe","none"))))
    )
    dev.off()
  }
}