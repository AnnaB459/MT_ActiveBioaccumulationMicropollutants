# Boxplot Recoveries Anna
# 30.6.22

# install libraries
library(readr) # to read csv-file
library(tidyverse) # contains data science packages
library(naniar) # used for handling missing values
library(fs) # used for handling files 
library(dplyr) # contains many data manipulation functions, useful for dataframes
library(tidyr) # even more functions for data analysis

# read file
filename_ISTD <- "Rdata\\recoveries\\ISTD_recovery.csv"
ISTD_recovery <- read_csv(file = filename_ISTD, show_col_types = TRUE)
filename_stda <- "Rdata\\recoveries\\stda_recovery.csv"
stda_recovery <- read_csv(file = filename_stda, show_col_types = TRUE)

stda_only_IS <- read_csv(file="Rdata\\recoveries\\stda_recovery_only_IS.csv", show_col_types = TRUE)

# absolute recovery both methods ---------------------------------------------------------------------------------

# make a boxplot containing all points
all_ar <- c(ISTD_recovery$ar_mean,stda_recovery$ar_mean)
boxplot(all_ar,
        main='Absolute recoveries',
        ylab='absolute recovery [%]',
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.7)

stripchart(all_ar,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8)

# make a boxplot for absolute recovery
boxplot(ISTD_recovery$ar_mean, stda_recovery$ar_mean,
        main="Absolute recovery comparing methods",
        ylab='absolute recovery [%]',
        names=c("ISTD method","standard addition method"),
        col=c("orange","yellow"),
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.7
        )

means <- c(mean(ISTD_recovery$ar_mean),mean(stda_recovery$ar_mean))

stripchart(ISTD_recovery$ar_mean,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8)


stripchart(stda_recovery$ar_mean,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8, at=2)
points(x=1:2, means, col="red", pch=19)

# make some statistical tests
# normal distribution in samples?
ISTD_ar_normal <- shapiro.test(ISTD_recovery$ar_mean) # answer is yes
stda_ar_normal <- shapiro.test(stda_recovery$ar_mean) # answer is yes
# variance similarly distributed in samples?
var_ar <- var.test(ISTD_recovery$ar_mean, stda_recovery$ar_mean, alternative = "two.sided") # answer is yes
# therefore, can check if the sample means are different
ar_same_mean <- t.test(ISTD_recovery$ar_mean, stda_recovery$ar_mean) # yes, they have the same mean


# relative recovery both methods ---------------------------------------------------------------------------------
# make a boxplot for relative recovery
boxplot(ISTD_recovery$rr_mean, stda_recovery$rr_mean,
        main="Relative recovery comparing methods",
        ylab='relative recovery [%]',
        #ylim=c(-40,150),
        names=c("ISTD method","standard addition method"),
        col=c("orange","yellow"),
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.7
)

stripchart(ISTD_recovery$rr_mean,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8)

stripchart(stda_recovery$rr_mean,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8, at=2)
abline(100,0)

# make some statistical tests
# normal distribution in samples?
ISTD_rr_normal <- shapiro.test(ISTD_recovery$rr_mean) # answer is no
stda_rr_normal <- shapiro.test(stda_recovery$rr_mean) # answer is no
# variance similarly distributed in samples?
var_rr <- var.test(ISTD_recovery$rr_mean, stda_recovery$rr_mean, alternative = "two.sided") # answer is yes
# therefore, we cannot apply a t-test to check if the sample means are different, apply a Wilcoxon text
rr_same_median <- wilcox.test(ISTD_recovery$rr_mean, stda_recovery$rr_mean, paired=FALSE, var.equal = FALSE)

# relative recovery only substances which have an IS, compared with both methods

# make a boxplot
boxplot(ISTD_recovery$rr_mean, stda_only_IS$rr_mean,
        main="Relative recovery comparing same set of substances",
        ylab='relative recovery [%]',
        #ylim=c(-40,150),
        names=c("ISTD method","standard addition method"),
        col=c("orange","yellow"),
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.7
)

stripchart(ISTD_recovery$rr_mean,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8)

stripchart(stda_only_IS$rr_mean,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8, at=2)
abline(100,0)

rr_only_IS_same_median <- wilcox.test(ISTD_recovery$rr_mean, stda_only_IS$rr_mean, paired=TRUE, var.equal = FALSE)


# absolute recovery biomass comparison ---------------------------------------------------------------------------------
# make a boxplot for absolute recovery, biofilm against activated sludge
all_AS_ar <- as.data.frame(rbind(ISTD_recovery %>% filter(str_detect(Compound, 'AS')) %>% select('ar_mean'),stda_recovery %>% filter(str_detect(Compound, 'AS')) %>% select('ar_mean')))
colnames(all_AS_ar) <- "activated sludge"
all_BF_ar <- as.data.frame(rbind(ISTD_recovery %>% filter(str_detect(Compound, 'BF')) %>% select('ar_mean'),stda_recovery %>% filter(str_detect(Compound, 'BF')) %>% select('ar_mean')))
colnames(all_BF_ar) <- "biofilm"
AS_BF_ar <- cbind(all_AS_ar,all_BF_ar)
boxplot(AS_BF_ar,
        main="Absolute recovery comparing biomass type",
        ylab='absolute recovery [%]',
        col=c("tan4","green4"),
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.7
)

stripchart(all_AS_ar,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8)

stripchart(all_BF_ar,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8, at=2)
abline(100,0)

# make a paired t-test
# first check if difference is normally distributed
all_ar_dif <- all_AS_ar[,1]-all_BF_ar[,1] # difference between both groups
shapiro.test(all_ar_dif) #normality of difference -> no
# cannot use paired t-test -> use Wilcox test
if_dif_ar_as_bf <- wilcox.test(all_AS_ar[,1],all_BF_ar[,1],paired=TRUE,alternative="two.sided")


# relative recovery biomass comparison  ---------------------------------------------------------------------------------
# make a boxplot for relative recovery, biofilm against activated sludge
all_AS_rr <- as.data.frame(rbind(ISTD_recovery %>% filter(str_detect(Compound, 'AS')) %>% select('rr_mean'),stda_recovery %>% filter(str_detect(Compound, 'AS')) %>% select('rr_mean')))
colnames(all_AS_rr) <- "activated sludge"
all_BF_rr <- as.data.frame(rbind(ISTD_recovery %>% filter(str_detect(Compound, 'BF')) %>% select('rr_mean'),stda_recovery %>% filter(str_detect(Compound, 'BF')) %>% select('rr_mean')))
colnames(all_BF_rr) <- "biofilm"
AS_BF_rr <- cbind(all_AS_rr,all_BF_rr)
boxplot(AS_BF_rr,
        main="Relative recovery comparing biomass type",
        ylab='relative recovery [%]',
        col=c("tan4","green4"),
        cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.7
)

stripchart(all_AS_rr,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8)

stripchart(all_BF_rr,
           vertical = TRUE, 
           method = "jitter", add = TRUE, pch = 20, cex=0.8, at=2)
abline(100,0)


# make a paired t-test
# first check if difference is normally distributed
all_rr_dif <- all_AS_rr[,1]-all_BF_rr[,1] # difference between both groups
shapiro.test(all_rr_dif) #normality of difference -> no
# cannot use paired t-test, use Wilcoxon test
if_dif_rr_as_bf <- wilcox.test(all_AS_rr[,1],all_BF_rr[,1],paired=TRUE,alternative="two.sided")
