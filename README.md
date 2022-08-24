# MT_ActiveBioaccumulationMicropollutants

## Master thesis
This is the repository of all my code files I used for my master thesis.

Abstract of the master thesis: Synthetic micropollutants are omnipresent in wastewaters and natural waters where they can have detrimental effects on ecosystems. Microbial communities such as biofilm and activated sludge interact with micropollutants in various ways, including accumulation: a term encompassing both passive sorption and active bioaccumulation. In this master thesis, the QuEChERS extraction method was used in order to determine the amount of actively bioaccumulated compounds in biofilm and activated sludge. The firstly performed recovery experiment showed that out of 65 substances, 74% could be extracted with sufficient recoveries between 70% and 130% from both tested microbial communities. The following bioaccumulation experiment conducted with biofilm grown in river Ticino and activated sludge from wastewater treatment plant Airolo revealed that 20 of these substances accumulated to at least 10\% of the initial amount in at least one of the two microbial communities (passive sorption and / or active accumulation). Further, 10 of these substances showed active bioaccumulation to at least 5\%: Four substances actively bioaccumulated only in biofilm, five substances only in activated sludge and one substance in both microbial communities. Several of the bioaccumulated compounds have already shown bioaccumulation in activated sludge and/or biofilm in other studies. All actively bioaccumulated compounds contain amine moieties.

## Link to master thesis
... will put it here later

## Description of files
07_12_recovery_analysis_with_ISTD.R: calculation of recoveries for substances with an internal standard

Recovery_analysis_with_STDA.R: calculation of recoveries for substances evaluated with the standard addition method

boxplot_recoveries_Anna.R: graphical representation of recoveries in boxplots, plus statistical tests to determine differences between subsamples

07_13_biotransf_quechers.R: calculation and graphic representation of the bioaccumulation, needs function-file functions_for_quechers.R to work properly. 

functions_for_quechers.R: all functions needed in the file 07_13_biotransf_quechers.R

mass_balance.R: calculation of a mass balance of the bioaccumulation experiment
