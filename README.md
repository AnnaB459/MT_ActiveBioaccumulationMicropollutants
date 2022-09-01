# MT_ActiveBioaccumulationMicropollutants

## Master thesis
This is the repository of all my code files I used for my master thesis.

### Abstract of the master thesis: 
Synthetic micropollutants are omnipresent in wastewaters and natural waters where they can have detrimental effects on ecosystems. Microbial communities such as biofilm and activated sludge interact with micropollutants in various ways, including accumulation: a term encompassing both passive adsorption and active bioaccumulation. In this master thesis, the QuEChERS extraction method was used in order to determine the amount of actively bioaccumulated compounds in biofilm and activated sludge. The firstly performed recovery experiment showed that out of 65 substances, 74% could be extracted with sufficient recoveries between 70% and 130% from both tested microbial communities. The following bioaccumulation experiment conducted with biofilm grown in river Ticino and activated sludge from wastewater treatment plant Airolo revealed that 20 of these substances accumulated to at least 10% of the initial amount in at least one of the two microbial communities (passive adsorption and/or active accumulation). Further, eleven of these substances showed active bioaccumulation to at least 5%: Four substances actively bioaccumulated only in biofilm, six substances only in activated sludge and one substance in both microbial communities. Several of the bioaccumulated compounds have already shown bioaccumulation in activated sludge and/or biofilm in other studies. Eight of the actively bioaccumulated compounds contain aliphatic amine moieties, which could facilitate active bioaccumulation. 


Les micropolluants de synthèse sont omniprésents dans les eaux polluées et les eaux de surface où ils peuvent avoir des effets néfastes sur les écosystèmes. Les communautés microbiennes comme les biofilms ou les boues activées interagissent avec les micropolluants de différentes manières, y compris en les accumulant. Le terme accumulation inclut la sorption passive ainsi que la bioaccumulation active. Dans ce projet de master, la méthode d'extraction "QuEChERS" a été mise en œuvre afin de déterminer la quantité de micropolluants bioaccumulés activement dans les biofilms et les boues activées. Les expériences de récupération ont montré que, parmi 65 substances, 74% ont pu être extraites des deux communautés microbiennes avec des taux de récupération allant de 70% à 130%. Ensuite, des expériences de bioaccumulation ont été effectuées avec du biofilm de la rivière du Tessin et avec de la boue activée de la station d'épuration d'Airolo. Les résultats ont montré que 20 substances ont été accumulées avec des taux d'au moins 10% de la quantité initiale dans au moins une des communautés microbiennes (sorption passive et/ou bioaccumulation active). De plus, onze de ces substances ont montré une bioaccumulation active d'au moins 5%: quatre substances ont été accumulées activement uniquement dans le biofilm, six substances ont été accumulées activement uniquement dans les boues activées et une substance a été accumulée activement dans les deux communautés microbiennes. Plusieurs de ces substances ont déjà montré une bioaccumulation dans d'autres études. Huit substances parmi les onze substances accumulées activement contiennent des groupes amines aliphatiques ce qui pourrait faciliter l'accumulation active.


## Link to master thesis
... will put it here later

## Description of files
07_12_recovery_analysis_with_ISTD.R: calculation of recoveries for substances with an internal standard

Recovery_analysis_with_STDA.R: calculation of recoveries for substances evaluated with the standard addition method

boxplot_recoveries_Anna.R: graphical representation of recoveries in boxplots, plus statistical tests to determine differences between subsamples

07_13_biotransf_quechers.R: calculation and graphic representation of the bioaccumulation, needs function-file functions_for_quechers.R to work properly. 

functions_for_quechers.R: all functions needed in the file 07_13_biotransf_quechers.R

mass_balance.R: calculation of a mass balance of the bioaccumulation experiment
