# terrestrial-aquatic-LCT

R code for "A terrestrial-aquatic model reveals cross-scale interactions regulate lateral carbon transport from terrestrial ecosystems"
By Talbot et al. 

Code for custom functions needed to run analyses 
Functions-folder
TAM_MSVersion.R #Model code 
SummaryTAM_MSVersion.R #Calculate annual summary data from model output
genWeather_MSVersion.R #Weather generator
genWeatherMOD_MSVersion.R #Modified weather generator that takes mean monthly air temp as input
gddSum_MSVersion.R #Calculate phenology with weather data
findInitials_MSVersion.R #Find initial conditions for state variables that take longer than run years to reach equilibrium
evaporationFunc_MSVersion.R #Calculate aquatic evaporation 

Code for model runs/analyses/figures
SensitivityAnalysis-folder
SensitivityAnalysisBase_MSVersion.R #Run model over climate gradients w/base parameter$
SA_multiCor_MSVersion.R #Run model over climate gradients twice for each parameter (-1$
Analyze_SA_MSVersion.R #Summarize sensitivity analysis results and calculate mean perc$

ModelExperiments-folder
TheoryPlots_MSVersion.R #Create plots for LCT, NPP, water drainage, Soil C, Soil DOC, $
runTraits_MSVersion.R #Expanded sensitivity analysis for plant trait values that span $

Validation-folder
HUCcharacteristics_MSVersion.R  #Collate data from NHD & StreamCat for each HUC
CompareToButman_MSVersion.R #Run model for each pft for each HUC, summarize, make 1:1 plots and LCT/precip plot 

