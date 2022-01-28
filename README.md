# terrestrial-aquatic-LCT

R code for "A terrestrial-aquatic model reveals cross-scale interactions regulate lateral carbon transport from terrestrial ecosystems"
By Talbot et al. 

Functions #Code for custom functions needed to run analyses 
  Code
    TAM_MSVersion.R #Model code 
    SummaryTAM_MSVersion.R #Calculate annual summary data from model output
    genWeather_MSVersion.R #Weather generator for model experiments
    genWeatherMOD_MSVersion.R #Modified weather generator that takes mean monthly air temp as input and generates HUC-2 weather
    gddSum_MSVersion.R #Calculate phenology with weather data based on degree days and light
    findInitials_MSVersion.R #Find initial conditions for state variables that take longerto reach equilibrium (i.e., wood C, passive soil C)  
    evaporationFunc_MSVersion.R #Calculate evaporation from aquatic surface (Not implemented in this MS) 
    
SensitivityAnalysis #Code for model runs/analyses/figures
  Code
    SensitivityAnalysisBase_MSVersion.R #Run model over climate gradients w/base parameterization 
    SA_multiCor_MSVersion.R #Run model over climate gradients twice for each parameter (+/- 10%) with parallel computing
    Analyze_SA_MSVersion.R #Summarize sensitivity analysis results and calculate mean percent change in each variable
    MakeWeqather_MSVersion.R #create weather timeseries from climate gradients
  OutputFiles
    Out #folder with simulations for each parameter and wood C/passive soil C initial conditions (uses old parameter names!)
    dfOut.csv #weather generator output for climate gradients
    phenolDF.csv #phenology based on weather generator output for climate gradients
    SASoil_base.csv #model simulations with base parameterization
    SASoil_trait_i0_cp.csv #initial conditions for passive soil C
    SASoil_trait_i0.csv #initial conditions for wood C
    
ModelExperiments #model experiment code and output files
  Code
    TheoryPlots_MSVersion.R #Create plots for LCT, NPP, water drainage, Soil C, Soil DOC, ..
    runTraits_MSVersion.R #Expanded sensitivity analysis for plant trait values that span range in PFT
  OutputFiles
    Amax_i0_cp_NewS.csv #file with passive soil C initial conditions for all values of Amax
    Amax_i0_cw_NewS.csv #file with wood C initial conditions for all values of Amax
    Amax_wide.csv #model simulations for all values of Amax
    ..same filed for each of Kwue, Topt, and aw...
    
Validation
  Code
    HUCcharacteristics_MSVersion.R  #Collate data from NHD & StreamCat for each HUC
    HUCcharacteristicsRO_MSVersion  #Collate data from NHD & StreamCat for each HUC/version with runoff ratio included
    CompareToButman_MSVersion.R #Run model for each pft for each HUC, summarize, make 1:1 plots and LCT/precip plot 
  OutputFiles
    ButmanTable.csv #table of C fluxes from Butman et al. 2016
    CompareLCT.csv #Modeled LCT fluxes vs. Butman et al. 2016
    CompareNPP.csv #Modeled NPP fluxes vs. Butman et al. 2016
    ModeledArealLCT.csv #modeled areal LCT fluxes
    PFT_i0_cf2.csv #initial conditions for various C pools for each CONUS HUC-02
    PFT_i0_cm.csv
    PFT_i0_co.csv
    PFT_i0_cp.csv
    PFT_i0_cr.csv
    PFT_i0.csv
    regionsArea.csv #area of each CONUS HUC-02
    regionsPFT.csv #area weighted PFT coverage in each CONUS HUC-02
    regionsPhenology.csv #Phenology for each PFT in each HUC-02
    regionsRunoff.csv #area weighted runoff and runoff ratio for each HUC-02
    regionsWeather.csv #weather timeseries for each HUC based on mean T and area weighted total annual precip
    RunoffRatios.csv #modeled and calculated runoff ratios
    storedRegions.csv #model outputs for each var in each HUC

