# terrestrial-aquatic-LCT

R code for "A terrestrial-aquatic model reveals cross-scale interactions regulate lateral carbon transport from terrestrial ecosystems"
By Talbot et al. 

Functions #Code for custom functions needed to run analyses <br />
  Code<br />
    TAM_MSVersion.R #Model code <br />
    SummaryTAM_MSVersion.R #Calculate annual summary data from model output<br />
    genWeather_MSVersion.R #Weather generator for model experiments<br />
    genWeatherMOD_MSVersion.R #Modified weather generator that takes mean monthly air temp as input and generates HUC-2 weather<br />
    gddSum_MSVersion.R #Calculate phenology with weather data based on degree days and light<br />
    findInitials_MSVersion.R #Find initial conditions for state variables that take longerto reach equilibrium (i.e., wood C, passive soil C)  <br />
    evaporationFunc_MSVersion.R #Calculate evaporation from aquatic surface (Not implemented in this MS) <br />
    <br />
SensitivityAnalysis #Code for model runs/analyses/figures<br />
  Code<br />
    SensitivityAnalysisBase_MSVersion.R #Run model over climate gradients w/base parameterization <br />
    SA_multiCor_MSVersion.R #Run model over climate gradients twice for each parameter (+/- 10%) with parallel computing <br />
    Analyze_SA_MSVersion.R #Summarize sensitivity analysis results and calculate mean percent change in each variable<br />
    MakeWeqather_MSVersion.R #create weather timeseries from climate gradients<br />
  OutputFiles<br />
    Out #folder with simulations for each parameter and wood C/passive soil C initial conditions (uses old parameter names!)<br />
    dfOut.csv #weather generator output for climate gradients<br />
    phenolDF.csv #phenology based on weather generator output for climate gradients<br />
    SASoil_base.csv #model simulations with base parameterization<br />
    SASoil_trait_i0_cp.csv #initial conditions for passive soil C<br />
    SASoil_trait_i0.csv #initial conditions for wood C<br />
    
ModelExperiments #model experiment code and output files<br />
  Code<br />
    TheoryPlots_MSVersion.R #Create plots for LCT, NPP, water drainage, Soil C, Soil DOC, ..<br />
    runTraits_MSVersion.R #Expanded sensitivity analysis for plant trait values that span range in PFT<br />
  OutputFiles<br />
    Amax_i0_cp_NewS.csv #file with passive soil C initial conditions for all values of Amax <br />
    Amax_i0_cw_NewS.csv #file with wood C initial conditions for all values of Amax<br />
    Amax_wide.csv #model simulations for all values of Amax<br />
    ..same filed for each of Kwue, Topt, and aw...<br />
    
Validation<br />
  Code<br />
    HUCcharacteristics_MSVersion.R  #Collate data from NHD & StreamCat for each HUC<br />
    HUCcharacteristicsRO_MSVersion  #Collate data from NHD & StreamCat for each HUC/version with runoff ratio included<br />
    CompareToButman_MSVersion.R #Run model for each pft for each HUC, summarize, make 1:1 plots and LCT/precip plot <br />
  OutputFiles<br />
    ButmanTable.csv #table of C fluxes from Butman et al. 2016<br />
    CompareLCT.csv #Modeled LCT fluxes vs. Butman et al. 2016<br />
    CompareNPP.csv #Modeled NPP fluxes vs. Butman et al. 2016<br />
    ModeledArealLCT.csv #modeled areal LCT fluxes<br />
    PFT_i0_cf2.csv #initial conditions for various C pools for each CONUS HUC-02<br />
    PFT_i0_cm.csv<br />
    PFT_i0_co.csv<br />
    PFT_i0_cp.csv<br />
    PFT_i0_cr.csv<br />
    PFT_i0.csv<br />
    regionsArea.csv #area of each CONUS HUC-02<br />
    regionsPFT.csv #area weighted PFT coverage in each CONUS HUC-02<br />
    regionsPhenology.csv #Phenology for each PFT in each HUC-02<br />
    regionsRunoff.csv #area weighted runoff and runoff ratio for each HUC-02<br />
    regionsWeather.csv #weather timeseries for each HUC based on mean T and area weighted total annual precip<br />
    RunoffRatios.csv #modeled and calculated runoff ratios<br />
    storedRegions.csv #model outputs for each var in each HUC<br />

