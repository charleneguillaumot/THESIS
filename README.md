# R scripts for Southern Ocean SDM
Provide codes for applying Species Distribution Models and Bayesian integrated approach in Southern Ocean case studies 
Refers to codes used in published papers from my phD thesis (2017-2021)

# Folder "Integrated approach" 
From Guillaumot C, Buba Y, Belmaker J, Fourcy D, Danis B, Dubois P, Saucède T (submission ongoing). Simple or hybrid ? Next generation ecological models to study the distribution of Southern Ocean marine species. Ecology Letters.
The codes to run the different methods developed in the paper are provided in this folder: simple SDM (GLM), spatial projection of the DEB, integration of the DEB within a SDM, and Bayesian integrated approach
All files are available (environment, presence records...) to entirely run the different models.


# Folder "Extrapolation paper" 
From Guillaumot C, Moreau C, Danis B, Saucède T (2020). Extrapolation in species distribution modelling. Application to Southern Ocean marine species. Progress in Oceanography. 188, 102438. https://doi.org/10.1016/j.pocean.2020.102438

It contains lines of codes and data to run BRT model and calculate the proportion of model predictions that is model extrapolation. It also contains a part in which the percentage of contribution of each environmental layer to extrapolation is calculated. Please send an email at charleguillaumot21@gmail.com if you want the files bias_OK.asc and stack_intermediate.grd, that are missing to run the code properly.


# Folder "Choice of predictors" 
From Guillaumot C, Danis B, Saucède T (2020). Selecting environmental descriptors is critical to modelling the distribution of Antarctic benthic species. Polar Biology. 1-19.
https://doi.org/10.1007/s00300-020-02714-2

It contains the two files that enable to calculate "extreme events" as described in the supplementary material of the paper.


# Folder "Spatial cross validation procedures"
From Guillaumot C, Artois J, Saucède T, Demoustier L, Moreau C, Eléaume M, Agüera A, Danis B (2019). Broad-scale species distribution models applied to data-poor areas. Progress in Oceanography, 175, 198-207. https://doi.org/10.1016/j.pocean.2019.04.007

CV.BRT and Function_ja.R codes
-> To calibrate random and spatial cross-validation procedures
CLOCK-2 function 
-> Details and describes how to generate the spatial 2-fold CLOCK cross-validation procedure. See illustration B in the jpeg file also contained in the folder. 
CLOCK-3 function 
-> Details and describes how to generate the spatial 3-fold CLOCK cross-validation procedure. See illustration C in the jpeg file also contained in the folder. 
CLOCK-4 function
-> Details and describes how to generate the spatial 4-fold CLOCK cross-validation procedure. 
CLOCK-6 function
-> Details and describes how to generate the spatial 6-fold CLOCK cross-validation procedure. See illustration D in the jpeg file also contained in the folder. 

# delim.area file (Guillaumot et al. 2016)
From Guillaumot C, A Martin, M Eléaume & T Saucède (2016) ‘SDMPlay’: Species Distribution Modelling Playground, CRAN. https://cran.r-project.org/web/packages/SDMPlay 
This function helps to crop the environmental descriptors dataset used to calibrate BRT models (RasterStack) in latitude, longitude and bathymetry interval.


