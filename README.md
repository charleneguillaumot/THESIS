# R scripts for Southern Ocean SDM
Provide codes for applying Species Distribution Models in Southern Ocean case studies 

Refers to published papers
---------------------------
Guillaumot C, A Martin, M Eléaume & T Saucède (2016) ‘SDMPlay’: Species Distribution Modelling Playground, CRAN. https://cran.r-project.org/web/packages/SDMPlay 

Guillaumot C, Artois J, Saucède T, Demoustier L, Moreau C, Eléaume M, Agüera A, Danis B. submitted. Species distribution models in a data-poor and broad scale context. Progress in Oceanography


# delim.area (Guillaumot et al. 2016)
This function helps to crop the environmental descriptors dataset used to calibrate the BRT models (RasterStack) in latitude, longitude and bathymetry interval.

# CV.BRT and Function_ja.R codes
Details of the codes used in Guillaumot et al. (submitted in Progress in Oceanography) to calibrate random and spatial cross-validation procedures.

# CLOCK-2 function 
Detail and description of the script used to generate the spatial 2-fold CLOCK cross-validation procedure (Guillaumot et al. submitted). See illustration B in the jpeg file. 

# CLOCK-3 function 
Detail and description of the script used to generate the spatial 3-fold CLOCK cross-validation procedure (Guillaumot et al. submitted). 
See illustration C in the jpeg file. 

# CLOCK-6 function 
Detail and description of the script used to generate the spatial 6-fold CLOCK cross-validation procedure (Guillaumot et al. submitted). 
See illustration D in the jpeg file. 

# CLOCK-4 function
A similar code was added to generate a 4-fold CLOCK spatial cross-validation procedure

# Folder "Extrapolation paper" 
Refers to example code of the paper Guillaumot et al. (2020) "Extrapolation in species distribution modelling. Application to Southern Ocean marine species" submitted to Progress in Oceanography. It contains lines of codes and data to run BRT model and calculate the proportion of model predictions that is actually extrapolation. It also contains a part in which the percentage of contribution of each environmental layer to extrapolation is calculated.  
