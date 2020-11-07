CODE SUMMARY:

This directory contains two MATLAB script files that are used to generate consumer resource modeling data:

1. solveCRM.m: for experimentally-tested organisms. Parameterizes C matrix based on experimentally-derived nutrient utilization data, and parameterizes D matrix based on flux-balance predictions using genome scale models. When run as presented, this script will simulate the growth dynamics of a 4-species community made up of B. subtilis, M. extorquens, S. oneidensis, and P. aeruginosa in all 63 nutrient combinations tested in the 'com4' experiment. 

2. solveCRMGeneric.m: for statistical ensembles of simulated organisms. Parameterizes C matrix probabilistically based on generalist-specialist designations for each organism (see Methods), and D matrix based on probability of transforming a given nutrient. Script is set to simulate community dynamics with increasing pool of available secreted metabolites, as well as with increasing number of provided primary nutrients. When run as presented, this script will run simulations for a 13-species community under the assumptions of 'CRM-A' (i.e. no organisms are generalists). 

These scripts depend on the following files:
	- CRM.m: Forms the mathematical basis for the consumer resource model.
	- generateCMatrix.m: Generates species-specific C matrices, which define each organisms' nutrient utilization capabilities
	- generateDMatrix.m: Generates species-specific D matrices using flux-balance analysis.
	- defineMedium.m: used by 'generateDMatrix.m' to define medium conditions for flux balance analysis.

SYSTEM REQUIREMENTS:
	Software: 
		- MATLAB R2020a or newer.
		- The COnstraint-Based Reconstruction and Analysis (COBRA) Toolbox for MATLAB (For FBA simulations), available at https://opencobra.github.io/cobratoolbox/stable/. 
	Tested on: Mac OS X Catalina and Mojave
	Non-standard hardware requirements: None

INSTALLATION GUIDE:
	Instructions: Please refer to the MATLAB installation guide at www.mathworks.com/help/install/install-products.html, as well as the COBRA toolbox installation guide at opencobra.github.io/cobratoolbox/latest/installation.html for detailed instructions and for typical installation times.

DEMO:
	Instructions: Both 'solveCRM.m' and 'solveCRMGeneric.m' are designed to run as standalone functions within the 'scriptsConsumerResource' directory.
	Expected outputs:
		- For solveCRM.m: A matrix ('relAbus') containing relative abundance data of each condition. Additionally, if plotIndividuals is set to 1, a figure showing the time-dependent trajectory of organism and resource abundances, as well as final species abundances for the last tested nutrient condition (similar to data in Supplementary Figure 24).
		- For solveCRMGeneric.m: A data structure ('epiCRM'), which contains the simulated condition-specific community growth and taxonomic data. This structure can be saved and visualized within 'analyzeTaxonomy.m' to compare growth and taxonomy data with experimental results.
	Expected run times:
		- For solveCRM.m: ~60-90 seconds for all 63 nutrient conditions of com4.
		- For solveCRMGeneric.m:  ~30-45 minutes for up to 10 secreted metabolites with 50 random simulations each.

INSTRUCTIONS FOR USE:
	Please refer to 'DEMO' section above for instruction on reproducing results in the manuscript, as well as the header of each script for more specific details on data inputs.