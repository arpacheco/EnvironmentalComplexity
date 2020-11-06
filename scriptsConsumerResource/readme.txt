There are two MATLAB script files that are used to generate consumer resource modeling data:

1. solveCRM.m: for experimentally-tested organisms. Parameterizes C matrix based on experimentally-derived nutrient utilization data, and parameterizes D matrix based on flux-balance predictions using genome scale models.

2. solveCRMGeneric.m: for statistical ensembles of simulated organisms. Parameterizes C matrix probabilistically based on generalist-specialist designations for each organism (see Methods), and D matrix based on probability of transforming a given nutrient. Script is set to simulate community dynamics with increasing pool of available secreted metabolites, as well as with increasing number of provided primary nutrients.

Both these scripts depend on CRM.m, which forms the mathematical basis for the consumer resource model.