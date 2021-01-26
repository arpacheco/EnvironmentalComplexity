function C = generateCMatrix(expDataFile,selectSpecies,selectNutrients)

load(expDataFile);

C = zeros(length(selectSpecies),length(selectNutrients));

for s = 1:length(selectSpecies)
    C(s,:) = stockGrowthCont.(selectSpecies{s})(find(ismember(nutrientsStock,selectNutrients)));
end