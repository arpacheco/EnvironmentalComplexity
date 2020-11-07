function [C,data] = generateCMatrix(expDataFile,selectSpecies,nutrients,selectNutrients)

%% Initialize C matrix

C = zeros(length(selectSpecies),length(selectNutrients));

%% Load and process growth data

sheetNames = {};
for i = 1:length(selectSpecies)
    sheetNames = [sheetNames,[selectSpecies{i} '-A']];
end

for s = 1:length(sheetNames)
    [~, ~, raw] = xlsread(expDataFile,sheetNames{s});
    cells = raw(size(raw,1)-7:size(raw,1),2);
    dataTemp = raw(size(raw,1)-7:size(raw,1),3:14);
    dataTemp = reshape([dataTemp{:}],size(dataTemp));
    dataTemp(find(dataTemp) == 0) = NaN;
    
    dataTemp = dataTemp-min(min(dataTemp)); % Subtract min OD
    dataTemp(find(dataTemp < 0.006)) = 0; % Correct for no growth
    
    sp = split(sheetNames{s},'-');

    if contains(sheetNames{s},'A')
        data.(char(sp(1))) = dataTemp;
    else
        data.(char(sp(1))) = [data.(char(sp(1)));dataTemp];
    end

end

%% Load and process nutrient data

nutrientCoordVec = [];
for i = 1:length(nutrients)
nutrientCoordVec = [nutrientCoordVec repelem(i,3)];
end
nutrientCoordMat = flipud(rot90(reshape(nutrientCoordVec,12,8)));
nutrientMat = nutrients(nutrientCoordMat);

%% Make C matrix

for s = 1:length(selectSpecies)
    for n = 1:length(selectNutrients)
        nutrientCoords = find(ismember(nutrientMat,selectNutrients{n}));
        C(s,n) = mean(data.(selectSpecies{s})(nutrientCoords));
    end
end