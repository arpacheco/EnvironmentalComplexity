% analyzeTaxonomy
%
% Processes 16s read data from 'com13' experiment to analyze and plot taxonomic
% diversity metrics, containing the following sections:
%
% 1. Define experimental data to load
%   - User inputs .xlsx file with processed 16s read data, an .xlsx file
%   mapping experimental coordinates to nutrient conditions, and an .xlsx
%   file containing OD600 growth data
%
% 2. Load and process 16s data
%   - Maps experimental coordinates in 16s .xlsx file to appropriate
%   nutrient conditions and processes organism names for visualization
%
% 3. Load and process OD data
%   - Maps experimental coordinates in OD600 .xlsx file to appropriate
%   nutrient conditions and determines outlying data points using MAD (see
%   Methods). Outputs OD600-measured community yields to determine the
%   degree to which each nutrient condition resulted in growth.
%
% 4. Calculate and plot mean and standard error of taxonomic classes by nutrient types
%   - Takes in a list of nutrients and determines which nutrient type
%   (amino acid, carbohydrate, or organic acid) each belongs to. Plots
%   organism relative abundances according to nutrient types
%
% 5. Calculate and plot coefficients of variation between replicates
%   - Calculates condition-specific coefficients of variation (CoVs) of
%   taxonomic distributions as a metric for determinism in community
%   assembly. Plots distribution of CoVs as well as CoV values according to
%   number of nutrients.
%
% 6. Plot taxonomic distributions
%   - Plots taxonomic distributions of organisms according to nutrient
%   condition (all replicates and averages).
%
% 7. Calculate and plot species richness and Shannon entropy
%
% 8. Manually calculate significance in trends for species richness or Shannon
%   - For quantifying increases in species richness or Shannon entropy with
%   increasing environmental complexity. User manually inputs the number of
%   carbon sources to compare.
%
% 9. Analyze and plot species richness and Shannon entropy comparison to CRM
%   - Takes in .mat data for consumer resource model (CRM) simulations and
%   compares simulated distributions of species richness and Shannon
%   entropy to experimental distributions.
%
% 10. Analyze and plot species richness and Shannon entropy epistasis
%   - Takes in same .mat data for CRM simulations and compares epistasis
%   metrics (see Methods) to those of experiments
%
% 11. Visualize epistasis scenario types
%
% Alan R. Pacheco, 05/19/20, modified 09/03/2020

%% 1. Define experimental data to load

clearvars
close all

taxDataFile = 'dataExperimental/level-7_Taxonomy.xlsx'; % Taxonomic distribution data
sheetName = 'level-7';
speciesResolution = 1;

mapDataFile = 'dataExperimental/coordNutrientMap.xlsx'; % Coordinate-nutrient mapping data
yieldDataFile = 'dataExperimental/com3_4_13_OD_288h.xlsx'; % Total yield data
normalized = 1;
groThreshold = 0.01;

% Map coordinates to nutrient sources
[~, ~, raw] = xlsread(mapDataFile,'Sheet1');
map.coords = raw(1:end,1);
map.nutrients = raw(1:end,2);

%% 2. Load and process 16s data
[~, ~, raw] = xlsread(taxDataFile,sheetName);
coordsAll = raw(2:end,1);
colNames = raw(1,2:end);
seqData = raw(2:end,2:end);
seqData = reshape([seqData{:}],size(seqData));
dataNorm = seqData./sum(seqData,2);

% Format taxa names
taxaNames = cell(1,length(colNames));

for i = 1:length(colNames)
    s1 = split(split(colNames(i),';'));
    if length(strfind(colNames{i},'Unassigned')) > 0
        taxaNames{i} = char(s1(1));
    else
        if speciesResolution
            s2 = split(split(s1(end-1),'__'));
            genus = s2{end};
            s3 = split(split(s1(end),'__'));
            taxaNames{i} = [char([genus(1),'. ']) char(s3(end))];
        else
            s2 = split(split(s1(end),'__'));
            taxaNames{i} = char(s2(end));
        end
    end
end

% Format coordinates
coords = coordsAll;

clearvars -except map coords seqData taxaNames yieldDataFile normalized sheetName

%% 3. Load and process OD data

sheetNames = {'com13-A','com13-B'};

for s = 1:length(sheetNames)
    [~, ~, raw] = xlsread(yieldDataFile,sheetNames{s});
    cells = raw(size(raw,1)-7:size(raw,1),2);
    dataTemp = raw(size(raw,1)-7:size(raw,1),3:14);
    dataTemp = reshape([dataTemp{:}],size(dataTemp));
    dataTemp(find(dataTemp) == 0) = NaN;
    
    sp = split(sheetNames{s},'-');

    if contains(sheetNames{s},'A')
        yieldData.(char(sp(1))) = dataTemp;
    else
        yieldData.(char(sp(1))) = [yieldData.(char(sp(1)));dataTemp];
    end

end

conditions1 = {'Glucose','Glutamate','Alanine','Pyruvate','Proline','Trehalose','Threonine','Mannitol','Fructose','GlcNAc','Sucrose','Glycerol','Gluconate','Propionate','Mannose','Succinate','Maltose','Malate','Ribose','Lactate','Cellobiose','Acetate','Galacturonate','Citrate','Arabinose','Xylose','DSerine','Galactose','G6P','Lactose','Formate','Sorbitol'};
[conditions2,conditions4,conditions8,conditions16] = deal({});

for i = 1:1:length(conditions1)/2
    conditions2 = vertcat(conditions2,strjoin([conditions1(i),conditions1(end-i+1)],'-'));
end
for i = 1:2:length(conditions2)
    conditions4 = vertcat(conditions4,strjoin(conditions2(i:i+1),'-'));
end
for i = 1:2:length(conditions4)
    conditions8 = vertcat(conditions8,strjoin(conditions4(i:i+1),'-'));
end
for i = 1:2:length(conditions8)
    conditions16 = vertcat(conditions16,strjoin(conditions8(i:i+1),'-'));
end
conditions32 = strjoin(conditions1(:),'-');

conditionsUnique = ['Neg',conditions1,conditions2',conditions4',conditions8',conditions16',conditions32]';
conditions = repelem(conditionsUnique,3);

coordsRows = [16;16;16;...
    ones(12,1)*1;ones(12,1)*2;ones(12,1)*3;ones(12,1)*4;ones(12,1)*5;ones(12,1)*6;ones(12,1)*7;ones(12,1)*8;...
    ones(12,1)*9;ones(12,1)*10;ones(12,1)*11;ones(12,1)*12;...
    ones(12,1)*13;ones(12,1)*14;...
    ones(12,1)*15;
    ones(6,1)*16;...
    ones(3,1)*16];
coordsCols = [10;11;12;repmat(1:12,1,15)';(1:9)'];

coords = [coordsRows,coordsCols];

numCSources = vertcat(0,ones(length(conditions1),1)*1,ones(length(conditions2),1)*2,ones(length(conditions4),1)*4,ones(length(conditions8),1)*8,ones(length(conditions16),1)*16,32);

% Compute growth
fieldNames = fieldnames(yieldData);
for s = 1:length(fieldNames)
    [minGro.(fieldNames{s}),maxGro.(fieldNames{s})] = deal(zeros(length(conditions),1));
    
    for i = 1:length(conditions)
        maxGro.(fieldNames{s})(i,:) = yieldData.(fieldNames{s})(coords(i,1),coords(i,2));
        minGro.(fieldNames{s})(i,:) = min(min(yieldData.(fieldNames{s})));
    end
    
    difference = maxGro.(fieldNames{s})-minGro.(fieldNames{s});
    diffGro.(fieldNames{s}) = reshape(difference(1:length(conditions) - mod(length(conditions), 3)), 3, [])';
end

% Remove outliers
for s = 1:length(fieldNames)
  
    outliersRemoved=[];
    for i = 1:size(diffGro.(fieldNames{s}),1)
        if median(diffGro.(fieldNames{s})(i,:)) > 0.01 && length(unique(diffGro.(fieldNames{s})(i,:))) == length(diffGro.(fieldNames{s})(i,:)) %If there is growth (doesn't eliminate outliers in no growth samples) and if there is variation (don't remove if values are identical, leads to Inf MAD)
            diffMads = 0.6745*abs(diffGro.(fieldNames{s})(i,:) - median(diffGro.(fieldNames{s})(i,:)))/mad(diffGro.(fieldNames{s})(i,:),1);
            if ~isempty(find(diffMads > 3.5))
                diffGro.(fieldNames{s})(i,find(diffMads > 3.5)) = NaN;
                outliersRemoved = [outliersRemoved,i];
            end
        end
    end
end

% Get mean of growth and eliminate taxonomic data for conditions that have no growth
for s = 1:length(fieldNames)
    meanDiffGro.(fieldNames{s}) = mean(diffGro.(fieldNames{s}),2);
end

% Average data over replicates
meanSeqData = zeros((size(seqData,1)-1)/3,size(seqData,2));
for i = 1:size(meanSeqData,1)
    meanSeqData(i,:) = mean(seqData(3*i-2:3*i,:),1);
end
meanSeqData(end+1,:) = seqData(end,:)';
meanSeqData(find(isnan(meanSeqData))) = 0;

clearvars -except map numCSources coords conditions conditionsUnique data diffGro outliersRemoved seqData taxaNames yieldDataFile normalized meanSeqData sheetName

%% 4. Calculate and plot mean and standard error of taxonomic classes by nutrient types
load dataExperimental/32CSNutrients.mat

if strcmp(sheetName,'level-5')
    taxaClasses = {'Other','Acinetobacter','Other','Other','Other','Other','Other','Pseudomonas','Other','Other','Other'};
    uniqueTaxaClasses = {'Acinetobacter','Pseudomonas','Other'};
else
    taxaClasses = taxaNames;
    uniqueTaxaClasses = taxaNames;
end

uniqueNutrientClasses = {'A','C','O','AC','AO','CC','CO','OO','ACO'};
[meanAbundancesByClass,stdAbundancesByClass] = deal(zeros(length(uniqueNutrientClasses),length(uniqueTaxaClasses)));
for c = 1:length(uniqueNutrientClasses)
    indicesNC = find(ismember(nutrientClasses,uniqueNutrientClasses{c}));
    
    for s = 1:length(uniqueTaxaClasses)
        indicesTC = find(ismember(taxaClasses,uniqueTaxaClasses{s}));
        meanAbundancesByClass(c,s) = mean(sum(meanSeqData(indicesNC,indicesTC),2)./sum(meanSeqData(indicesNC,:),2));
        stdAbundancesByClass(c,s) = std(sum(meanSeqData(indicesNC,indicesTC),2)./sum(meanSeqData(indicesNC,:),2))/sqrt(length(indicesNC));
    end
end

figure
if strcmp(sheetName,'level-5')
colors = [31,82,211;255,103,0;152,151,151]/255;
figure
b = bar(meanAbundancesByClass);
for k = 1:size(meanAbundancesByClass,2)
    b(k).FaceColor = colors(k,:);
end
else
bar(meanAbundancesByClass);
end
hold on
ngroups = length(uniqueNutrientClasses);
nbars = length(uniqueTaxaClasses);
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, meanAbundancesByClass(:,i), stdAbundancesByClass(:,i), '.','LineWidth',1);
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
end
hold off
set(gca,'FontSize',16,'XTickLabel',uniqueNutrientClasses)
legend(uniqueTaxaClasses)
xlabel('Nutrient category')
ylabel('Relative abundance')
title('Taxonomy by nutrient type')

%%  5. Calculate and plot coefficients of variation between replicates
relAbus = seqData./sum(seqData,2);
relAbus(find(isnan(relAbus))) = 0;
CoVs = zeros(floor(size(relAbus,1)/3),1);
count = 1;
for i = 1:length(CoVs)
    
    covSamp = zeros(1,size(relAbus,2)); % records the coefficients of variation for each organism in the sample
    for j = 1:size(relAbus,2)
        sample = relAbus([i*3-2:i*3],j);
        covSamp(j) = std(sample)/mean(sample);
    end
    
    covSamp(find(sum(relAbus([i*3-2:i*3],:),1) == 0)) = []; % Remove organisms that had zero abundance, so as to not skew the metric

    CoVs(count) = mean(covSamp);
    count = count + 1;
end

figure
subplot(1,2,1)
histogram(CoVs,'BinWidth',0.2)
set(gca,'FontSize',16)
xlabel('Average coefficient of variation')
ylabel('Frequency')
title('Average coefficients of variation')

subplot(1,2,2)
scatter(numCSources(2:end),CoVs,40,'Filled')
set(gca,'FontSize',16)
ylabel('Average coefficient of variation')
xlabel('Number of carbon sources')
title('Average CoV by NCS')

[spearmanCoV,pValCoV] = corr(numCSources(2:end),CoVs,'Type','Spearman'); % Calculate Spearman correlation coefficient of coefficient of variation to number of carbon sources

%% 6. Plot taxonomic distributions

% Plot relative abundances by nutrient source, all replicates
figure
bar(seqData./sum(seqData,2),'stacked')
set(gca,'FontSize',14)
set(gca,'XTick',[2:3:length(map.nutrients),length(map.nutrients)]);
set(gca,'XTickLabel',[map.nutrients(2:3:length(map.nutrients));'';'Initial']);
set(gca,'XTickLabelRotation',45)
ylim([0,1])
xlim([0,193])
legend(taxaNames)
title('Taxonomic distribution, all replicates')
xlabel('Nutrient condition')
ylabel('Relative abundance')

% Plot averaged data
figure
bar(meanSeqData(1:63,:)./sum(meanSeqData(1:63,:),2),'stacked')
set(gca,'FontSize',14)
set(gca,'XTick',1:1:length(meanSeqData)-1);
set(gca,'XTickLabel',[map.nutrients(2:3:length(map.nutrients))]);
set(gca,'XTickLabelRotation',45)
ylim([0,1])
xlim([0,64])
legend(taxaNames)
title('Mean taxonomic distribution')
xlabel('Nutrient condition')
ylabel('Relative abundance')

%% 7. Calculate and plot species richness and Shannon entropy

% Plot alpha diversity by number of nutrients
numCSourcesAllReplicates = repelem(numCSources(2:end),3);
speciesRichness = zeros(size(meanSeqData,1)-1,1);
for i = 1:length(speciesRichness)
    speciesRichness(i) = length(find(meanSeqData(i,:)));
end

figure
subplot(1,2,1)
boxplot(speciesRichness,numCSources(2:end),'Notch','off')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'FontSize',16)
xlabel('Number of nutrients')
ylabel('Species richness')
title('Species richness by number of nutrients')

% Plot Shannon entropy by number of nutrients
shannon = zeros(size(meanSeqData,1)-1,1);

for i = 1:length(shannon)
    s = sum(meanSeqData(i,:));
    
    for j = 1:size(relAbus,2)
        if meanSeqData(i,j) > 0
            shannon(i) = shannon(i) - (meanSeqData(i,j)/s)*log2(meanSeqData(i,j)/s);
        end
    end
end

subplot(1,2,2)
boxplot(shannon,numCSources(2:end),'Notch','off')
set(findobj(gca,'type','line'),'linew',2)
set(gca,'FontSize',16)
xlabel('Number of nutrients')
ylabel('Shannon entropy')
title('Shannon entropy by number of nutrients')

%% 8. Manually calculate significance in trends for species richness or Shannon
numCarbonSourcesToCompare = [16,32]; % Manually input number of carbon sources to compare

sample1 = speciesRichness(find(numCSources == numCarbonSourcesToCompare(1))-1);
sample2 = speciesRichness(find(numCSources == numCarbonSourcesToCompare(2))-1);
[hSpeciesRichness,pSpeciesRichness] = ttest2(sample1,sample2)

sample1 = shannon(find(numCSources == numCarbonSourcesToCompare(1))-1);
sample2 = shannon(find(numCSources == numCarbonSourcesToCompare(2))-1);
[hShannon,pShannon] = ttest2(sample1,sample2)

%% 9. Analyze and plot species richness and Shannon entropy comparison to CRM
uniqueNumCSources = unique(numCSources);
lineColors = [149,178,176;51, 101, 138;0, 21, 36]./255;

dataCRM_NoGeneralists = 'dataConsumerResource/epistasisS13_N32_G0_RS50_nSM10_062320.mat';
dataCRM_WithGeneralists = 'dataConsumerResource/epistasisS13_N32_GV_RS50_nSM10_062320.mat';

figure
subplot(1,2,1)
load(dataCRM_NoGeneralists)
errorbar(epiCRM.uniqueNumNutrientsByCase(2:end),mean(epiCRM.meanGrowingSpeciesByNumNutrients(2:end,2:end),2),mean(epiCRM.stdErrorGrowingSpeciesByNumNutrients(2:end,2:end),2),mean(epiCRM.stdErrorGrowingSpeciesByNumNutrients(2:end,2:end),2),'LineWidth',4,'Color',lineColors(1,:))
hold on
load(dataCRM_WithGeneralists)
errorbar(epiCRM.uniqueNumNutrientsByCase(2:end),mean(epiCRM.meanGrowingSpeciesByNumNutrients(2:end,2:end),2),mean(epiCRM.stdErrorGrowingSpeciesByNumNutrients(2:end,2:end),2),mean(epiCRM.stdErrorGrowingSpeciesByNumNutrients(2:end,2:end),2),'LineWidth',4,'Color',lineColors(2,:))
hold on
[meanRichnessByNumCSources,stdErrorRichnessByNumCSources] = deal(zeros(length(uniqueNumCSources)-1,1));
for i = 2:length(uniqueNumCSources)
    indicesCurr = find(numCSources == uniqueNumCSources(i));
    meanRichnessByNumCSources(i-1) = mean(speciesRichness(indicesCurr-1));
    stdErrorRichnessByNumCSources(i-1) = std(speciesRichness(indicesCurr-1))/length(speciesRichness(indicesCurr-1));
end
errorbar(uniqueNumCSources(2:end),meanRichnessByNumCSources,stdErrorRichnessByNumCSources,stdErrorRichnessByNumCSources,'LineWidth',4,'Color',lineColors(3,:))
set(gca,'FontSize',16,'XTick',epiCRM.uniqueNumNutrientsByCase(2:end))
xlabel('Number of carbon sources')
ylabel('Species richness')
legend('CRM-A','CRM-B','Experiment')
title('Species richness comparison to CRM')
ylim([0.5,6.5])

subplot(1,2,2)
lineColors = [252, 159, 91;228,63,111;86,54,53]./255;
load(dataCRM_NoGeneralists)
errorbar(epiCRM.uniqueNumNutrientsByCase(2:end),mean(epiCRM.meanShannonByNumNutrients(2:end,2:end),2),mean(epiCRM.stdErrorShannonByNumNutrients(2:end,2:end),2),mean(epiCRM.stdErrorShannonByNumNutrients(2:end,2:end),2),'LineWidth',4,'Color',lineColors(1,:))
hold on
load(dataCRM_WithGeneralists)
errorbar(epiCRM.uniqueNumNutrientsByCase(2:end),mean(epiCRM.meanShannonByNumNutrients(2:end,2:end),2),mean(epiCRM.stdErrorShannonByNumNutrients(2:end,2:end),2),mean(epiCRM.stdErrorShannonByNumNutrients(2:end,2:end),2),'LineWidth',4,'Color',lineColors(2,:))
hold on
[meanShannonByNumCSources,stdErrorShannonByNumCSources] = deal(zeros(length(uniqueNumCSources)-1,1));
for i = 2:length(uniqueNumCSources)
    indicesCurr = find(numCSources == uniqueNumCSources(i));
    meanShannonByNumCSources(i-1) = mean(shannon(indicesCurr-1));
    stdErrorShannonByNumCSources(i-1) = std(shannon(indicesCurr-1)/length(shannon(indicesCurr-1)));
end
errorbar(uniqueNumCSources(2:end),meanShannonByNumCSources,stdErrorShannonByNumCSources,stdErrorShannonByNumCSources,'LineWidth',4,'Color',lineColors(3,:))
set(gca,'FontSize',16,'XTick',epiCRM.uniqueNumNutrientsByCase(2:end))
xlabel('Number of carbon sources')
ylabel('Shannon entropy')
legend('CRM-A','CRM-B','Experiment')
title('Shannon entropy comparison to CRM')
ylim([0,2.4])

%% 10. Analyze and plot species richness and Shannon entropy epistasis

% Load CRM file for comparison
load(dataCRM_NoGeneralists)

[richnessDiff,shannonDiff] = deal(zeros(length(speciesRichness),1));
count = length(find(sum(nutrientComboMap,2)==1))+1;
for i = 1:length(uniqueNumCSources)-2
    CSourceCurr = uniqueNumCSources(i+2);
    
    indicesHigherAll = find(sum(nutrientComboMap,2) == uniqueNumCSources(i+2)); % Indices of more complex combos
    indicesLowerAll = find(sum(nutrientComboMap,2) == uniqueNumCSources(i+1)); % Indices of next lowest combo
    for j = 1:length(indicesHigherAll)
        if i > 1
            indicesLower = indicesLowerAll(j*2-1:j*2); % Indices of specific lower combo
        else
            indicesLower = find(nutrientComboMap(indicesHigherAll(j),:));
        end
        
        readsHigh = meanSeqData(indicesHigherAll(j),:);
        readsLowSum = sum(meanSeqData(indicesLower,:),1); % Sum of abundances in constituent lower combos
        
        readsLowFromEach = sum(meanSeqData(indicesLower,:)~=0,2);
        
        richnessDiff(count) = length(find(readsHigh)) - max(readsLowFromEach);
        
        [shannonExpected,shannonReal] = deal(0);
        s = sum(readsHigh);
        for k = 1:length(readsHigh)
            if readsHigh(k) > 0
                shannonReal = shannonReal - (readsHigh(k)/s)*log2(readsHigh(k)/s);
            end
        end
        s = sum(readsLowSum);
        for k = 1:length(readsLowSum)
            if readsLowSum(k) > 0
                shannonExpected = shannonExpected - (readsLowSum(k)/s)*log2(readsLowSum(k)/s);
            end
        end
        
        shannonDiff(count) = shannonReal-shannonExpected;
        
        count = count + 1;
    end
end

figure
subplot(1,2,1)
histogram(richnessDiff(33:end))
xlabel('Species Richness Epistasis')
ylabel('Frequency')
set(gca,'FontSize',16)
ylim([0,12])

subplot(1,2,2)
histogram(shannonDiff(33:end))
xlabel('Shannon Entropy Epistasis')
ylabel('Frequency')
set(gca,'FontSize',16)
ylim([0,12])

lineColors = [0, 21, 36;149,178,176;51, 101, 138]./255;

figure
binWidth = 0.5;
load(dataCRM_NoGeneralists)
toplot = mean(epiCRM.alphaDiff(34:end,2:end),2);
low = floor(min(toplot))+floor((min(toplot)-floor(min(toplot)))/(binWidth/2))*(binWidth/2);
high = floor(max(toplot))+ceil((max(toplot)-floor(max(toplot)))/(binWidth/2))*(binWidth/2);
histogram(toplot,'BinEdges',[low-binWidth/2:binWidth:high+binWidth/2],'DisplayStyle','stairs','LineWidth',4,'EdgeColor',lineColors(2,:))
hold on
load(dataCRM_WithGeneralists)
toplot = mean(epiCRM.alphaDiff(34:end,2:end),2);
low = floor(min(toplot))+floor((min(toplot)-floor(min(toplot)))/(binWidth/2))*(binWidth/2);
high = floor(max(toplot))+ceil((max(toplot)-floor(max(toplot)))/(binWidth/2))*(binWidth/2);
histogram(toplot,'BinEdges',[low-binWidth/2:binWidth:high+binWidth/2],'DisplayStyle','stairs','LineWidth',4,'EdgeColor',lineColors(3,:))
hold on
binWidth = 1;
toplot = richnessDiff(33:end);
low = floor(min(toplot))+floor((min(toplot)-floor(min(toplot)))/(binWidth/2))*(binWidth/2);
high = floor(max(toplot))+ceil((max(toplot)-floor(max(toplot)))/(binWidth/2))*(binWidth/2);
histogram(toplot,'BinEdges',[low-binWidth/2:binWidth:high+binWidth/2],'DisplayStyle','stairs','LineWidth',4,'EdgeColor',lineColors(1,:))

legend({'CRM-A','CRM-B','Experiment'})
set(gca,'FontSize',16)
xlim([-6,6])
ylim([0,40])
title('Species richness epistasis comparison')

% Get statistics on mean,variance of species richness epistasis
% XLimits = get(gca,'XLim');
% figure
% alphaDistSpan = [mean(richnessDiff(33:end))-std(richnessDiff(33:end)):1e-3:mean(richnessDiff(33:end))+std(richnessDiff(33:end))];
% plot(alphaDistSpan,ones(length(alphaDistSpan),1))
% hold on
% scatter(mean(richnessDiff(33:end)),1)
% load(dataCRM_NoGeneralists)
% CRMDistSpan = [mean(mean(epiCRM.alphaDiff(34:end,2:end),2))-std(mean(epiCRM.alphaDiff(34:end,2:end),2)):1e-3:mean(mean(epiCRM.alphaDiff(34:end,2:end),2))+std(mean(epiCRM.alphaDiff(34:end,2:end),2))];
% plot(CRMDistSpan,2*ones(length(CRMDistSpan),1))
% hold on
% scatter(mean(mean(epiCRM.alphaDiff(34:end,2:end),2)),2)
% hold on
% load(dataCRM_WithGeneralists)
% CRMDistSpan = [mean(mean(epiCRM.alphaDiff(34:end,2:end),2))-std(mean(epiCRM.alphaDiff(34:end,2:end),2)):1e-3:mean(mean(epiCRM.alphaDiff(34:end,2:end),2))+std(mean(epiCRM.alphaDiff(34:end,2:end),2))];
% plot(CRMDistSpan,2*ones(length(CRMDistSpan),1))
% hold on
% scatter(mean(mean(epiCRM.alphaDiff(34:end,2:end),2)),2)
% ylim([0,3])
% xlim(XLimits)

% Plot shannon epistasis
lineColors = [86,54,53;252, 159, 91;228,63,111]./255;
figure
binWidth = 0.25;
load(dataCRM_NoGeneralists)
toplot = mean(epiCRM.shannonDiff(34:end,2:end),2);
low = floor(min(toplot))+floor((min(toplot)-floor(min(toplot)))/(binWidth/2))*(binWidth/2);
high = floor(max(toplot))+ceil((max(toplot)-floor(max(toplot)))/(binWidth/2))*(binWidth/2);
histogram(toplot,'BinEdges',[low:binWidth:high],'DisplayStyle','stairs','LineWidth',4,'EdgeColor',lineColors(2,:))
hold on
load(dataCRM_WithGeneralists)
toplot = mean(epiCRM.shannonDiff(34:end,2:end),2);
low = floor(min(toplot))+floor((min(toplot)-floor(min(toplot)))/(binWidth/2))*(binWidth/2);
high = floor(max(toplot))+ceil((max(toplot)-floor(max(toplot)))/(binWidth/2))*(binWidth/2);
histogram(toplot,'BinEdges',[low-binWidth/2:binWidth:high+binWidth/2],'DisplayStyle','stairs','LineWidth',4,'EdgeColor',lineColors(3,:))
hold on
binWidth = 0.5;
toplot = shannonDiff(33:end);
low = floor(min(toplot))+floor((min(toplot)-floor(min(toplot)))/(binWidth/2))*(binWidth/2);
high = floor(max(toplot))+ceil((max(toplot)-floor(max(toplot)))/(binWidth/2))*(binWidth/2);
histogram(toplot,'BinEdges',[low:binWidth:high],'DisplayStyle','stairs','LineWidth',4,'EdgeColor',lineColors(1,:))

legend({'CRM-A','CRM-B','Experiment'})
set(gca,'FontSize',16)
xlim([-3,3])
ylim([0,40])
title('Shannon entropy epistasis comparison')

% Get statistics on mean, variance of shannon epistasis
% XLimits = get(gca,'XLim');
% figure
% shannonDistSpan = [mean(shannonDiff(33:end))-std(shannonDiff(33:end)):1e-3:mean(shannonDiff(33:end))+std(shannonDiff(33:end))];
% plot(shannonDistSpan,ones(length(shannonDistSpan),1))
% hold on
% scatter(mean(shannonDiff(33:end)),1)
% load(dataCRM_NoGeneralists)
% CRMDistSpan = [mean(mean(epiCRM.shannonDiff(34:end,2:end),2))-std(mean(epiCRM.shannonDiff(34:end,2:end),2)):1e-3:mean(mean(epiCRM.shannonDiff(34:end,2:end),2))+std(mean(epiCRM.shannonDiff(34:end,2:end),2))];
% plot(CRMDistSpan,2*ones(length(CRMDistSpan),1))
% hold on
% scatter(mean(mean(epiCRM.shannonDiff(34:end,2:end),2)),2)
% hold on
% load(dataCRM_WithGeneralists)
% CRMDistSpan = [mean(mean(epiCRM.shannonDiff(34:end,2:end),2))-std(mean(epiCRM.shannonDiff(34:end,2:end),2)):1e-3:mean(mean(epiCRM.shannonDiff(34:end,2:end),2))+std(mean(epiCRM.shannonDiff(34:end,2:end),2))];
% plot(CRMDistSpan,2*ones(length(CRMDistSpan),1))
% hold on
% scatter(mean(mean(epiCRM.shannonDiff(34:end,2:end),2)),2)
% ylim([0,3])
% xlim(XLimits)

%% 11. Visualize epistasis scenario types
									
biologProportions = [0.708333333333333,0.885416666666667,0.427083333333333,0.645833333333333,0.770833333333333,0.416666666666667,0.906250000000000,0.833333333333333,0.895833333333333,0.458333333333333,0.822916666666667,0.468750000000000,0.697916666666667]; % Proportions of BIOLOG carbon sources on which growth was observed
meanRelAbus = meanSeqData(1:63,:)./sum(meanSeqData(1:63,:),2);

type1Indices = [33,38,39,42,50,53,60];
type2Indices = [54,57,61,62,63];
type3Indices = [35,36,37,41,43,44,46,48,49,52,55,56,58];
type4Indices = [34,40,45,47,51,59];

typeLengths = [length(type1Indices),length(type2Indices),length(type3Indices),length(type4Indices)];

figure
bar(typeLengths/sum(typeLengths),'FaceColor',[101, 82, 77]/255)
xlim([0.5,4.5])
set(gca,'FontSize',16)
ylabel('Relative frequency')
xlabel('Epistasis type')
title('Epistasis types')