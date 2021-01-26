% analyzeYield
%
% Processes OD600 growth data from 'com3', 'com4', and 'com13' experiments
% to analyze and plot biomass yields according to environmental complexity,
% containing the following sections:
%
% 1. Define experimental data to load
%   - User inputs .xlsx file with processed OD600 growth data and an .xlsx file
%   mapping experimental coordinates to nutrient conditions.
%
% 2. Load and process OD data
%   - Maps experimental coordinates in OD600 .xlsx file to appropriate
%   nutrient conditions and determines outlying data points using MAD (see
%   Methods). Outputs OD600-measured community yields to determine the
%   degree to which each nutrient condition resulted in growth.
%
% 3. Calculate and plot mean and standard error of taxonomic classes by nutrient types
%   - Takes in a list of nutrients and determines which nutrient type
%   (amino acid, carbohydrate, or organic acid) each belongs to. Plots
%   organism relative abundances according to nutrient types
%
% 4. Manually calculate significance in trends for growth yields
%   - For quantifying changes in yield with increasing environmental
%   complexity. User manually inputs the number of carbon sources to
%   compare.
%
% 5. Plot yields in all nutrient conditions
%
% 6. Plot epistasis distributions
%   - Takes in .mat data for consumer resource model (CRM) simulations and
%   compares simulated distributions of yields to experimental
%   distributions.
%
% Alan R. Pacheco, 05/29/20, modified 09/03/2020

%% 1. Define experimental data to load

clearvars
close all

expDataFile = 'dataExperimental/com3_4_13_OD_288h.xlsx';
sheetNames = {'com3-A','com4-A','com3-B','com4-B','com13-A','com13-B'};
groThreshold = 0.01;
colors = [255,188,66;4,150,255;57,62,65]./255;

%% 2. Load and process OD data

for s = 1:length(sheetNames)
    [~, ~, raw] = xlsread(expDataFile,sheetNames{s});
    cells = raw(size(raw,1)-7:size(raw,1),2);
    dataTemp = raw(size(raw,1)-7:size(raw,1),3:14);
    dataTemp = reshape([dataTemp{:}],size(dataTemp));
    dataTemp(find(dataTemp) == 0) = NaN;
    
    sp = split(sheetNames{s},'-');

    if contains(sheetNames{s},'A')
        data.(char(sp(1))) = dataTemp;
    else
        data.(char(sp(1))) = [data.(char(sp(1)));dataTemp];
    end

end

clearvars -except data groThreshold sheetNames colors

% Map conditions to coordinates
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

clearvars -except numCSources coords conditions conditionsUnique data groThreshold sheetNames colors

% Compute growth
fieldNames = fieldnames(data);
for s = 1:length(fieldNames)
    [minGro.(fieldNames{s}),maxGro.(fieldNames{s})] = deal(zeros(length(conditions),1));
    
    for i = 1:length(conditions)
        maxGro.(fieldNames{s})(i,:) = data.(fieldNames{s})(coords(i,1),coords(i,2));
        minGro.(fieldNames{s})(i,:) = min(min(data.(fieldNames{s})));
    end
    
    difference = maxGro.(fieldNames{s})-minGro.(fieldNames{s});
    diffGro.(fieldNames{s}) = reshape(difference(1:length(conditions) - mod(length(conditions), 3)), 3, [])';
end

clearvars -except numCSources coords conditions conditionsUnique data maxGro diffGro fieldNames groThreshold sheetNames colors

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

for s = 1:length(fieldNames)

    meanDiffGro.(fieldNames{s}) = nanmean(diffGro.(fieldNames{s}),2);
    stdDiffGro.(fieldNames{s}) = nanstd(diffGro.(fieldNames{s})')';
end

clearvars -except numCSources coords conditions conditionsUnique data maxGro meanDiffGro diffGro fieldNames groThreshold colors

%% 3. Plot OD600 by number of carbon sources

uniqueCSourceNums = unique(numCSources);

figure
for s = 1:length(fieldNames)
    subplot(1,length(fieldNames),s)
    boxplot(meanDiffGro.(fieldNames{s}),numCSources)
    set(findobj(gca,'type','line'),'linew',2)
    set(gca,'FontSize',16)
    xlabel('Number of carbon sources')
    ylabel('Biomass (OD600)')
    title([fieldNames{s}])
    ylim([0,0.55])  

end

%% 4. Manually calculate significance in trends for growth yields

s = 2;
numCarbonSourcesToCompare = [1,32];

sample1 = meanDiffGro.(fieldNames{s})(find(numCSources == numCarbonSourcesToCompare(1)));
sample2 = meanDiffGro.(fieldNames{s})(find(numCSources == numCarbonSourcesToCompare(2)));

[h,p] = ttest2(sample1,sample2)

%% 5. Plot yields in all nutrient conditions

colors = [255,188,66;4,150,255;57,62,65]./255; % For communities
% colors = [48,44,129;66,154,203;255,184,111;170,188,116]./255; % For individual organisms

[~, ~, raw] = xlsread('dataExperimental/coordNutrientMap.xlsx','Sheet1');
map.coords = raw(1:end,1);
map.nutrients = raw(1:end,2);

ctrs = [1:length(conditionsUnique)];
[plotData,b] = deal([]);
for s = 1:length(fieldNames)
    plotData = [plotData,nanmean(diffGro.(fieldNames{s}),2)];
    b = [b,nanstd(diffGro.(fieldNames{s})')'];
end

figure
hBar = bar(ctrs, plotData);
colormap(colors)
for k1 = 1:size(plotData,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, [hBar(k1).XOffset]');
    ydt(k1,:) = hBar(k1).YData;
end
hold on
errorbar(ctr, ydt, b', '.k')
hold off
legend(fieldNames)
set(gca,'FontSize',16,'XTick',[1:length(conditionsUnique)],'XTickLabel',['Neg';map.nutrients(2:3:length(map.nutrients)-1)],'XTickLabelRotation',45)
xlabel('Nutrients')
ylabel('Biomass (OD600)')
xlim([0,length(conditionsUnique)+1])
ylim([0,.55])

%% 6. Plot epistasis distributions

% Load CRM data for comparison
load dataConsumerResource/epistasisS13_N32_G0_RS50_nSM10.mat
load dataExperimental/32CSNutrients.mat

yieldDiff = struct();
for s = 1:length(fieldNames)
    yieldDiff.(fieldNames{s}) = zeros(size(nutrientComboMap,1),1);
    count = length(find(sum(nutrientComboMap,2)==1))+2;
    
    for i = 1:length(uniqueCSourceNums)-2
        CSourceCurr = uniqueCSourceNums(i+2);

        indicesHigherAll = find(sum(nutrientComboMap,2) == uniqueCSourceNums(i+2)); % Indices of more complex combos
        indicesLowerAll = find(sum(nutrientComboMap,2) == uniqueCSourceNums(i+1)); % Indices of next lowest combo
        for j = 1:length(indicesHigherAll)
            if i > 1
                indicesLower = indicesLowerAll(j*2-1:j*2); % Indices of specific lower combo
            else
                indicesLower = find(nutrientComboMap(indicesHigherAll(j),:));
            end

            yieldsHigh = meanDiffGro.(fieldNames{s})(indicesHigherAll(j)+1,:);
            yieldsLow = sum(meanDiffGro.(fieldNames{s})(indicesLower+1,:),1); % Sum of abundances in constituent lower combos

            yieldDiff.(fieldNames{s})(count) = yieldsHigh - yieldsLow/2;

            count = count + 1;
        end
    end
end

% Plot individual distributions of yield epistasis
colors = [255,188,66;4,150,255;57,62,65]./255; % For communities
figure
binWidth = 0.025;
toplot = mean(epiCRM.yieldDiff(34:end,2:end),2);
low = floor(min(toplot))+floor((min(toplot)-floor(min(toplot)))/(binWidth/2))*(binWidth/2);
high = floor(max(toplot))+ceil((max(toplot)-floor(max(toplot)))/(binWidth/2))*(binWidth/2);
histogram(toplot,'BinEdges',[low:binWidth:high+binWidth],'DisplayStyle','stairs','LineWidth',2,'EdgeColor',[229, 190, 237]./255)
hold on
for s = 1:length(fieldNames)
    toplot = yieldDiff.(fieldNames{s})(34:end);
    low = floor(min(toplot))+floor((min(toplot)-floor(min(toplot)))/(binWidth/2))*(binWidth/2);
    high = floor(max(toplot))+ceil((max(toplot)-floor(max(toplot)))/(binWidth/2))*(binWidth/2);
    if strcmp(fieldNames{s},'com3')
        histogram(toplot,'BinEdges',[low:binWidth:high+binWidth],'DisplayStyle','stairs','LineWidth',2,'EdgeColor',colors(s,:))
    else
        histogram(toplot,'BinEdges',[low-binWidth/2:binWidth:high+binWidth/2],'DisplayStyle','stairs','LineWidth',2,'EdgeColor',colors(s,:))
    end
    hold on
end
xlabel('Yield Epistasis')
ylabel('Frequency')
set(gca,'FontSize',16)
legend(['CRM';fieldNames])
ylim([0,40])
XLimits = get(gca,'XLim');

% Get statistics on mean, variance of epistasis
% figure
% com3DistSpan = [mean(yieldDiff.com3(34:end))-std(yieldDiff.com3(34:end)):1e-3:mean(yieldDiff.com3(34:end))+std(yieldDiff.com3(34:end))];
% plot(com3DistSpan,ones(length(com3DistSpan),1))
% hold on
% scatter(mean(yieldDiff.com3(34:end)),1)
% hold on
% com4DistSpan = [mean(yieldDiff.com4(34:end))-std(yieldDiff.com4(34:end)):1e-3:mean(yieldDiff.com4(34:end))+std(yieldDiff.com4(34:end))];
% plot(com4DistSpan,2*ones(length(com4DistSpan),1))
% hold on
% scatter(mean(yieldDiff.com4(34:end)),2)
% hold on
% com13DistSpan = [mean(yieldDiff.com13(34:end))-std(yieldDiff.com13(34:end)):1e-3:mean(yieldDiff.com13(34:end))+std(yieldDiff.com13(34:end))];
% plot(com13DistSpan,3*ones(length(com13DistSpan),1))
% hold on
% scatter(mean(yieldDiff.com13(34:end)),3)
% hold on
% CRMDistSpan = [mean(mean(epiCRM.yieldDiff(34:end,2:end),2))-std(mean(epiCRM.yieldDiff(34:end,2:end),2)):1e-3:mean(mean(epiCRM.yieldDiff(34:end,2:end),2))+std(mean(epiCRM.yieldDiff(34:end,2:end),2))];
% plot(CRMDistSpan,4*ones(length(CRMDistSpan),1))
% hold on
% scatter(mean(mean(epiCRM.yieldDiff(34:end,2:end),2)),4)
% xlim(XLimits)
% ylim([0,5])