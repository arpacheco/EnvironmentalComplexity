% solveCRMGeneric.m 
%
% Runs statistical ensembles of consumer resource models from CRM.m based
% on species-specific nutrient utilization probabilities
%
% Alan R. Pacheco, 05/07/2020, modified 09/03/2020

%% Define filenames, species, and nutrients

clearvars
close all

numSpecies = 13; % Number of organisms

% Uncomment below for all-specialist population (CRM-A)
isGeneralist = zeros(1,numSpecies);
specialistConsumptionProb = 0.5*ones(1,numSpecies);

% Uncomment below for experimentally-derived species-specific nutrient utilization probabilities (CRM-B)
% load biologGrowthProbsJun2019.mat
% specialistConsumptionProb = biologGrowthProbs/100;
% isGeneralist = zeros(1,numSpecies);

% Modeling parameters
maxSecMets = 10; % Max number of available secreted metabolites to test
transfProb = 0.5; % Probability of each nutrient being secreted as another
numRandSims = 50; % Number of random simulations per number of secreted metabolites
monod = 1; % Use Monod dynamics, default = 1
sequentialConsumption = 0; % Have the organisms exhaust one nutrient before moving to the next, default = 0

%% Run CRM
load binaryNutrients32.mat

[growingSpecies,yields,shannon] = deal(zeros(size(nutrientMat,1),maxSecMets+1,numRandSims));
relAbus = zeros(size(nutrientMat,1),numSpecies,maxSecMets+1,numRandSims);
for sm = 1:maxSecMets + 1

    numSecMets = sm - 1
    
    for q = 2:size(nutrientMat,1)

        selectNutrientCombo = q;

        relAbusPerCombo = zeros(numRandSims,numSpecies);
        [yieldsPerCombo,growingSpeciesPerCombo,shannonPerCombo] = deal(zeros(numRandSims,1));

        for r = 1:numRandSims
                        
            % Get nutrient concentrations
            selectNutrients = nutrients(find(nutrientMat(selectNutrientCombo,:)));

            selectNutrientsIndices = find(ismember(nutrients,selectNutrients));


            % Set simulation parameters

            global S M C D g w l m k d kR monodTerm sequential

            % N is CFU/mL, R is g/mL

            S = numSpecies;
            g = ones(1,S).*1; % Conversion factor from energy uptake to growth rate (1/energy) SET TO 1
            m = ones(1,S).*.05; % Minimal energy uptake for maintenance of each species (energy/hour)

            M = length(selectNutrients); % Number of limiting nutrients

            monodTerm = monod;
            sequential = sequentialConsumption;
            
            % Use the D matrix to get all the metabolites in the simulation
            if numSecMets > 0

                alphabet = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};
                transformedNutrients = alphabet(1:maxSecMets);

                D = zeros(length(selectNutrients)+length(transformedNutrients),length(selectNutrients)+length(transformedNutrients),S);
                for k = 1:S
                    DCurr = D(:,:,k);
                    for i = length(selectNutrients)+1:length(DCurr)
                        for j = 1:length(DCurr)
                            if rand(1) < transfProb && i ~= j
                                DCurr(i,j) = 1;
                            end
                        end
                    end
                    DCurr = DCurr./sum(DCurr,1); % Rows must add to 1
                    D(:,:,k) = DCurr;
                end
                D(find(isnan(D))) = 0;

                simulationNutrients = [string(num2cell(selectNutrients)),transformedNutrients];

            else
                simulationNutrients = string(num2cell(selectNutrients));
                D = zeros(length(simulationNutrients),length(simulationNutrients),S);
            end

            w = ones(1,length(simulationNutrients))*1e9; % Energy content of resource a (energy/g)
            l = ones(1,length(simulationNutrients)).*0.25; % Leakage fraction for resource a (unitless) SET TO 0.8

            % Create C matrix based on generalist/specialist designation
            C = zeros(S,length(simulationNutrients)); % Uptake rate per unit concentration of resource a by species i (mL/hour)
            for s = 1:S
                randPrefs = rand(1,length(simulationNutrients));
                if isGeneralist(s)
                    C(s,find(randPrefs >= 1-generalistConsumptionProb(s))) = (1);
                else
                    C(s,find(randPrefs >= 1-specialistConsumptionProb(s))) = rand(1);
                end
            end

            C = C.*1e-5;

            d = 1/(10/(300*48)); % Timescale for resource dilution 

            % Set initial conditions and solve system of ODEs

            N0 = ones(1,S).*6e6; % approximately OD 0.5 diluted to 5µl/300µl

            R0 = ones(1,length(simulationNutrients)).*1e-9; % Baseline must be greater than 0 to allow for D matrix turnover
            % R0(find(ismember(simulationNutrients,selectNutrients))) = nutrientMassPerVolume(selectNutrientsIndices)'./length(selectNutrients); % g/L, scaled by number of nutrients
            totalNutrient = 1.5;
            R0(find(ismember(simulationNutrients,string(num2cell(selectNutrients))))) = totalNutrient./length(selectNutrients);  % g/L, scaled but otherwise equal since CRM does not care about number of carbon atoms
            k = R0/48; % External supply of resource a (grams/mL/hour), equals R0 divided by 48 hour replenishment period
            kR = ones(1,length(simulationNutrients)).*250;

            I0 = [N0 R0];

            timeRange = [0:.1:288];
            [T,X]=ode15s('CRM',timeRange,I0);

            % Get relative abundances

            N = X(:,1:S);
            R = X(:,S+1:end);

            abundances = N(end,:);
            abundances(find(abundances < N0)) = 0;

            relAbu = abundances./sum(abundances);

            relAbusPerCombo(r,:) = relAbu;
            yieldsPerCombo(r) = sum(abundances)/8e8;
            growingSpeciesPerCombo(r) = length(find(abundances > N0));
            for a = 1:length(relAbu)
                if relAbu(a) > 0
                    shannonPerCombo(r) = shannonPerCombo(r) - relAbu(a).*log2(relAbu(a));
                end
            end
        end
        
        relAbus(q,:,sm,:) = relAbusPerCombo';
        yields(q,sm,:) = yieldsPerCombo;
        growingSpecies(q,sm,:) = growingSpeciesPerCombo;  
        shannon(q,sm,:) = shannonPerCombo;
    end

relAbus(find(isnan(relAbus))) = 0;

%% Plot overall phenotypes
close all

numNutrientsByCase = sum(nutrientMat,2);
uniqueNumNutrientsByCase = unique(numNutrientsByCase);

figure

% Plot species richness
[meanGrowingSpeciesByNumNutrients,stdErrorGrowingSpeciesByNumNutrients] = deal(zeros(length(uniqueNumNutrientsByCase),maxSecMets+1));
for i = 1:length(uniqueNumNutrientsByCase)
    meanGrowingSpeciesByNumNutrients(i,:) = mean(mean(growingSpecies(find(ismember(numNutrientsByCase,uniqueNumNutrientsByCase(i))),:,:),3),1);
    if i > 1
        stdErrorGrowingSpeciesByNumNutrients(i,:) = std(squeeze(mean(squeeze(growingSpecies(find(ismember(numNutrientsByCase,uniqueNumNutrientsByCase(i))),:,:)),1))');
    else
        stdErrorGrowingSpeciesByNumNutrients(i,:) = mean(mean(growingSpecies(find(ismember(numNutrientsByCase,uniqueNumNutrientsByCase(i))),:,:),3),1);
    end
end
stdErrorGrowingSpeciesByNumNutrients = stdErrorGrowingSpeciesByNumNutrients./sqrt(numRandSims);
imagesc(flipud(meanGrowingSpeciesByNumNutrients))
colormap(summer)
xlabel('Number of secreted metabolites')
ylabel('Number of nutrients')
set(gca,'FontSize',16,'YTick',[1:length(uniqueNumNutrientsByCase)],'YTickLabel',flipud(uniqueNumNutrientsByCase),'XTick',[1:maxSecMets+1],'XTickLabel',[0:maxSecMets])
title(['Species richness; S = ' num2str(numSpecies) ', generalists = ' num2str(length(find(isGeneralist))) ', tP = ' num2str(transfProb) ', randSims = ' num2str(numRandSims)])
colorbar
caxis([0 7])

% Plot shannon entropy
[meanShannonByNumNutrients,stdErrorShannonByNumNutrients] = deal(zeros(length(uniqueNumNutrientsByCase),maxSecMets+1));
maxShannon = numSpecies*(1/numSpecies)*log2(1/numSpecies)*-1;
for i = 1:length(uniqueNumNutrientsByCase)
    meanShannonByNumNutrients(i,:) = mean(mean(shannon(find(ismember(numNutrientsByCase,uniqueNumNutrientsByCase(i))),:,:),3),1);
    if i > 1
        stdErrorShannonByNumNutrients(i,:) = std(squeeze(mean(squeeze(shannon(find(ismember(numNutrientsByCase,uniqueNumNutrientsByCase(i))),:,:)),1))');
    else
        stdErrorShannonByNumNutrients(i,:) = mean(mean(shannon(find(ismember(numNutrientsByCase,uniqueNumNutrientsByCase(i))),:,:),3),1);
    end
end
stdErrorShannonByNumNutrients = stdErrorShannonByNumNutrients./sqrt(numRandSims);

figure
imagesc(flipud(meanShannonByNumNutrients)/maxShannon)
colormap(autumn)
xlabel('Number of secreted metabolites')
ylabel('Number of nutrients')
set(gca,'FontSize',16,'YTick',[1:length(uniqueNumNutrientsByCase)],'YTickLabel',flipud(uniqueNumNutrientsByCase),'XTick',[1:maxSecMets+1],'XTickLabel',[0:maxSecMets])
title(['Relative shannon entropy; S = ' num2str(numSpecies) ', generalists = ' num2str(length(find(isGeneralist))) ', tP = ' num2str(transfProb) ', randSims = ' num2str(numRandSims)])
colorbar
caxis([0 0.65])

lineColors = [16/255,37/255,66/255;repelem([248/255,112/255,96/255],maxSecMets,1)];
figure
subplot(1,2,1)
for i = 1:maxSecMets+1
    errorbar(uniqueNumNutrientsByCase(2:end),meanGrowingSpeciesByNumNutrients(2:end,i),stdErrorGrowingSpeciesByNumNutrients(2:end,i),stdErrorGrowingSpeciesByNumNutrients(2:end,i),'LineWidth',2,'Color',lineColors(i,:))
    hold on
end
set(gca,'FontSize',16,'XTick',uniqueNumNutrientsByCase(2:end))
xlabel('Number of carbon sources')
ylabel('Species richness')
ylim([1,7])

subplot(1,2,2)
for i = 1:maxSecMets+1
    errorbar(uniqueNumNutrientsByCase(2:end),meanShannonByNumNutrients(2:end,i),stdErrorShannonByNumNutrients(2:end,i),stdErrorShannonByNumNutrients(2:end,i),'LineWidth',2,'Color',lineColors(i,:))
    hold on
end
set(gca,'FontSize',16,'XTick',uniqueNumNutrientsByCase(2:end))
xlabel('Number of carbon sources')
ylabel('Shannon entropy')
ylim([0.4,2.4])
legend({'Num Secreted Metabolites = 0','Num Secreted Metabolites > 0'})

% Analyze epistasis of species richness and shannon entropy
uniqueNumCSources = unique(sum(nutrientMat,2));
yieldsForEpi = squeeze(mean(yields,3));
[alphaDiff,shannonDiff,yieldDiff] = deal(zeros(size(nutrientMat,1),maxSecMets+1));
for p = 1:maxSecMets+1
    count = length(find(sum(nutrientMat,2)==1))+2;
    for i = 1:length(uniqueNumCSources)-2
        CSourceCurr = uniqueNumCSources(i+2);

        indicesHigherAll = find(sum(nutrientMat,2) == uniqueNumCSources(i+2)); % Indices of more complex combos
        indicesLowerAll = find(sum(nutrientMat,2) == uniqueNumCSources(i+1)); % Indices of next lowest combo
        for j = 1:length(indicesHigherAll)
            if i > 1
                indicesLower = indicesLowerAll(j*2-1:j*2); % Indices of specific lower combo
            else
                indicesLower = find(nutrientMat(indicesHigherAll(j),:))+1;
            end

            readsHigh = squeeze(relAbus(indicesHigherAll(j),:,p,:)); % Rows = species, cols = samples
            
            alphaHighFromEach = mean(sum(readsHigh>0,1));
            alphaHighFirstRand = sum(readsHigh(:,1)>0); % Take only first random sampling
            [alphaLowFromEach,shannonLowFromEach] = deal(zeros(length(indicesLower),1));
            [alphaLowFirstRand,shannonLowFirstRand] = deal(zeros(length(indicesLower),1));
            for mm = 1:length(indicesLower)
                alphaLowFromEach(mm) = mean(sum(squeeze(relAbus(indicesLower(mm),:,p,:)>0),1));
                
                alphaFirstRand = sum(squeeze(relAbus(indicesLower(mm),:,p,:)>0),1); % Take only first random sampling
                alphaLowFirstRand(mm) = alphaFirstRand(1);
                
                readsLowCurr = squeeze(relAbus(indicesLower(mm),:,p,:));
                shannonRand = zeros(size(readsLowCurr,2),1);
                for nn = 1:size(readsLowCurr,2)
                    for n = 1:size(readsLowCurr,1)
                        if relAbus(indicesLower(mm),n,p,nn) > 0
                            shannonRand(nn) = shannonRand(nn) - relAbus(indicesLower(mm),n,p,nn)*log2(relAbus(indicesLower(mm),n,p,nn));
                        end
                    end
                end
                shannonLowFromEach(mm) = mean(shannonRand);
                shannonLowFirstRand(mm) = shannonRand(1);
            end
            
            readsLowSum = squeeze(sum(relAbus(indicesLower,:,p,:),1));
            [alphaHigh,alphaLowSum] = deal(zeros(size(readsHigh,2),1));
            for jj = 1:length(alphaHigh)
                alphaHigh(jj) = length(find(readsHigh(:,jj)));
                alphaLowSum(jj) = length(find(readsLowSum(:,jj)));
            end

            alphaDiff(count,p) = alphaHighFirstRand - max(alphaLowFirstRand);
            
            [shannonExpected,shannonReal] = deal(zeros(size(readsHigh,2),1));
            for nn = 1:size(readsHigh,2)
                readsHighCurr = readsHigh(:,nn);
                s = sum(readsHighCurr);
                for n = 1:length(readsHighCurr)
                    if readsHighCurr(n) > 0
                        shannonReal(nn) = shannonReal(nn) - (readsHighCurr(n)/s)*log2(readsHighCurr(n)/s);
                    end
                end
            end
            
            shannonDiff(count,p) = shannonReal(1) - max(shannonLowFirstRand);
            
            yieldsHigh = squeeze(yields(indicesHigherAll(j),p,:));
            yieldsLow = squeeze(yields(indicesLower,p,:));
            yieldDiff(count,p) = mean(yieldsHigh) - mean(sum(yieldsLow,1)/2);

            count = count + 1;
        end
    end
end

[meanAlphaEpiByCSource,stdAlphaEpiByCSource,meanShannonEpiByCSource,stdShannonEpiByCSource,meanYieldEpiByCSource,stdYieldEpiByCSource] = deal(zeros(length(uniqueNumCSources(3:end)),maxSecMets+1));
for i = 1:length(uniqueNumCSources(3:end))
    meanAlphaEpiByCSource(i,:) = mean(alphaDiff(find(sum(nutrientMat,2) == uniqueNumCSources(i+2)),:),1);
    stdAlphaEpiByCSource(i,:) = std(alphaDiff(find(sum(nutrientMat,2) == uniqueNumCSources(i+2)),:),1);
    meanShannonEpiByCSource(i,:) = mean(shannonDiff(find(sum(nutrientMat,2) == uniqueNumCSources(i+2)),:),1);
    stdShannonEpiByCSource(i,:) = std(shannonDiff(find(sum(nutrientMat,2) == uniqueNumCSources(i+2)),:),1);
    meanYieldEpiByCSource(i,:) = mean(yieldDiff(find(sum(nutrientMat,2) == uniqueNumCSources(i+2)),:),1);
    stdYieldEpiByCSource(i,:) = std(yieldDiff(find(sum(nutrientMat,2) == uniqueNumCSources(i+2)),:),1);
end

% Compile data to export manually
epiCRM = struct();
epiCRM.alphaDiff = alphaDiff;
epiCRM.shannonDiff = shannonDiff;
epiCRM.yieldDiff = yieldDiff;
epiCRM.uniqueNumNutrientsByCase = uniqueNumNutrientsByCase;
epiCRM.meanShannonByNumNutrients = meanShannonByNumNutrients;
epiCRM.stdErrorShannonByNumNutrients = stdErrorShannonByNumNutrients;
epiCRM.meanGrowingSpeciesByNumNutrients = meanGrowingSpeciesByNumNutrients;
epiCRM.stdErrorGrowingSpeciesByNumNutrients = stdErrorGrowingSpeciesByNumNutrients;

end

%% Plot yields with environmental complexity

figure
yi = mean(yields(:,2:end,:),3);
boxplot(reshape(yi',size(yi,2)*size(yi,1),1),repelem(sum(nutrientMat,2),size(yi,2)))
set(findobj(gca,'type','line'),'linew',2)
set(gca,'FontSize',16)
xlabel('Number of carbon sources (nCS)')
ylabel(['Biomass (OD600)'])
ylim([0,0.8])

yieldsForCoop = mean(squeeze(mean(yields,2)),2)*8e8;
colors = parula(10);