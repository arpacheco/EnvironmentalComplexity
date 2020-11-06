% solveCRM.m 
%
% Runs consumer resource models from CRM.m based on experimentally-derived
% species-specific nutrient utilization probabilities and FBA-predicted
% exchange. Requires COBRA toolbox for FBA predictions.
%
% Alan R. Pacheco, 04/23/2020

%% Define filenames and select species and nutrients

clearvars

selectSpecies = {'Bs','Me','So','Pa'}; % B. subtilis, M. extorquens, S. oneidensis, and P. aeruginosa (corresponds to 'com4')
crossfeeding = 1; % If 0, D matrix is all zeros
plotIndividuals = 1; % If 1, will plot simulation data for last individual community
monod = 1; % If 1, uses Monod dynamics in species growth and metabolite consumption
sequentialConsumption = 0; % If 1, organisms consume metabolites sequentially
timestep = 0.01;

% Load monoculture data
monoGrowthDataFile = 'BsMePaSo_PROld.xlsx';
load 32CSNutrients.mat
allPossibleSpecies = {'Bs','Me','So','Pa'};

% Generate D matrix from model files
DFull = generateDMatrix('modelsBsMeSoPa.mat','nutrients32FBA.mat','minMed.mat');

[growingSpecies,yields] = deal(zeros(64,1));
relAbus = zeros(size(nutrientComboMap,1),length(selectSpecies));
for q = 1:size(nutrientComboMap,1)-1

selectNutrientCombo = q;
disp([num2str(q),'/',num2str(size(nutrientComboMap,1))])
    
%% Get nutrient concentrations
selectNutrients = nutrients(find(nutrientComboMap(selectNutrientCombo,:)))
selectNutrientsIndices = find(ismember(nutrients,selectNutrients));
nutrientMassPerVolume = MW.*nutrientConcs; % g/L

%% Set simulation parameters

global S M C D g w l m k d kR monodTerm sequential initConcs
% N is CFU/mL, R is g/mL

monodTerm = monod;
sequential = sequentialConsumption;

S = length(selectSpecies); % Number of species
g = ones(1,S).*1; % Conversion factor from energy uptake to growth rate (1/energy) SET TO 1
m = ones(1,S).*0.025; % Minimal energy uptake for maintenance of each species (energy/hour)

M = length(selectNutrients); % Number of limiting nutrients
w = ones(1,length(nutrients))*1e9; % Energy content of resource a (energy/g)
l = ones(1,length(nutrients)).*0.8; % Leakage fraction for resource a (unitless) SET TO 0.8

% Use the D matrix to get all the metabolites in the simulation
D = DFull;
if crossfeeding
    D(find(D>1)) = 1-(1/D(find(D>1))); % ratios should max out at 1
    
    transformedNutrients = nutrients(find(sum(sum(D(:,find(ismember(nutrients,selectNutrients)),:),3),2)));
    simulationNutrients = [selectNutrients;setdiff(transformedNutrients,selectNutrients)];
    simulationNutrientIndices = [find(ismember(nutrients,selectNutrients));find(ismember(nutrients,setdiff(transformedNutrients,selectNutrients)))];
    for s = 1:S
        DCurr = D(:,:,s);
        DCurr = DCurr(simulationNutrientIndices,simulationNutrientIndices)'; % Flip matrix since it's b converting to a
        DCurr = DCurr./sum(DCurr,1); % rows must add to 1
        if s == 1
            DNew = DCurr;
        else
            DNew(:,:,s) = DCurr;
        end
    end
    D = DNew;
    D(find(isnan(D))) = 0;
else
    simulationNutrients = selectNutrients;
    D = zeros(length(simulationNutrients),length(simulationNutrients),S);
end
simulationNutrientIndices = find(ismember(nutrients,simulationNutrients));
w = w(simulationNutrientIndices);
l = l(simulationNutrientIndices);

[C,growthData] = generateCMatrix(monoGrowthDataFile,selectSpecies,nutrients,simulationNutrients); % Uptake rate per unit concentration of resource a by species i (mL/hour)
C = C.*1e-5;

%% Set initial conditions and solve system of ODEs

N0 = ones(1,S).*6e6; % approximately OD 0.5 diluted to 5µl/300µl

R0 = ones(1,length(simulationNutrients)).*1e-1; % Baseline must be greater than 0 to allow for D matrix turnover
totalNutrient = 1.5; %g/L
R0(find(ismember(simulationNutrients,selectNutrients))) = totalNutrient./(length(selectNutrients));  % g/L, scaled but otherwise equal since CRM does not care about number of carbon atoms
k = R0/48; % External supply of resource a (grams/mL/hour), equals R0 divided by 48 hour replenishment period
d = 1/(10/(300*48)); % Timescale for dilution (hour). 10ul/300ul every 48h
kR = ones(1,length(simulationNutrients)).*3000; %250
initConcs = R0;

I0 = [N0 R0];

timeRange = [0:timestep:288];
[T,X]=ode15s('CRM',timeRange,I0);

%% Get relative abundances

N = X(:,1:S);
R = X(:,S+1:end);

abundances = N(end,:);
abundances(find(abundances < 0)) = 0;

relAbu = abundances./sum(abundances);

relAbus(q,:) = relAbu;
yields(q) = sum(abundances)/8e8;
growingSpecies(q) = length(find(abundances > N0));

end

relAbus(find(isnan(relAbus))) = 0;

%% Plot

close all

if plotIndividuals
    
c=[48,44,129;66,154,203;170,188,116;255,184,111]./255;

if crossfeeding && length(setdiff(transformedNutrients,selectNutrients)) > 0 
    numplots = 4; else numplots = 3; 
end
numplots=4;

figure
subplot(1,numplots,1)
for i = 1:S
    plot(T,N(:,i),'LineWidth',4,'Color',c(i,:))
    hold on
end

legend(selectSpecies)
xlabel('Time (h)')
ylabel('Species abundance (CFU/mL)')
yl = ylim;
yyaxis right
set(gca, 'YTick', round([0:yl(2)/8e8/5:yl(2)/8e8],2), 'ylim', [0, yl(2)/8e8])
ylabel('Species abundance (OD600)')
xl = xlim;
set(gca,'FontSize',16,'XTick',[0:48:xl(2)])
title('Growth yields')

subplot(1,numplots,2)
b = bar([relAbu;relAbu],'stacked');
b(1).Parent.Parent.Colormap = c(1:S,:);
xlim([0.5,1.5])
ylim([0,1])
legend(selectSpecies)
ylabel('Relative abundance')
set(gca,'FontSize',16,'XTickLabel',[])
title('Endpoint relative abundances')

subplot(1,numplots,3)
for a = 1:length(selectNutrients)
    plot(T,R(:,a),'LineWidth',4)
    hold on
end
legend(selectNutrients)
xlabel('Time (h)')
ylabel('Resource abundance (g/mL)')
xl = xlim;
set(gca,'FontSize',16,'XTick',[0:48:xl(2)])
title('Supplied resource abundance')

if crossfeeding && length(setdiff(transformedNutrients,selectNutrients)) > 0 
subplot(1,numplots,4)
for a = 1:length(setdiff(transformedNutrients,selectNutrients))
    plot(T,R(:,length(selectNutrients)+a),'LineWidth',4)
    hold on
end
legend(setdiff(transformedNutrients,selectNutrients))
xlabel('Time (h)')
ylabel('Resource abundance (g/mL)')
xl = xlim;
set(gca,'FontSize',16,'XTick',[0:48:xl(2)])
title('Transformed resource abundance')
end

end