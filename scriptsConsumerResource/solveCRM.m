% solveCRM.m 
%
% Runs consumer resource models from CRM.m based on experimentally-derived
% species-specific nutrient utilization probabilities and FBA-predicted
% exchange. Requires COBRA toolbox for FBA predictions.
%
% Outputs: 
%   -relAbus: an N x S matrix containing organism relative abundances
%    in each environmental condition, where N is the number of conditions
%    and S is the number of organisms. 
%
%   - if plotIndividuals is set to 1, a figure showing the time-dependent
%   trajectory of organism and resource abundances, as well as final
%   species abundances for the final tested nutrient condition. Here, since
%   the environmental conditions of com4 are tested, the figure will
%   display these data for the last condition, i.e. the condition
%   in which all 32 carbon sources are provided to the organisms.
%
% Alan R. Pacheco, 04/23/2020, modified 01/26/2021

%% Define filenames and select species and nutrients

clearvars

selectSpecies = {'Bs','Me','So','Pa'}; % B. subtilis, M. extorquens, S. oneidensis, and P. aeruginosa (corresponds to 'com4')
crossfeeding = 1; % If 0, D matrix is all zeros
premadeDMatrix = 1; % As a demo, a premade D Matrix has been generated for com4. If desired, set to 1 and see line 40
plotIndividuals = 1; % If 1, will plot simulation data for last individual community
monod = 1; % If 1, uses Monod dynamics in species growth and metabolite consumption
sequentialConsumption = 0; % If 1, organisms consume metabolites sequentially
timestep = 0.01; % Time resolution of simulation
totalExpLength = 288; % Total length of the experiment, in hours
dilutionTime = 48; % Frequency of dilutions, in hours

% Load monoculture and nutrient data
monoGrowthDataFile = '../dataExperimental/stockGrowthStrainsJan2021'; % Contains species-specific monoculture growth data for generating C matrix
load 32CSNutrients.mat % Contains list of nutrients and combinations
conditionsToTest = [1:size(nutrientComboMap,1)-1]; % Defines which nutrient combinations will be simulated. Here, all 63 condition from com4 are selected

% Generate D matrix, either by loading a premade example or by using FBA
if premadeDMatrix
    load 'DMatrixBsMeSoPa.mat';
else
    DFull = generateDMatrix('modelsBsMeSoPa.mat','nutrients32FBA.mat','minMed.mat');
end

%% Initialize data structures and run CRM
[growingSpecies,yields,rho] = deal(zeros(size(nutrientComboMap,1),1));
relAbus = zeros(size(nutrientComboMap,1),length(selectSpecies));
for q = 1:length(conditionsToTest)

selectNutrientCombo = conditionsToTest(q);
disp([num2str(q),'/',num2str(size(nutrientComboMap,1))])
    
% Get nutrient concentrations
selectNutrients = nutrients(find(nutrientComboMap(selectNutrientCombo,:)))
selectNutrientsIndices = find(ismember(nutrients,selectNutrients));
nutrientMassPerVolume = MW.*nutrientConcs; % g/L

% Set simulation parameters

global S M C D g w l m k d kR monodTerm sequential initConcs % N is CFU/mL, R is g/mL

monodTerm = monod;
sequential = sequentialConsumption;

S = length(selectSpecies); % Number of species
g = ones(1,S).*1; % Conversion factor from energy uptake to growth rate (1/energy) SET TO 1
m = ones(1,S).*0.025; % Minimal energy uptake for maintenance of each species (energy/hour)

M = length(selectNutrients); % Number of limiting nutrients
w = ones(1,length(nutrients))*2.5e9; % Energy content of resource a (energy/g)
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

C = generateCMatrix(monoGrowthDataFile,selectSpecies,simulationNutrients); % Uptake rate per unit concentration of resource a by species i (mL/hour)
C = C.*5e-5;
rho(q) = (mean(C(:))^2)/(mean(C(:))^2+var(C(:)));

% Set initial conditions and solve system of ODEs

N0 = ones(1,S).*6e6; % approximately OD 0.5 diluted to 5µl/300µl

R0 = ones(1,length(simulationNutrients)).*1e-10; % Baseline must be greater than 0 to allow for D matrix turnover
totalNutrient = 1.5; %g/L
R0(find(ismember(simulationNutrients,selectNutrients))) = totalNutrient./(length(selectNutrients));  % g/L, scaled but otherwise equal since CRM does not care about number of carbon atoms
k = R0.*0; % External continuous supply of resource a (grams/mL/hour)
d = Inf;
kR = ones(1,length(simulationNutrients)).*1e4;
initConcs = R0;

I0 = [N0 R0];

timeRange = [0:timestep:dilutionTime];
[T,X]=ode15s('CRM',timeRange,I0);

dilutionTimes = [0:dilutionTime:totalExpLength-dilutionTime];

[N,R] = deal([]);
N = [N;X(:,1:S)];
R = [R;X(:,S+1:end)];

for t = 2:length(dilutionTimes)
    
    Nn = X(end,1:S)*(10/300); % Dilute organisms 10µl into 300µl
    Rn = X(end,S+1:end)*(10/300) + R0; % Dilute remaining nutrient and add in R0
    
    timeRange = [T(end):timestep:T(end)+dilutionTime];
    [Tn,X]=ode15s('CRM',timeRange,[Nn Rn]);
    
    N = [N;X(:,1:S)];
    R = [R;X(:,S+1:end)];
    T = [T;Tn];
    
end

% Get relative abundances

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