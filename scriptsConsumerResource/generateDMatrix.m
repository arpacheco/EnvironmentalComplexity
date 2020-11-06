function D = generateDMatrix(modelsFileName,nutrientsFileName,minMedFileName)

%% Initialize toolbox and define variables
initCobraToolbox

load(modelsFileName)
load(nutrientsFileName)
load(minMedFileName)

modelNames = fieldnames(models);

D = zeros(length(nutrients),length(nutrients),length(modelNames));

%% Populate D matrix

for m = 1:length(modelNames)

	model = models.(modelNames{m});

	for n = 1:length(nutrients)
		modelsDefine.(modelNames{m}) = model;
		modelsNew = defineMedium([nutrients{n}; minMed],modelsDefine);
		modelNew = modelsNew.(modelNames{m});

		FBAsoln = optimizeCbModel(modelNew,'max','one'); %minimize taxicab norm: min |v| s.t.: S*v = b, c'v = f, lb <= v <= ub

		% If the model grew, match the secreted metabolites to reactions and get fluxes
		if FBAsoln.f >= 0
			excRxns = find(strncmp('EX_',modelNew.rxns,3));

			absRxn = intersect(excRxns,find(modelNew.S(find(ismember(modelNew.mets,nutrients{n})),:)));
			secRxns = intersect(excRxns,find(FBAsoln.x > 0));

			for r = 1:length(secRxns)
				secMetIndexInD = find(ismember(nutrients,modelNew.mets{find(modelNew.S(:,secRxns(r)))}));
				if ~isempty(secMetIndexInD)
					D(n,secMetIndexInD,m) = abs(FBAsoln.x(absRxn)/FBAsoln.x(secRxns(r)));
				end
			end
		end
	end
end