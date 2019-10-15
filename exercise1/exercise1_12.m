%% Exercise 1 - Parts 1 + 2

clear, clc, close all
load('ecoli_core_model.mat')

%% setup

% set biomass rxn as objective function
model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');

% set glucose exchange to -10 mmol/gDW/h
model = changeRxnBounds(model,'EX_glc(e)',-10.0,'l');

% substrates to consider
P1substrates = {'EX_glc(e)'; 'EX_lac-D(e)'; 'EX_fru(e)'; 'EX_etoh(e)'};
P2substrates = {'EX_glc(e)'; 'EX_lac-D(e)'; 'EX_gln-L(e)'; 'EX_glu-L(e)'; 'EX_fru(e)'; 'EX_pyr(e)'; 'EX_ac(e)'; 'EX_etoh(e)'};

%% aerobic

% define aerobic conditions
model = changeRxnBounds(model,'EX_o2(e)',-20.0,'l');

% solve for P1 substrates
[model, P1solution(1).data] = substrateSolver(model,P1substrates,-10.0);
P1solution(1).label = 'aerobic';

% solve for P2 substrates
[model, P2solution(1).data] = substrateSolver(model,P2substrates,-10.0);
P2solution(1).label = 'aerobic';

%% anaerobic

% define anaerobic conditions
model = changeRxnBounds(model,'EX_o2(e)',0,'l');

% solve for P1 substrates
[model, P1solution(2).data] = substrateSolver(model,P1substrates,-10.0);
P1solution(2).label = 'anaerobic';

% solve for P2 substrates
[model, P2solution(2).data] = substrateSolver(model,P2substrates,-10.0);
P2solution(2).label = 'aerobic';

%% Part 1 - export glucose substrate data

for cnt = 1:numel(P1solution)
    writetable(table(model.rxns,P1solution(cnt).data(1).solution.x), [pwd '/out/' P1solution(cnt).label '_flux_data.csv']);
end

%% functions

function [model, t] = substrateSolver(model,substrates,flux)

    % for each substrate
    for cnt = 1:numel(substrates)

        % create vector for substrate fluxes (all zero except one)
        fluxes = zeros(numel(substrates),1);
        fluxes(cnt) = flux;

        % set all fluxes
        model = changeRxnBounds(model,substrates,fluxes,'l');

        % solve optimisation problem
        t(cnt).solution = optimizeCbModel(model,'max');
        t(cnt).label = substrates(cnt);
        t(cnt).sol_f = t(cnt).solution.f;
        
        % perform gene knockout analysis
        [grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model);
        
        % genes for which the growth rate falls below 10%
        t(cnt).essentialGenes = grRatio < 0.1; % logical array
        
        % count essential genes
        t(cnt).essentialGenesCount = sum(t(cnt).essentialGenes);
    end
    
end
