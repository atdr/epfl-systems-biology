%% Exercise 2 - Part 1
clear, clc, close all

% load model
tmp = load('small_ecoli.mat');
model = tmp.model_red;
clear tmp

% load thermodynamics database
tmp = load('thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp

%% perform FBA (optimize for growth)

% oxygen flux for aerobic or anaerobic conditions
O2 = [-20 0];

% demand reactions
Rxn = {'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e'};

% loop through reactions
for i = 1:4
    
    % set carbon source boundary
    model = changeRxnBounds(model, Rxn{i}, -10, 'l');
    
    % loop through aerobic/anaerobic conditions
    for j = 1:2
        
        % set oxygen boundary
        model = changeRxnBounds(model, 'DM_o2_e', O2(j), 'l');
        
        % do FBA solution and save result
        FBAsolution = optimizeCbModel(model, 'max');
        YieldFBA(i,j) = FBAsolution.f;
    end
    
    % reset carbon source to 0
    model = changeRxnBounds(model, Rxn{i}, 0, 'l');
end

%% prepare model for TFA and convert

prepped_model = prepModelforTFA(model, ReactionDB, model.CompartmentData);

min_obj = 0.0;
tmp = convToTFA(prepped_model, ReactionDB, [], 'DGo', [], min_obj);

% Add net flux variables, which are equal to forwards flux - backwards flux
% The netflux variables work similar to the FBA fluxes
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp);

%% perform TFA (optimise for growth)

% loop through reactions
for i = 1:4
    
    % set carbon source boundary
    this_tmodel = changeTFArxnBounds(this_tmodel, Rxn{i}, -10, 'l');
    
    % loop through aerobic/anaerobic conditions
    for j = 1:2
        
        % set oxygen boundary
        this_tmodel = changeTFArxnBounds(this_tmodel, 'DM_o2_e', O2(j), 'l');
        
        % do TFA solution and save result
        TFAsolution = solveTFAmodelCplex(this_tmodel);
        YieldTFA(i,j) = TFAsolution.val;
    end
    
    % reset carbon source to 0
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, 0, 'l');
end
