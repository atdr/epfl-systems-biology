clc
addpath(genpath('C:\AAA\EPFL Y3\Principles and Application of Systems biology\MATLAB\matTFA-master'));
addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab'))
changeCobraSolver('cplex_direct');
%% Load the model
tmp = load('small_ecoli.mat');
model = tmp.model_red;
clear tmp
%% Load the thermodynamics database
tmp = load('thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp
%% Perfom an FBA on the model (optimize for growth)
O2 = [-20 0]; %aerobic or anaerobic reaction
Rxn = {'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e'};%vector of exchange reactions
for i = 1:4 %loop for 4 exchange reactions
    model = changeRxnBounds (model, Rxn{i}, -10, 'l'); %boundary setting for carbon source of interest
    for j = 1:2 %loop for aerobic/anaerobic reactions
        model = changeRxnBounds (model, 'DM_o2_e', O2(j), 'l'); %boudary setting for aerobic/anaerobic reaction
        FBAsolution = optimizeCbModel (model, 'max'); 
        YieldFBA (i,j) = FBAsolution.f;
    end
    model = changeRxnBounds (model, Rxn{i}, 0, 'l');%reinitialise boundary setting of carbon source
end

%% Prepare the model for TFA
prepped_model = prepModelforTFA(model, ReactionDB, model.CompartmentData);

%% Convert to TFA
min_obj = 0.0;
tmp = convToTFA(prepped_model, ReactionDB, [], 'DGo', [], min_obj);

% Add net flux variables, which are equal to forwards flux - backwards flux
% The netflux variables work similar to the FBA fluxes
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp);

%% Perfom TFA on the model (optimise for growth)
%soltFA = solveTFAmodelCplex(this_tmodel);
for i = 1:4 %loop for 4 exchange reactions
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, -10, 'l'); %boundary setting for carbon source of interest
    for j = 1:2 %loop for aerobic/anaerobic reactions
        this_tmodel = changeTFArxnBounds (this_tmodel, 'DM_o2_e', O2(j), 'l'); %boudary setting for aerobic/anaerobic reaction
        TFAsolution = solveTFAmodelCplex(this_tmodel);
        YieldTFA (i,j) = TFAsolution.val;
    end
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, 0, 'l');%reinitialise boundary setting of carbon source
end

% %% We add some generic data for concentrations
% %metNames = {'adp_c', 'atp_c'};
% C_lb = [1e-06, 1e-03]';
% C_ub = [7e-04, 5e-02]';
% LC_varNames = {'LC_adp_c', 'LC_atp_c'};
% % find the indices of these variables in the variable names of the tfa
% 
% id_LC_varNames = find_cell(LC_varNames, this_tmodel.varNames);
% % Set to the model these log-concentration values
% this_tmodel.var_lb(id_LC_varNames) = log(C_lb);
% this_tmodel.var_ub(id_LC_varNames) = log(C_ub);
% 
% %% Optimise for growth with concentration data 
% soltFA_w_concentrations = solveTFAmodelCplex(this_tmodel);