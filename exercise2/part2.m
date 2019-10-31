%% Exercise 2 - Part 2
clear, clc, close all

% change the solver 
changeCobraSolver('cplex_direct');

%% Load the model
% matTFA contains both a smallEcoli.mat and small_ecoli.mat model
% here we're explicitly using small_ecoli.mat from Moodle
tmp = load('small_ecoli.mat');
model = tmp.model_red;
clear tmp

%% Load the thermodynamics database
tmp = load('thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp

%% Perfom an FBA on the model (optimize for growth)
FBAsolution = optimizeCbModel(model);

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

O2 = [-20 0]; %aerobic or anaerobic reaction
Rxn = {'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e'};%vector of exchange reactions
for i = 1:4 %loop for 4 exchange reactions
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, -10, 'l'); %boundary setting for carbon source of interest
    for j = 1:2 %loop for aerobic/anaerobic reactions
        this_tmodel = changeTFArxnBounds (this_tmodel, 'DM_o2_e', O2(j), 'l'); %boudary setting for aerobic/anaerobic reaction
        TFAsolution = solveTFAmodelCplex(this_tmodel);
        YieldTFA (i,j) = TFAsolution.val;
    end
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, 0, 'l');%reinitialise boundary setting of carbon source
end

%% load concentration data

mets = readtable('metabolomics_data.csv');
% SD sometimes interpreted as string (depends on PC) -- convert if needed
if ischar(mets.StandardDeviation)
    mets.StandardDeviation = str2double(mets.StandardDeviation);
end

%% adjust model

% find the index for each metabolite
% NB: find_cell silently ignores non-matches; use this function instead
% https://uk.mathworks.com/matlabcentral/answers/142925-matching-texts-in-two-cell-arrays#answer_145977
%mets.modelIndex = cellfun(@(a) strmatch(a,this_tmodel.varNames),strcat('LC_', mets.modelID),'Uniform',false);
%mets.modelIndex = cell2mat(mets.modelIndex);
% actually this thing returns cells which are a bit sad

% find the model index
for met = 1:size(mets,1)
    try
        mets.modelIndex(met) = find_cell(strcat('LC_', mets.modelID(met)), this_tmodel.varNames);
    catch
        mets.modelIndex(met) = NaN;
    end  
end

% filter metabolites in and not in model
notmets = mets(isnan(mets.modelIndex),:);
mets = mets(~isnan(mets.modelIndex),:);

% save both lists
writetable(mets(:,1:2), [pwd '/out/mets.csv']);
writetable(notmets(:,1:2), [pwd '/out/notmets.csv']);

% calculate upper and lower bounds
for met = 1:size(mets,1)
    if isnan(mets.StandardDeviation(met))
        % when we don't have SD, use 50-150% of conc
        mets.C_ub(met) = 1.5 * mets.ConcentrationInMmol(met);
        mets.C_lb(met) = 0.5 * mets.ConcentrationInMmol(met);
    else
        % total range +/-1 SD around given conc
        mets.C_ub(met) = mets.ConcentrationInMmol(met) + mets.StandardDeviation(met);
        mets.C_lb(met) = mets.ConcentrationInMmol(met) - mets.StandardDeviation(met);
    end    
end

% set bounds
this_tmodel.var_lb(mets.modelIndex) = log(mets.C_lb);
this_tmodel.var_ub(mets.modelIndex) = log(mets.C_ub);

%% Perform TFA for growth with concentration data 

O2 = [-20 0]; %aerobic or anaerobic reaction
Rxn = {'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e'};%vector of exchange reactions
for i = 1:4 %loop for 4 exchange reactions
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, -10, 'l'); %boundary setting for carbon source of interest
    for j = 1:2 %loop for aerobic/anaerobic reactions
        this_tmodel = changeTFArxnBounds (this_tmodel, 'DM_o2_e', O2(j), 'l'); %boudary setting for aerobic/anaerobic reaction
        TFAsolution_w_con = solveTFAmodelCplex(this_tmodel);
        YieldTFA_w_con (i,j) = TFAsolution_w_con.val;
    end
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, 0, 'l');%reinitialise boundary setting of carbon source
end