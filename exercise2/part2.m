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

%% prepare model for TFA and convert

prepped_model = prepModelforTFA(model, ReactionDB, model.CompartmentData);

min_obj = 0.0;
tmp = convToTFA(prepped_model, ReactionDB, [], 'DGo', [], min_obj);

% Add net flux variables, which are equal to forwards flux - backwards flux
% The netflux variables work similar to the FBA fluxes
% NF_rxn = F_rxn - B_rxn
this_tmodel = addNetFluxVariables(tmp);

%% Perfom TFA on the model (without metabolomics data)

% oxygen flux for aerobic or anaerobic conditions
O2 = [-20 0];

% substrate demand reactions
Rxn = {'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e'};

% loop through substrates
for i = 1:numel(Rxn)
    
    % set carbon source boundary
    this_tmodel = changeTFArxnBounds(this_tmodel, Rxn{i}, -10, 'l');
    
    % loop through aerobic/anaerobic conditions
    for j = 1:numel(O2)
        
        % set oxygen boundary
        this_tmodel = changeTFArxnBounds(this_tmodel, 'DM_o2_e', O2(j), 'l');
        
        % do TFA solution and save result
        TFAsolution = solveTFAmodelCplex(this_tmodel);
        YieldTFA(i,j) = TFAsolution.val;
    end
    
    % reset carbon source to 0
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, 0, 'l');
end

%% load concentration data

mets = readtable('metabolomics_data.csv');
% SD sometimes interpreted as string / cell (depends on MATLAB version)
% convert if needed
if ~isa(mets.StandardDeviation,'double')
    mets.StandardDeviation = str2double(mets.StandardDeviation);
end

%% input metabolomics data

% find the index for each metabolite
% NB: find_cell silently ignores non-matches; use this function instead
% https://uk.mathworks.com/matlabcentral/answers/142925-matching-texts-in-two-cell-arrays#answer_145977
%mets.modelIndexLC = cellfun(@(a) strmatch(a,this_tmodel.varNames),strcat('LC_', mets.modelID),'Uniform',false);
%mets.modelIndexLC = cell2mat(mets.modelIndexLC);
% actually this thing returns cells which are a bit sad

% find index of LC variable for each metabolite
for met = 1:size(mets,1)
    try
        mets.modelIndexLC(met) = find_cell(strcat('LC_', mets.modelID(met)), this_tmodel.varNames);
    catch
        mets.modelIndexLC(met) = NaN;
    end
end

% filter metabolites in and not in model
notmets = mets(isnan(mets.modelIndexLC),:);
mets = mets(~isnan(mets.modelIndexLC),:);

% save both lists (just first 2 columns)
writetable(mets(:,1:2), [pwd '/out/mets.csv']);
writetable(notmets(:,1:2), [pwd '/out/notmets.csv']);

% calculate upper and lower concentration bounds
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

% save model before applying metabolomics data
tmodel_no_metabolomics = this_tmodel;

% set metabolite concentration bounds
this_tmodel.var_lb(mets.modelIndexLC) = log(mets.C_lb);
this_tmodel.var_ub(mets.modelIndexLC) = log(mets.C_ub);

%% Perform TFA for growth with concentration data

% loop through substrates
for i = 1:numel(Rxn)
    
    % set carbon source boundary
    this_tmodel = changeTFArxnBounds(this_tmodel, Rxn{i}, -10, 'l');
    
    % loop through aerobic/anaerobic conditions
    for j = 1:numel(O2)
        
        % set oxygen boundary
        this_tmodel = changeTFArxnBounds(this_tmodel, 'DM_o2_e', O2(j), 'l');
        
        % do TFA solution and save result
        TFAsolution_w_con = solveTFAmodelCplex(this_tmodel);
        YieldTFA_w_con(i,j) = TFAsolution_w_con.val;
    end
    
    % reset carbon source to 0
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, 0, 'l');
end

%% metabolite concentration variability analysis

% with metabolomics data
metConcLims_with = doMetConcVarAnalysis(this_tmodel, Rxn, O2, YieldTFA_w_con, mets);
% NB: these are the log bounds that we *imposed*
% (compare with values from optimisation)
metConcLims_with_imposed = real([log(mets.C_lb) log(mets.C_ub)]); % use real() since log(0) = -Inf + 0i according to Matlab

% without metabolomics data
metConcLims_without = doMetConcVarAnalysis(tmodel_no_metabolomics, Rxn, O2, YieldTFA, mets);
% NB: these are the log bounds that were originally in the model
metConcLims_without_imposed = [tmodel_no_metabolomics.var_lb(mets.modelIndexLC) tmodel_no_metabolomics.var_ub(mets.modelIndexLC)];

%% analyse data

% find ranges with and without metabolomics data
for c = 1:numel(Rxn)
    metConcLims_with_without_ranges{c,1} = [metConcLims_with{c,1}(:,2)-metConcLims_with{c,1}(:,1),...
        metConcLims_without{c,1}(:,2)-metConcLims_without{c,1}(:,1)];
    
    % export for Escher
    writetable(table(mets.modelID,metConcLims_with_without_ranges{c,1}), [pwd '/out/metConcLims_with_without_ranges_' Rxn{c} '.csv']);
end

% calculate difference between limits and bounds
for c = 1:numel(Rxn)
    metConcLims_with_change{c,1} = metConcLims_with{c,1} - metConcLims_with_imposed;
    metConcLims_without_change{c,1} = metConcLims_without{c,1} - metConcLims_without_imposed;
    
    % export table
    writetable(table(mets.modelID,metConcLims_with_change{c,1},metConcLims_without_change{c,1}), [pwd '/out/metConcLims_change_' Rxn{c} '.csv']);
end

%% function space

function metConcLims = doMetConcVarAnalysis(this_tmodel, Rxn, O2, Yields, mets)
% performs metabolite concentration variability analysis

% previous objective variable (where 1 is in this_tmodel.f)
objVar = 'F_Ec_biomass_iJO1366_WT_53p95M';
objVarIndex = find_cell(objVar, this_tmodel.varNames);

% remove this objective
this_tmodel.f(objVarIndex) = 0;

% data output cell
metConcLims = {};

% loop through reactions
for i = 1:numel(Rxn)
    
    % set carbon source boundary
    this_tmodel = changeTFArxnBounds(this_tmodel, Rxn{i}, -10, 'l');
    
    % aerobic conditions only
    j = 1;
    this_tmodel = changeTFArxnBounds(this_tmodel, 'DM_o2_e', O2(j), 'l');
    
    % set upper and lower bounds of previous objective variable to
    % optimum previously found
    this_tmodel.var_ub(objVarIndex) = Yields(i,j);
    this_tmodel.var_lb(objVarIndex) = Yields(i,j);
    
    % initialise output for metabolite concentrations
    metConcLims{i,j} = [];
    
    % loop through metabolites
    for k = 1:size(mets,1)
        %for k = 1:1
        
        % set objective to maximise LC of the metabolite
        this_tmodel.f(mets.modelIndexLC(k)) = 1;
        % solve
        TFAsolution_max = solveTFAmodelCplex(this_tmodel);
        
        % set objective to minimise LC of the metabolite
        this_tmodel.f(mets.modelIndexLC(k)) = -1;
        % solve
        TFAsolution_min = solveTFAmodelCplex(this_tmodel);
        
        % save the results
        metConcLims{i,j} = [metConcLims{i,j}; TFAsolution_min.x(mets.modelIndexLC(k)) TFAsolution_max.x(mets.modelIndexLC(k))];
        
        % now we're done, remove this metabolite as objective
        this_tmodel.f(mets.modelIndexLC(k)) = 0;
        
    end
    
    % reset carbon source to 0
    this_tmodel = changeTFArxnBounds (this_tmodel, Rxn{i}, 0, 'l');
end

end
