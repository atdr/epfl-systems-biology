%% Exercise 2 - Part 3
clear, clc, close all
% requires the Parallel Computing Toolbox (!)

% change the solver
changeCobraSolver('cplex_direct');

%% Load the model
tmp = load('small_ecoli.mat');
model = tmp.model_red;
clear tmp

%% Data for the Flux Analysis

O2 = [-20 0];
O2_label = {'aero', 'anaero'};
Rxn = {'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e' };

% pre-calculate length of for loops
end_i = length(Rxn);
end_j = length(O2);
end_k = length(model.rxns);

%% FBA Flux Analysis

% initialise output matrices
Fluxmin_FBA = NaN(end_i, end_j, end_k);
Fluxmax_FBA = NaN(end_i, end_j, end_k);

% grab environment to pass to parallel workers
environment = getEnvironment();

parfor i = 1:end_i
    
    % set up environment for each worker
    restoreEnvironment(environment);
    changeCobraSolver('cplex_direct');
    
    % copy model
    model1 = model;

    % Set the substrate flux
    model1 = changeRxnBounds (model1, Rxn(i), -10, 'l');
    
    for j = 1:end_j
        
        % Set the O2 flux
        model1 = changeRxnBounds (model1, 'DM_o2_e', O2(j), 'l');
        
        % Determine the maximum biomass production and set it as constant
        FBAsolution = optimizeCbModel (model1, 'max');
        model1 = changeRxnBounds (model1, 'Ec_biomass_iJO1366_WT_53p95M', FBAsolution.f, 'b');
        
        % Flux analysis
        for k = 1:end_k
            
            % Set the objective reaction
            model1 = changeObjective (model1, model1.rxns (k));
            
            % Perform the flux analysis
            FBAsolutionmin = optimizeCbModel (model1, 'min');
            Fluxmin_FBA (i, j, k) = FBAsolutionmin.f;
            FBAsolutionmax = optimizeCbModel (model1, 'max');
            Fluxmax_FBA (i, j, k) = FBAsolutionmax.f;

        end
        
        % print an update
        fprintf('Finished FBA flux analysis for:\t%s\t%s\n',Rxn{i},O2_label{j});
    end
    
    model1 = changeRxnBounds (model1, Rxn(i), 0, 'l');

end

%% Load the thermodynamics database

tmp = load('thermo_data.mat');
ReactionDB = tmp.DB_AlbertyUpdate;
clear tmp

%% Prep model for TFA analysis

prepped_model = prepModelforTFA(model, ReactionDB, model.CompartmentData);
min_obj = 0.0;
tmp = convToTFA(prepped_model, ReactionDB, [], 'DGo', [], min_obj);
TFA_model = addNetFluxVariables(tmp);
clear tmp

%% TFA Flux Analysis

% initialise output matrices
Fluxmin_TFA = NaN(end_i, end_j, end_k);
Fluxmax_TFA = NaN(end_i, end_j, end_k);

% grab environment to pass to parallel workers
environment = getEnvironment();

parfor i = 1:end_i
    
    % set up environment for each worker
    restoreEnvironment(environment);
    changeCobraSolver('cplex_direct');
    
    % copy model
    model2 = TFA_model;

    % Set the substrate flux to -10
    model2 = changeTFArxnBounds (model2, Rxn(i), -10, 'l');
    
    for j = 1:end_j
        
        % Set the O2 flux
        model3 = changeTFArxnBounds (model2, 'DM_o2_e', O2(j), 'l');
        
        % Change the objective rxn biomass to 1
        objVar = 'F_Ec_biomass_iJO1366_WT_53p95M';
        objVarIndex = find_cell(objVar, model3.varNames);
        model3.f(objVarIndex) = 1;
        
        % Determine the maximum biomass production and set it as constant
        TFAsolution = solveTFAmodelCplex (model3);
        model3 = changeTFArxnBounds (model3, 'Ec_biomass_iJO1366_WT_53p95M', TFAsolution.val, 'b');
        
        % Change the objective rxn biomass to 0
        model3.f (objVarIndex) = 0;
        
        % Flux analysis
        for k = 1:end_k
            k_ = k+3478; % index of NF variable
            
            % Change the objective rxn to 1
            model3.f (k_) = 1;
            
            try
                % Daniel's approach
%                 % Set direction to minimise & solve
%                 model3.f (k_) = -1;
%                 TFAsolutionmin = solveTFAmodelCplex(model3);
%                 Fluxmin_TFA (i, j, k) = -TFAsolutionmin.val;
%                 
%                 % Set direction to maximise & solve
%                 model3.f (k) = 1;
%                 TFAsolutionmax = solveTFAmodelCplex(model3);
%                 Fluxmax_TFA (i, j, k) = TFAsolutionmax.val;

                % Andreas' approach
                % Set direction to minimise & solve
                model3.objtype = +1;
                TFAsolutionmin = solveTFAmodelCplex(model3);
                Fluxmin_TFA(i, j, k) = TFAsolutionmin.val;
                
                % Set direction to maximise & solve
                model3.objtype = -1;
                TFAsolutionmax = solveTFAmodelCplex(model3);
                Fluxmax_TFA(i, j, k) = TFAsolutionmax.val;
            catch
                % Fluxmin_TFA and Fluxmax_TFA are initialised as NaN, so no
                % need to re-assign here
            end
            
            % Change the objective rxn to 0
            model3.f (k) = 0;
        end
        
        % print an update
        fprintf('Finished TFA flux analysis for:\t%s\t%s\n',Rxn{i},O2_label{j});        
        
    end
    
    % Set the substrate flux to 0
    model2 = changeTFArxnBounds (model2, Rxn(i), 0, 'l');
    
end


%% Compare FBA and TFA fluxes

Comp_Struct.rxns = TFA_model.rxns;
Comp_Struct.number = NaN(length(TFA_model.rxns), length (O2), length (Rxn));

for i = 1:end_i
    for j = 1:end_j  
        for k = 1:end_k

            % Check if the TFA rxn is forward
            if Fluxmax_TFA(i,j,k) > 0 && Fluxmin_TFA (i,j,k) > 0

                % Check if the FBA rxn is forward
                if Fluxmax_FBA(i,j,k) > 0 && Fluxmin_FBA (i,j,k) > 0
                    Comp_Struct.number (k,j,i) = 0;
                else
                    Comp_Struct.number (k,j,i) = 1;
                end

            % Check if the TFA rxn is backward
            elseif Fluxmax_TFA(i,j,k) < 0 && Fluxmin_TFA (i,j,k) < 0

                % Check if the FBA rxn is backward
                if Fluxmax_FBA(i,j,k) < 0 && Fluxmin_FBA (i,j,k) < 0
                    Comp_Struct.number (k,j,i) = 0;
                else
                    Comp_Struct.number (k,j,i) = -1;
                end

            % Else the TFA rxn is bidirectional    
            else
                Comp_Struct.number (k,j,i) = 0;
            end
            
        end
    end
end


%% Export data to csv file

% Export glc model

Model_Glc_Aero = table (Comp_Struct.rxns, Comp_Struct.number(:,1,1));
writetable(Model_Glc_Aero, [pwd '/Ex2/Model_Glc_Aero.csv']);
Model_Glc_Anaero = table (Comp_Struct.rxns, Comp_Struct.number(:,2,1));
writetable(Model_Glc_Anaero, [pwd '/Ex2/Model_Glc_Anaero.csv']);

% Export lac model

Model_Lac_Aero = table (Comp_Struct.rxns, Comp_Struct.number(:,1,2));
writetable(Model_Lac_Aero, [pwd '/Ex2/Model_Lac_Aero.csv']);
Model_Lac_Anaero = table (Comp_Struct.rxns, Comp_Struct.number(:,2,2));
writetable(Model_Lac_Anaero, [pwd '/Ex2/Model_Lac_Anaero.csv']);

% Export ac model

Model_Ac_Aero = table (Comp_Struct.rxns, Comp_Struct.number(:,1,3));
writetable(Model_Ac_Aero, [pwd '/Ex2/Model_Ac_Aero.csv']);
Model_Ac_Anaero = table (Comp_Struct.rxns, Comp_Struct.number(:,2,3));
writetable(Model_Ac_Anaero, [pwd '/Ex2/Model_Ac_Anaero.csv']);

% Export et model

Model_Et_Aero = table (Comp_Struct.rxns, Comp_Struct.number(:,1,4));
writetable(Model_Et_Aero, [pwd '/Ex2/Model_Et_Aero.csv']);
Model_Et_Anaero = table (Comp_Struct.rxns, Comp_Struct.number(:,2,4));
writetable(Model_Et_Anaero, [pwd '/Ex2/Model_Et_Anaero.csv']);

%% short export

for i = 1:end_i
    for j = 1:end_j
        % create table of reactions
        T = table(Comp_Struct.rxns, Comp_Struct.number(:,j,i));
        % export all
        writetable(T, [pwd '/out/' Rxn{i} '_' O2_label{j} '.csv']);
        % export only ones that have changed
        writetable(T(T{:,2}~=0,:), [pwd '/out/' Rxn{i} '_' O2_label{j} '_changed.csv']);
    end
end