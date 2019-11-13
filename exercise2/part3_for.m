clear all
% Change the solver 
changeCobraSolver('cplex_direct');

%% Load the model
tmp = load('small_ecoli.mat');
model = tmp.model_red;
clear tmp

%% Data for the Flux Analysis

O2 = [-20 0];
Rxn = {'DM_glc_e','DM_lac-D_e','DM_ac_e','DM_etoh_e' };

%% FBA Flux Analysis

% Copy the model
FBA_model = model;

% Save the flux data
Fluxmin_FBA = NaN (length(Rxn), length (O2), length(FBA_model.rxns));
Fluxmax_FBA = NaN (length(Rxn), length (O2), length(FBA_model.rxns));

for i = 1:length(Rxn)

    % Set the substrate flux
    FBA_model = changeRxnBounds (FBA_model, Rxn(i), -10, 'l');
    
    for j = 1:length(O2)
        
        % Set the O2 flux
        model2 = changeRxnBounds (FBA_model, 'DM_o2_e', O2(j), 'l');
        
        % Determine the maximum biomass production and set it as constant
        model2 = changeObjective (model2, 'Ec_biomass_iJO1366_WT_53p95M');
        FBAsolution = optimizeCbModel (model2, 'max');
        model2 = changeRxnBounds (model2, 'Ec_biomass_iJO1366_WT_53p95M', FBAsolution.f, 'b');
        
        % Flux analysis
        for k = 1:length(model2.rxns)
            
            % Set the objective reaction
            model2 = changeObjective (model2, model2.rxns (k));
            
            % Perform the flux analysis
            FBAsolutionmin = optimizeCbModel (model2, 'min');
            Fluxmin_FBA (i, j, k) = FBAsolutionmin.f;
            FBAsolutionmax = optimizeCbModel (model2, 'max');
            Fluxmax_FBA (i, j, k) = FBAsolutionmax.f;
        end
    end
    
    FBA_model = changeRxnBounds (FBA_model, Rxn(i), 0, 'l');

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

% Copy model
model2 = TFA_model;

% Save the flux data
Fluxmin_TFA = NaN (length(Rxn), length (O2), 599);
Fluxmax_TFA = NaN (length(Rxn), length (O2), 599);

for i = 1:length(Rxn)
    
    % Set the substrate flux to -10
    model2 = changeTFArxnBounds (model2, Rxn(i), -10, 'l');
    
    for j = 1:length(O2)
        
        % Set the O2 flux
        model3 = changeTFArxnBounds (model2, 'DM_o2_e', O2(j), 'l');
        
        % Change the objective rxn biomass to 1
        objVar = 'F_Ec_biomass_iJO1366_WT_53p95M';
        objVarIndex = find_cell(objVar, model3.varNames);
        model3.f(objVarIndex) = 1;
        
        % Determine the maximum biomass production and set it as constant
        try
            TFAsolution = solveTFAmodelCplex (model3);
            model3 = changeTFArxnBounds (model3, 'Ec_biomass_iJO1366_WT_53p95M', TFAsolution.val, 'b');
            
            % Change the objective rxn biomass to 0
            model3.f (objVarIndex) = 0;
            
            % Flux analysis
            for k = 3479:4077
                
                % Change the objective rxn to 1
                model3.f (k) = 1;
                
                % Set direction to minimise & solve
                model3.objtype = 1;
                TFAsolutionmin = solveTFAmodelCplex(model3);
                Fluxmin_TFA (i, j, (k-3478)) = TFAsolutionmin.val;
                
                % Set direction to maximise & solve
                model3.objtype = -1;
                TFAsolutionmax = solveTFAmodelCplex(model3);
                Fluxmax_TFA (i, j, (k-3478)) = TFAsolutionmax.val;

                % Change the objective rxn to 0
                model3.f (k) = 0;
            end
            
        catch
            Fluxmin_TFA (i, j, :) = NaN;
            Fluxmax_TFA (i, j, :) = NaN;
        end
    end
    
    % Set the substrate flux to 0
    model2 = changeTFArxnBounds (model2, Rxn(i), 0, 'l');

end


%% Compare FBA and TFA fluxes

Comp_Struct.rxns = TFA_model.rxns;
Comp_Struct.number = NaN(length(TFA_model.rxns), length (O2), length (Rxn));

for i = 1:length(Rxn)
    
    for j = 1:length (O2)
        
        for k = 1: length (TFA_model.rxns)

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

            % Else if the TFA rxn is bidirectional    
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

