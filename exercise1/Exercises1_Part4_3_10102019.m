clear, clc, close all
load('ecoli_core_model.mat')

sol = struct;

% flux limit in each dimension
upperFlux = 30;

% list of substrates to consider
subs = {'EX_glc(e)'; 'EX_pyr(e)'; 'EX_succ(e)'};
subsName = {'glucose'; 'pyruvate'; 'succinate';};

%% optimisation process

% go through each substrate
for subIndex = 1:numel(subs)
    sub = subs{subIndex};
    
    % initialise array layer
    sol(subIndex).label = sub;
    sol(subIndex).growthRates = zeros(upperFlux);
    sol(subIndex).shadowPrices = zeros(upperFlux);
    sol(subIndex).fluxes = cell(upperFlux);
    
    % change the substrate flux
    for subFlux = 0:upperFlux
        model = changeRxnBounds(model,sub,-subFlux,'b');
        
        % change the oxygen flux
        for oxFlux = 0:upperFlux
            model = changeRxnBounds(model,'EX_o2(e)',-oxFlux,'b');
            
            % do the optimistion
            FBAsolution = optimizeCbModel(model,'max');
            
            % save the growth rate
            sol(subIndex).growthRates(subFlux+1,oxFlux+1) = FBAsolution.f;
            
            % save the shadow price
            try
                sol(subIndex).shadowPrices(subFlux+1,oxFlux+1) = FBAsolution.w(28);
            catch
                sol(subIndex).shadowPrices(subFlux+1,oxFlux+1) = NaN;
            end
            
            % save the fluxes
            try
                sol(subIndex).fluxes{subFlux+1,oxFlux+1} = FBAsolution.x;
            catch
                sol(subIndex).fluxes(subFlux+1,oxFlux+1) = cell(1);
            end
            
        end
    end
end

%% create plots

charPts = [1,3,6,10,20];
% if (i == 10 & ismember(j,charPts))
%     writetable(table(model.rxns,FBAsolution.x), ['ptype_glucose_10_' num2str(j) '_flux_data.csv']);
% end

% go through each substrate
for subIndex = 1:numel(subs)
    sub = subs{subIndex};
    subName = subsName{subIndex};
    
    % plot growth yield as a function of substrate uptake (for oxFlux = 20)
    sz = size(sol(subIndex).growthRates);
    figure
    plot([0:sz(1)-1], sol(subIndex).growthRates(:,20+1))
    xlabel(sprintf('%s uptake rate / mmol.g_{DW}^{-1}.h^{-1}',subName))
    ylabel('growth yield / h^{-1}')
    figExport(8,8,sprintf('%s-yield-ox20',subName))
    
    % plot phenotypic phase plane (3D)
    figure
    surfl([0:sz(2)-1],[0:sz(1)-1], sol(subIndex).growthRates)
    xlabel('oxygen uptake rate / mmol.g_{DW}^{-1}.h^{-1}')
    ylabel(sprintf('%s uptake rate / mmol.g_{DW}^{-1}.h^{-1}',subName))
    zlabel('growth yield / h^{-1}')
    figExport(16,16,sprintf('%s-3D',subName))
    
    % plot shadow price regions
    figure
    pcolor([0:sz(2)-1],[0:sz(1)-1], sol(subIndex).shadowPrices)
    xlabel('oxygen uptake rate / mmol.g_{DW}^{-1}.h^{-1}')
    ylabel(sprintf('%s uptake rate / mmol.g_{DW}^{-1}.h^{-1}',subName))
    figExport(8,8,sprintf('%s-shadow-price-regions',subName))
    
end
