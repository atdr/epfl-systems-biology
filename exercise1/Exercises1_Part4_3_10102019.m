clear, clc, close all
load('ecoli_core_model.mat')

sol = struct;

% flux limit in each dimension
upperFlux = 30;

% list of substrates to consider
subs = {'EX_glc(e)'; 'EX_pyr(e)'; 'EX_succ(e)'};
subsName = {'glucose'; 'pyruvate'; 'succinate';};
% subs = {'EX_succ(e)'};
% subsName = {'succinate';};

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
        model = changeRxnBounds(model,sub,-subFlux,'l');
        
        % change the oxygen flux
        for oxFlux = 0:upperFlux
            model = changeRxnBounds(model,'EX_o2(e)',-oxFlux,'l');
            
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
    
    % now we're done with this substrate, reset its flux to 0
    model = changeRxnBounds(model,sub,0,'l');
end

%% create plots

% select points (subFlux,oxFlux) to export data for flux profile
% glucose
charPts{1} = [10, 1; 10, 3; 10, 6; 10, 10; 10, 20];
% pyruvate
charPts{2} = [15, 10; 15, 25];
% succinate
charPts{3} = [0, 9; 0, 24];

% go through each substrate
for subIndex = 1:numel(subs)
    sub = subs{subIndex};
    subName = subsName{subIndex};
    
    % plot growth yield as a function of substrate uptake (for oxFlux = 20)
    sz = size(sol(subIndex).growthRates);
    figure
    plot([0:sz(1)-1], sol(subIndex).growthRates(:,20+1))
    xlabel('uptake rate')
    ylabel('growth yield / h^{-1}')
    figExport(5,5,sprintf('%s-yield-ox20',subName))
    
    % plot phenotypic phase plane (3D)
    figure
    surfl([0:sz(2)-1],[0:sz(1)-1], sol(subIndex).growthRates)
    xlabel('oxygen')
    ylabel(subName)
    zlabel('growth yield / h^{-1}')
    figExport(8,7,sprintf('%s-3D',subName))
    
    % plot shadow price regions
    figure
    pcolor([0:sz(2)-1],[0:sz(1)-1], sol(subIndex).shadowPrices)
    xlabel('oxygen uptake rate')
    ylabel(sprintf('%s uptake rate',subName))
    figExport(6,6,sprintf('%s-shadow-price-regions',subName))
    
    % go through each characteristic point for this substrate
    for charPtIndex = 1:size(charPts{subIndex},1)
        charPtFluxes = charPts{subIndex}(charPtIndex,:);
        
        % export flux data
        writetable(table(model.rxns,sol(subIndex).fluxes{charPtFluxes(1)+1,charPtFluxes(2)+1}),  [pwd '/out/' sprintf('%s_sub%s_ox%s_flux_data.csv', subName, num2str(charPtFluxes(1)), num2str(charPtFluxes(2)))] );
    end
end
