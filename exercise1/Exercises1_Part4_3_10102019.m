clear, clc, close all
load('ecoli_core_model.mat')

charPts = [1,3,6,10,20];
sol = struct;

% flux limit in each dimension
upperFlux = 30;

% list of substrates to consider
subs = {'EX_glc(e)'; 'EX_pyr(e)'; 'EX_succ(e)'};

for subIndex = 1:numel(subs)
    sub = subs(subIndex);
    sol(subIndex).label = sub;
    sol(subIndex).growthRates = zeros(upperFlux);
    sol(subIndex).shadowPrices = zeros(upperFlux);

% substrate flux
for subFlux = 0:upperFlux
    model = changeRxnBounds(model,sub,-subFlux,'b');
    
    % oxygen flux
    for oxFlux = 0:upperFlux
        model = changeRxnBounds(model,'EX_o2(e)',-oxFlux,'b');
        
        % do the optimistion
        FBAsolution = optimizeCbModel(model,'max');
        
        % save the growth rate
        sol(subIndex).growthRates(subFlux+1,oxFlux+1) = FBAsolution.f;
        
        % save the shadow price (catch errors)
        try
            sol(subIndex).shadowPrices(subFlux+1,oxFlux+1) = FBAsolution.w(28);
        catch
            sol(subIndex).shadowPrices(subFlux+1,oxFlux+1) = NaN;
        end
        
    end
end

end

% if (i == 10 & ismember(j,charPts))
%     writetable(table(model.rxns,FBAsolution.x), ['ptype_glucose_10_' num2str(j) '_flux_data.csv']);
% end


figure (3)
pcolor([0:1:30],[0:1:30], shadowPrices)
xlabel('Glucose uptake rate (mmol.gDW^-^1.hr^-^1)')
ylabel('Oxygen uptake rate (mmol.gDW^-^1.hr^-^1)')
