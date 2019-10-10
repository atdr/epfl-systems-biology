clc
clear all
addpath(genpath('C:\AAA\EPFL Y3\Principles and Application of Systems biology\MATLAB\matTFA-master'));
addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab'))
changeCobraSolver('cplex_direct');
load('C:\AAA\EPFL Y3\Principles and Application of Systems biology\MATLAB\ecoli_core_model.mat')

charPts = [1,3,6,10,20];

growthRates = zeros(21);
for
for i = 0:30
for j = 0:30
model = changeRxnBounds(model,'EX_glc(e)',-i,'b');
model = changeRxnBounds(model,'EX_o2(e)',-j,'b');
FBAsolution = optimizeCbModel(model,'max');
growthRates(i+1,j+1) = FBAsolution.f;
try
   shadowPrices(i+1,j+1) = FBAsolution.w(28);
catch
   shadowPrices(i+1,j+1) = NaN;
end

if (i == 10 & ismember(j,charPts))
    writetable(table(model.rxns,FBAsolution.x), ['ptype_glucose_10_' num2str(j) '_flux_data.csv']);
end
end
end
end

figure (3)
pcolor([0:1:30],[0:1:30], shadowPrices)
xlabel('Glucose uptake rate (mmol.gDW^-^1.hr^-^1)')
ylabel('Oxygen uptake rate (mmol.gDW^-^1.hr^-^1)')
