clc
clear all
addpath(genpath('C:\AAA\EPFL Y3\Principles and Application of Systems biology\MATLAB\matTFA-master'));
addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab'))
changeCobraSolver('cplex_direct');
load('C:\AAA\EPFL Y3\Principles and Application of Systems biology\MATLAB\ecoli_core_model.mat')

growthRates = zeros(21);
for i = 0:30
for j = 0:30
model = changeRxnBounds(model,'EX_glc(e)',-i,'b');
model = changeRxnBounds(model,'EX_o2(e)',-j,'b');
FBAsolution = optimizeCbModel(model,'max');
growthRates(i+1,j+1) = FBAsolution.f;
end
end

figure (2)
surfl([0:1:30],[0:1:30], growthRates)
xlabel('Glucose uptake rate (mmol.gDW^-^1.hr^-^1)')
ylabel('Oxygen uptake rate (mmol.gDW^-^1.hr^-^1)')
zlabel('Growth yield (hr^-^1)')
