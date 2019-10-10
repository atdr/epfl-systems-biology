clc
clear all
addpath(genpath('C:\AAA\EPFL Y3\Principles and Application of Systems biology\MATLAB\matTFA-master'));
addpath(genpath('C:\Program Files\IBM\ILOG\CPLEX_Studio129\cplex\matlab'))
changeCobraSolver('cplex_direct');
load('C:\AAA\EPFL Y3\Principles and Application of Systems biology\MATLAB\ecoli_core_model.mat')

model = changeRxnBounds (model, 'EX_o2(e)', -20, 'b'); %boudary setting for aerobic/anaerobic reaction
growthRates = zeros(21,1)
for i = 0:20 %loop for 4 exchange reactions
    model = changeRxnBounds (model, 'EX_glc(e)', -i, 'b'); %boundary setting for carbon source of interest
    FBAsolution = optimizeCbModel (model, 'max'); 
    growthRates(i+1) = FBAsolution.f;
end

figure (1)
plot([-20:1:0], growthRates)
xlabel('Glucose uptake rate (mmol.gDW^-^1.hr^-^1)')
ylabel('Growth yield (hr^-^1)')
