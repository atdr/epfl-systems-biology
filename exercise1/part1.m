%% Exercise 1 - Part 1

clear, clc, close all
load('ecoli_core_model.mat')

model = changeObjective(model,'Biomass_Ecoli_core_N(w/GAM)-Nmet2');%define objective function

O2 = [-20 0]; %aerobic or anaerobic reaction
EX_rxns = ["EX_glc(e)","EX_lac-D(e)","EX_fru(e)","EX_etoh(e)"];%vector of exchange reactions
for i = 1:4 %loop for 4 exchange reactions
    model = changeRxnBounds (model, EX_rxns(i), -10, 'l'); %boundary setting for carbon source of interest
    for j = 1:2 %loop for aerobic/anaerobic reactions
        model = changeRxnBounds (model, 'EX_o2(e)', O2(j), 'l'); %boudary setting for aerobic/anaerobic reaction
        FBAsolution = optimizeCbModel (model, 'max');
        Yield (i,j) = FBAsolution.f;
    end
    model = changeRxnBounds (model, EX_rxns(i), 0, 'l');%reinitialise boundary setting of carbon source
end
