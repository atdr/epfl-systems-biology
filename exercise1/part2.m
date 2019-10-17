%% Exercise 1 - Part 2

clear, clc, close all
load('ecoli_core_model.mat')

EX_rxns = ["EX_glc(e)","EX_lac-D(e)","EX_gln-L(e)","EX_glu-L(e)","EX_fru(e)","EX_pyr(e)","EX_ac(e)","EX_etoh(e)"];%vector of exchange reactions
Gene=zeros(8,2);
O2 = [-20 0]; %aerobic or anaerobic reaction

for i = 1:length(EX_rxns)%loop for 8 reactions studied
    model = changeRxnBounds(model,EX_rxns(i),-10,'l'); %boundary setting for carbon source of interest
    for j = 1:2 %loop for aerobic/anaerobic reactions
        model = changeRxnBounds (model, 'EX_o2(e)', O2(j), 'l'); %boudray setting for aerobic/anaerobic reaction
        [grRatio,grRateKO,grRateWT,hasEffect,delRxns,fluxSolution] = singleGeneDeletion(model);
        Gene(i,j) = sum(grRatio <=0.1); %verify difference in growth rate to asses if essential gene
    end
    model = changeRxnBounds(model,EX_rxns(i),0,'l'); %reinitialise boundary setting of carbon source
end
