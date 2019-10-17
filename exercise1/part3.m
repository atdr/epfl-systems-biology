%% Exercise 1 - Part 3

clear, clc, close all
load('ecoli_core_model.mat')

model = changeRxnBounds (model, "EX_glc(e)", -10, 'l'); % boundary setting for carbon source

% initialise output vectors
minFlux = zeros (5,2);
maxFlux = zeros (5,2);


optPercentage = 80;
O2 = [-20 0]; % aerobic or anaerobic
EX_rxns = {"EX_ac(e)","EX_co2(e)","EX_lac-D(e)","ACONTb","MDH"}; % vector of reactions

for j = 1:2 % loop for aerobic/anaerobic reactions
    model = changeRxnBounds (model, 'EX_o2(e)', O2(j), 'l'); % boundary setting for aerobic/anaerobic reaction
    [minFlux(:,j),maxFlux(:,j),Vmin,Vmax] = fluxVariability(model,optPercentage,[],EX_rxns);
end
