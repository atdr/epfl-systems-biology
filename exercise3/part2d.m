%% Exercise 3 - Part 2d
clear, clc, close all

% nomenclature of chemical species
% w(1) MAPKKK
% w(2) MAPKKK-P
% w(3) MAPKK
% w(4) MAPKK-P
% w(5) MAPKK-PP
% w(6) MAPK
% w(7) MAPK-P
% w(8) MAPK-PP

% selected time scale to reach steady-state
timerange = [0:4e3]; % s

% data from brief and initial conditions
w0 = [100 0 300 0 0 300 0 0]; % nM

Vmax1 = [2.5,1.0,0.5,0.2,0.1];
for i = 1:length(Vmax1) % nM/s
    [t,w] = ode45(@(t,w)balances(t,w,Vmax1(i)),timerange,w0); % ODE solver
    figure (i)
    plot(t,w(:,1),'b-',t,w(:,2),'c-',t,w(:,3),'y-',t,w(:,4),'m-',t,w(:,5),'r-',t,w(:,6),'k-',t,w(:,7),'g-',t,w(:,8))
    xlabel('time / s')
    ylabel('concentration / nM')
    %title('Concentrations as a function of time')
    if i == 1
        legend('[MAPKKK]','[MAPKKK-P]','[MAPKK]','[MAPKK-P]','[MAPKK-PP]','[MAPK]','[MAPK-P]','[MAPK-PP]')
        figExport(15,10,sprintf('part2D_%d_08112019',10*Vmax1(i)))
    else
        figExport(7,5,sprintf('part2D_%d_08112019',10*Vmax1(i)))
    end
end

function [dwdt] = balances(t,w,Vmax1)
% define given values from the brief
Vmax = [Vmax1 0.25 0 0 0.75 0.75 0 0 0.5 0.5]; % nM/s
kcat = 0.025; % s^-1
Km = [10 8 15]; % nM
Ka = 0;
Ki = 0.1;

% define rates of reaction
v1 = Vmax(1)*w(1)*(1+Ka*w(8))/((Km(1)+w(1))*(1+Ki*w(8)));
v2 = Vmax(2)*w(2)/(Km(2)+w(2));
v3 = kcat*w(2)*w(3)/(Km(3)+w(3));
v4 = kcat*w(2)*w(4)/(Km(3)+w(4));
v5 = Vmax(5)*w(5)/(Km(3)+w(5));
v6 = Vmax(6)*w(4)/(Km(3)+w(4));
v7 = kcat*w(5)*w(6)/(Km(3)+w(6));
v8 = kcat*w(5)*w(7)/(Km(3)+w(7));
v9 = Vmax(9)*w(8)/(Km(3)+w(8));
v10 = Vmax(10)*w(7)/(Km(3)+w(7));

% define differential equations for material balances
dwdt(1) = v2-v1;
dwdt(2) = v1-v2;
dwdt(3) = v6-v3;
dwdt(4) = v3-v4+v5-v6;
dwdt(5) = v4-v5;
dwdt(6) = v10-v7;
dwdt(7) = v7+v9-v10-v8;
dwdt(8) = v8-v9;
dwdt = dwdt(:); % force dwdt to be a column vector to use ode15s
end
