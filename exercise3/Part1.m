clear,clc,close all

%define constants
Vmax= 50; %mmol/min
KD = 0.25; %mM
KI = 1; %mM

%define values for [I] and [S]
I = [0, 1, 2, 3]; %mM
S = (0:0.01:20); %mM

for i=1:4 %for each [I]
    for k=1:length(S) %for each [S]
    v(i,k) = (Vmax*S(k))/(KD+S(k)+(S(k)*I(i)/KI)+(I(i)/(KI*KD)));
    end
end 

figure (1)
plot(S(:),v(1,:),S,v(2,:),S,v(3,:),S,v(4,:))
xlabel('[S] (mM)')
ylabel('v (mmol/min)')
legend('[I]=0 mM','[I]=1 mM','[I]=2 mM','[I]=3 mM')

figExport(12,8,'part1_07112019')