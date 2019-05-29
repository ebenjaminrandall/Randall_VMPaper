function [values,indices,ratios,alphaplot] = clinicalratios(data)

Tnew  = data.Tdata; 
Pnew  = data.Pdata; 
Hnew  = data.Hdata;
Pbar  = data.Pbar;
ts    = data.ts;
te    = data.te; 
tr    = data.tr; 
t4    = data.t4; 
dt    = data.dt; 

%% Extract quantities from SP data

%Find max SBP in phase I of VM
[SPmaxI,k_1] = max(Pnew(ts:ts+round(5/dt))); 
k_SPmaxI = ts + k_1 - 1; 

%Min SBP during phase II of VM (mmHg)
[SPminII,k_2] = min(Pnew(k_SPmaxI:ts+round(10/dt)));
k_SPminII = k_SPmaxI + k_2 - 1; 

%If SPmaxI < Pbar or SPminII > Pbar, find average of max and min points
if SPmaxI < Pbar || SPminII > Pbar
    k_SPendI = round((k_SPmaxI + k_SPminII)/2) - 1; 
else
    k_3 = find(Pnew(k_SPmaxI:k_SPminII) >= Pbar,1,'last');
    k_SPendI = k_SPmaxI + k_3 - 1; 
end 
SPendI = Pnew(k_SPendI); 

%Max SBP during late phase II of VM (mmHg)
[SPmaxII,k_4] = max(Pnew(k_SPminII:te-1)); 
k_SPmaxII = k_SPminII + k_4 -1; 

%Max SBP during phase IV of VM (mmHg)
[SPmaxIV,k_5] = max(Pnew(te:end)); 
k_SPmaxIV = te + k_5 - 1; 

%% Extract quantities from HR data

%Max HR during phase III of VM (ms)
k_7 = find(Hnew(te:tr) == max(Hnew(te:tr)),1,'first'); 
k_HmaxIII = te + k_7 - 1; 
HmaxIII = Hnew(k_HmaxIII); 
RRminIII = 60000/HmaxIII; 

%Min HR on interval before SBP returns to baseline 
[HminIV,k_8] = min(Hnew(k_HmaxIII:t4)); 
k_HminIV = k_HmaxIII + k_8 - 1; 
RRmaxIV = 60000/HminIV; 

%% Clinical ratios 

%VASOCONSTRICTIVE CAPACITY
%Formulated as the linear regression between the min BP in phase II and the
%max BP in late phase II
x_alpha = Tnew(k_SPminII:k_SPmaxII); 
y_alpha = Pnew(k_SPminII:k_SPmaxII);
if length(x_alpha) == 1
    alpha = 0;
    p_alpha = [];
    yfit_alpha = [];
else
    p_alpha = polyfit(x_alpha,y_alpha,1);
    yfit_alpha = polyval(p_alpha,x_alpha); 
    alpha = p_alpha(1); 
end 

%BARORECEPTOR SENSITIVITY
%Formulated as difference in RR-interval versus the difference in SBP
%overshoot 
beta = (RRmaxIV - RRminIII)/(SPmaxIV - Pbar); 

%VALSALVA RATIO
%Ratio between maximum HR in phase III to minimum HR in phase IV 
gamma = HmaxIII/HminIV; 

%% Outputs 

values  = [SPmaxI; SPendI; SPminII; 
    SPmaxII; SPmaxIV; 
    HmaxIII; HminIV]; 
indices = [k_SPmaxI; k_SPendI; k_SPminII; 
    k_SPmaxII; k_SPmaxIV; 
    k_HmaxIII; k_HminIV]; 
ratios = [alpha; beta; gamma];
ratios = round(ratios,1); 

alphaplot.x_alpha = x_alpha;
alphaplot.yfit_alpha = yfit_alpha;
alphaplot.p_alpha = p_alpha; 

