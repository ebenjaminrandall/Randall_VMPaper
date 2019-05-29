function [pars, low, hi] = load_global(data)

age    = data.age;  

%Initial mean values
Pbar   = data.Pbar;
HminR  = data.HminR; 
HmaxR  = data.HmaxR; 
Hbar   = data.Hbar; 

%% PARAMETERS  

A    = 5;           %1  %dimensionless 
B    = .5;            %2  %dimensionless

Kb   = .1;            %3  %dimensionless
Kpb  = 5;
Kpr  = 1; 
Ks   = 5; 

taub  = .9;            %6  %sec
taupb = 1.8;           %7  %sec
taupr = 6;
taus  = 10;            %8  %sec
tauH  = .5; 

qw   = .04;           %10 %mmHg^(-1)
qpb  = 10;           %11 %sec
qpr  = 1; 
qs   = 10;           %12 %sec

Ds   = 3;             %19 sec

%% Patient specific parameters

sw  = Pbar;            %13 %mmHg
spr = data.Pthbar;  

%Intrinsic HR
HI = 118 - .57*age;   %16 %bpm
if HI < Hbar 
    HI = Hbar;
end 
%Maximal HR
HM = 208 - .7*age;    %bpm     
Hs = (1/Ks)*(HM/HI - 1); 

%% Calculate sigmoid shifts

Pc_ss  = data.Pbar; 
ewc_ss = 1 - sqrt((1 + exp(-qw*(Pc_ss - sw)))/(A + exp(-qw*(Pc_ss - sw)))); 
ebc_ss = Kb*ewc_ss; 
ec_ss  = ewc_ss - ebc_ss; 

Pa_ss  = data.Pbar - data.Pthbar; 
ewa_ss = 1 - sqrt((1 + exp(-qw*(Pa_ss - sw)))/(A + exp(-qw*(Pa_ss - sw)))); 
eba_ss = Kb*ewa_ss; 
ea_ss  = ewa_ss - eba_ss; 

n_ss   = B*ec_ss + (1 - B)*ea_ss;

Tpb_ss = .8;
Ts_ss  = .2; 

%Steady-state sigmoid shifts 
spb = n_ss + log(Kpb/Tpb_ss - 1)/qpb;  %14 sec^(-1)
ss  = n_ss -  log(Ks/Ts_ss - 1)/qs;    %15 sec^(-1)

%% At end of expiration and inspiration

Gpr_ss = 1/(1 + exp(qpr*(data.Pthbar - spr)));

Tpr_ss = Kpr*Gpr_ss; 
Hpr = (HmaxR - HminR)/HI/Tpr_ss ;

Hpb = (1 - Hbar/HI + Hpr*Tpr_ss + Hs*Ts_ss)/Tpb_ss;

%% Outputs

pars = [A; B;              
    Kb; Kpb; Kpr; Ks;                           %Gains
    taub; taupb; taupr; taus; tauH;     %Time Constants
    qw; qpb; qpr; qs;                %Sigmoid Steepnesses
    sw; spb; spr; ss;               %Sigmoid Shifts
    HI; Hpb; Hpr; Hs;
    Ds];                 %Heart Rate Parameters 
    
pars = log(pars);

low     = pars - log(10);
hi      = pars + log(10);
low([2 21 22 23])  = log(.01); 
hi([2 21 22 23])   = log(1); 

%% Load in optimized parameter values

load optHR.mat pars_LM
pars = pars_LM; 

end 

