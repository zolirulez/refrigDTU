function xdot = Refrig_noRec(t,x,u,d)

global memo

% Separating
pGC = x(1);
hGC = x(2);
pLT = x(3);
hLT = x(4);
% n = u(1);
% Ain = u(2);
Tamb = d(1);
Tcab = d(2);
inp = [u; x];

% Parameters
VGC = 1; % TODO
VLT = 1; % TODO
AGC = 10; % TODO
ALT = 10; % TODO
alphaGC = 6000; %TODO
alphaLT = 100; % TODO

% Iteration of equations
opt = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective'); 
LB = [200e3 200e3 200e3 200e3 0];
UB = [600e3 600e3 600e3 600e3 10];
param = lsqnonlin(@(param) RefrigParam(param,inp),memo,LB,UB,opt)
h = param(1:4);
mdot = param(5);

% Gas Cooler
dddp_GC = CoolProp.PropsSI('d(D)/d(P)|H','P',pGC,'H',hGC,'CO2');
dddh_GC = CoolProp.PropsSI('d(D)/d(H)|P','P',pGC,'H',hGC,'CO2');
TGC = CoolProp.PropsSI('T','P',pGC,'H',hGC,'CO2');
dGC = CoolProp.PropsSI('D','P',pGC,'H',hGC,'CO2');
MGC = dGC*VGC;
GCdot = [...
    dddp_GC dddh_GC; ...    
    VGC     -MGC]\...       
    [0; ...
    mdot*(h(2)-h(3))+alphaGC*AGC*(Tamb-TGC)...
    ];

% Evaporator, TODO - modify to saturation equations
dddp_LT = CoolProp.PropsSI('d(D)/d(P)|H','P',pLT,'H',hLT,'CO2');
dddh_LT = CoolProp.PropsSI('d(D)/d(H)|P','P',pLT,'H',hLT,'CO2');
TLT = CoolProp.PropsSI('T','P',pLT,'H',hLT,'CO2');
dLT = CoolProp.PropsSI('D','P',pLT,'H',hLT,'CO2');
MLT = dLT*VLT;
LTdot = [...
    dddp_LT dddh_LT; ...    
    VLT     -MLT]\...       
    [0; ...
    mdot*(h(3)-h(4))+alphaLT*ALT*(Tcab-TLT)...
    ];

% Memory
memo = param;
% Grouping derivatives
xdot = [GCdot; LTdot];