function resid = RefrigParam(param,inp)
% This function is to decrease the residuals in the four main equations
% of the refrigeration system.

% Separating
h = param(1:4);
mdot = param(5);
u = inp(1:2);
x = inp(3:6);
n = u(1);
Ain = u(2);
pGC = x(1);
hGC = x(2);
pLT = x(3);

% Initialization of residual
resid = zeros(5,1);

% Constant parameters
VComp = 5e-5;
kappa = 1; % TODO
f = 2;
eta = 0.7;

% High pressure valve, Bernoulli
d3 = CoolProp.PropsSI('D','P',pGC,'H',h(3),'CO2');
resid(1) = mdot - Ain*sqrt(f/kappa*(pGC-pLT)*d3);

% Compressor, flow and isentropic efficiency
d1 = CoolProp.PropsSI('H','P',pLT,'H',h(1),'CO2');
resid(2) = d1 - mdot/(VComp*n);
s1 = CoolProp.PropsSI('S','P',pLT,'H',h(1),'CO2');
h2id = CoolProp.PropsSI('H','P',pGC,'S',s1,'CO2');
resid(3) = h(2) - (h2id-h(1))/eta+h(1);

% Gas Cooler, linear enthalpy distribution
resid(4) = h(3) - 2*hGC-h(2);

% Valve, isenthalpy
resid(5) = h(4) - h(3);