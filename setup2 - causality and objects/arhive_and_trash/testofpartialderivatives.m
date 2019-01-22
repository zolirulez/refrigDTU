d = 220;
p = 85e5;
Dh_Dp = CoolProp.PropsSI('d(H)/d(P)|D','D',d,'P',p,'CO2')
v = 1/d;
cp = 0.85e3;
R = 188.92; %8.3144598;
Tc = 304.25;
pc = 73.9e5;
a = 27*R*R*Tc*Tc/64/pc;
b = R*Tc/(8*pc);
% M = 44.01e-3;
% Vm = M*v;
kappa = cp/R-1;
v*(kappa+1)-kappa*b

Dh_Dd = CoolProp.PropsSI('d(H)/d(D)|P','D',d,'P',p,'CO2')
-(p*(kappa+1)+a/v/v*kappa*(2*b/v-1))/d/d
-p*(kappa+1)/d/d-a*kappa*(2*d*b-1)

% Note: for low pressures, and under the saturation curve it is not valid.
% It seems that there is some problem with the dhdv