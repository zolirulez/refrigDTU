d = 250;
p = 85e5;
Dh_Dp = CoolProp.PropsSI('d(H)/d(P)|D','D',d,'P',p,'CO2')
v = 1/d;
cp = 0.85e3;
R = 8.3144598;
Tc = 304.25;
pc = 73.9e5;
a = 27*R*R*Tc*Tc/64/pc;
b = R*Tc/(8*pc);
M = 44.01e-3;
Vm = M*v;
(cp/R-1)*(Vm-b)+Vm

Dh_Dd = CoolProp.PropsSI('d(H)/d(D)|P','D',d,'P',p,'CO2')
a/(Vm*Vm)*(2-cp/R+2*b/Vm)+p*cp/R
a/Vm/Vm+(cp-R)/R*(p-a/Vm/Vm+2*b*a/Vm/Vm/Vm)+p