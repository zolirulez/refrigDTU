function Dx = GasCoolerCells(t,x,u,param)
% Function to run the gas cooler

% Unwrap input arguments
p = x(1);
d = x(2);
Dm = u(1:2);
hDm = u(3:4);
DQ = u(5);
V = param(1);

% Mass balance
Dd = -V\diff(Dm);
% Energy balance and PH diagram
DdDp = CoolProp.PropsSI('d(D)/d(P)|H','D',d,'P',p,'CO2'); 
DdDh = CoolProp.PropsSI('d(D)/d(H)|P','D',d,'P',p,'CO2'); 
DpDh_vec = [-V d*V; DdDp DdDh]\[-diff(hDm)+DQ; Dd];

Dp = DpDh_vec(1);
Dx = [Dp; Dd];