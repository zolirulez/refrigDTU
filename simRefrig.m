clearvars
close all
%--------------------------------------------------------------
t0 =    0.0;        % [s] Initial time
tf = 1*60;         % [s] Final time
Ts = 1;            % [s] Sample Time
t = [t0:Ts:tf]';   % [s]Sampleinstants
N = length(t);
const_u = diag([25 7e-7]);
const_d = diag([273+35 273-15]);
u = const_u*ones(2,length(t));
d = const_d*ones(2,length(t));
%--------------------------------------------------------------
%SteadyState
%--------------------------------------------------------------
% Initialization
h = [450; 550; 300; 300]*1e3;
mdot = 0.05;
global memo
memo = [h; mdot];
x0 = [90e5; (h(2)+h(3))/2; 17e5; (h(1)+h(4))/2];

% Simulation
options = odeset('RelTol',1e-6,'AbsTol',1e-6);
nx = 10; nu = 2; 
X = x0;
T = t(1);
for k = 1:N-1
    k
    x = X(:,end);
    [Tk,Xk]=ode45(@Refrig_noRec,[t(k) t(k+1)],x,options,u(:,k),d(:,k));
    X = [X; Xk];
    T = [T; Tk];
end
k=N;