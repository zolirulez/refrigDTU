clearvars

% Simulation constraints
t0 =    0;          % [s] Initial time
tf = 60*60;         % [s] Final time
Ts = 60;            % [s] Sample Time
t = [t0:Ts:tf]';    % [s] Sample instants
N = length(t);

% Constants
V = 0.0192;
param = V;

% Inputs
h = [545; 300]*1e3;
Dm = [0.05; 0.05];
hDm = h.*Dm;
DQ = diff(hDm); % Chosen for steady state
inputs = [Dm; hDm; DQ];
u = diag(inputs)*ones(length(inputs),N);

% Initial conditions
p0 = 90e5;
d0 = CoolProp.PropsSI('D','H',mean(h),'P',p0,'CO2'); 
x0 = [p0; d0];
disp('The initial enthalpy value is')
mean(h)

% Simulation
X = zeros(length(x0),0);
T = zeros(length(x0),0);
x = x0;
for k = 1:N-1
    % Simulate from time t[k] to time t[k+1]
    [Tk,Xk]=ode15s(@GasCoolerCells,[t(k) t(k+1)],x(:,k),[],u(:,k),param);
    x(:,k+1) = Xk(end,:)';
    
    % Store the simulated results (for plotting)
    T = [T; Tk];
    X = [X; Xk];
end

% Plotting
p = x(1,:);
d = x(2,:);
figure(1)
subplot(211)
scatter(t/60,p)
title('Pressure')
subplot(212)
scatter(t/60,d)
title('Density')

% Final enthalpy value
disp('The final enthalpy value is')
CoolProp.PropsSI('H','P',p(end),'D',d(end),'CO2')