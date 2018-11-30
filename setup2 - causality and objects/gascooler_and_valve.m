clearvars
global bugnumber
bugnumber = 0;
% Initialization
givenVolume = 0.0192;
Parameters.InnerTubeDiameter = 0.0192;
Parameters.nParallelFlows = 20;
Parameters.OneTubelength = givenVolume/(Parameters.InnerTubeDiameter^2*pi/4);
Parameters.f1 = 0.01; % In case of just a few cells, sensitivity is pretty low
Parameters
gc = HeatExchanger;
p = [84.8e5 84.8e5-1.57e5];
h = [525e3 298e3];
tau = 0.1;
ODEoptions = [];
nCell = 5;
gc.initialize(nCell,p,h,Parameters,tau,ODEoptions)
Dm = 0.3186;
inputs.hInlet = 525e3;
DQ = -ones(nCell,1)*Dm*(inputs.hInlet-h(2))/nCell/Parameters.nParallelFlows;
inputs.DmInlet = Dm;
inputs.DmOutlet = Dm;
% Valve initialization
valve = Valve;
Kv = 0.9;
valve.initialize(Kv);
cR = 0.25; % capacity ratio
valveInputs.dInlet = gc.d(end);
valveInputs.pInlet = gc.p(end);
valveInputs.hInlet = gc.h(end);
valveInputs.pOutlet = 40e5;
% Gas Cooler input initialization
gcInputs.hInlet = h(1);
gcInputs.DmInlet = Dm;
% Controller intialization
PIvalve = PIController;
K = 1e-2;
Ti = 50;
mn = 0;
mx = 1;
neg = -1;
ref = 85e5;
timestep = 0.1;
initInt = 0;
PIvalve.initialize(K,Ti,initInt,mn,mx,neg,timestep);

% Simulation
itmax = 1000;
t_Dm = timestep:timestep:itmax*timestep;
Dm = zeros(gc.nCell+1,length(t_Dm));
tic
for it = 1:itmax
    % Controller
    cR = PIvalve.react(ref,gc.p(end));
    % Valve
    [valveDm, h] = valve.process(cR,valveInputs);
    gcInputs.DmOutlet = valveDm;
    % Gas cooler
    gcInputs.DQ = DQ;
    gc.timestep([(it-1)*timestep it*timestep],gcInputs);
    valveInputs.dInlet = gc.d(end);
    valveInputs.pInlet = gc.p(end);
    valveInputs.hInlet = gc.h(end);
    Dm(:,it) = gc.Dm;
end
toc
figure(1)
subplot(411)
plot(gc.record.t,gc.record.x(:,1:nCell))
subplot(412)
plot(gc.record.t,gc.record.x(:,1*nCell+1:2*nCell))
subplot(413)
plot(gc.record.t,gc.record.x(:,2*nCell+1:3*nCell))
subplot(414)
plot(t_Dm,Dm);
disp(['Number of CoolProp bugs were ' num2str(bugnumber)])