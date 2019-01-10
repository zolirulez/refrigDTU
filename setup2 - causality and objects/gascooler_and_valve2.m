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
Tau = 0.1;
ODEoptions = [];
nCell = 5;
gc.initialize(nCell,p,h,Parameters,Tau,ODEoptions)
Dm = 0.3186;
% Valve initialization
valve = Valve;
Kv = 0.9;
valve.initialize(Kv,Tau,Dm);
Inputs.valveInputs.dInlet = gc.d(end);
Inputs.valveInputs.pInlet = gc.p(end);
Inputs.valveInputs.hInlet = gc.h(end);
Inputs.valveInputs.pOutlet = 40e5;
Inputs.valveInputs.capacityRatio = 0.25;
% Gas Cooler input initialization
Inputs.gcInputs.hInlet = h(1);
Inputs.gcInputs.hInlet = 525e3;
DQ = -ones(nCell,1)*Dm*(Inputs.gcInputs.hInlet-h(2))/nCell/Parameters.nParallelFlows;
Inputs.gcInputs.DmInlet = Dm;
% Controller initialization
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

% Initialization of physical plant object
pp = PhysicalPlant;
Parts.gc = gc;
Parts.valve = valve;
Process = @process;
PostProcess = @postProcess; 
ODEoptions = [];
pp.initialize(Parts,Process,PostProcess,ODEoptions);

% Simulation
itmax = 1000;
tic
for it = 1:itmax
    % Controller
    cR = PIvalve.react(ref,gc.p(end));
    % Physical Plant
    pp.timestep(t(it-1:it),Inputs);
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

% Functions for Physical Plant
function Dx = process(pp,t,x,Inputs)
    % States
    xGC = x(1:pp.parts.gc.nCell*3);
    xValve = x(end);
    % Valve
    DValve = pp.parts.valve.process(t,xValve,Inputs.valveInputs);
    Inputs.gcInputs.DmOutlet = valveDm;
    % Gas cooler
    Inputs.gcInputs.DQ = DQ;
    DGC = pp.parts.gc.process(t,xGC,Inputs.gcInputs);
    Inputs.valveInputs.dInlet = pp.gc.d(end);
    Inputs.valveInputs.pInlet = pp.gc.p(end);
    Inputs.valveInputs.hInlet = pp.gc.h(end);
    % Derivative
    Dx = [DGC; DValve];
end
function postProcess(pp,X)
    pp.parts.gc.record.x = [pp.parts.gc.record.x(:,1:end-1); X];
    pp.parts.valve.record.x = [pp.parts.gc.record.Dm(:,end); X];
    pp.parts.gc.p = X(end,1:pp.gc.nCell)';
    pp.parts.gc.h = X(end,pp.gc.nCell+1:2*pp.gc.nCell)';
    pp.parts.gc.d = X(end,2*pp.gc.nCell+1:3*pp.gc.nCell)';
end