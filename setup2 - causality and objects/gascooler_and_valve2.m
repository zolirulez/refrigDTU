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
nCell = 2;
gc.initialize(nCell,p,h,Parameters,Tau,ODEoptions)
Dm = 0.3186;
% Valve initialization
valve = Valve;
Kv = 0.8;
valve.initialize(Kv,Tau,Dm);
Inputs.valveInputs.dInlet = gc.d(end);
Inputs.valveInputs.pInlet = gc.p(end);
Inputs.valveInputs.hInlet = gc.h(end);
Inputs.valveInputs.pOutlet = 40e5;
Inputs.valveInputs.capacityRatio = 0.01;
% Gas Cooler input initialization
Inputs.gcInputs.hInlet = h(1);
Inputs.gcInputs.hInlet = 525e3;
Inputs.gcInputs.DQ = -ones(nCell,1)*Dm*(Inputs.gcInputs.hInlet-h(2))/nCell/Parameters.nParallelFlows;
Inputs.gcInputs.DmInlet = Dm;
% Controller initialization
PIvalve = PIController;
K = 1e-5;
Ti = 50;
mn = 0;
mx = 1;
neg = -1;
ref = 85e5;
timestep = 1;
initInt = 0;
PIvalve.initialize(K,Ti,initInt,mn,mx,neg,timestep);

% Initialization of physical plant object
pp = PhysicalPlant;
Parts.gc = gc;
Parts.valve = valve;
Process = @process;
PostProcess = @postProcess; 
ODEoptions = [];
Initial = [gc.p; gc.h; gc.d; Dm];
pp.initialize(Parts,Process,PostProcess,Initial,ODEoptions);

% Simulation
itmax = 1000;
tic
for it = 1:itmax
    % Controller
    Inputs.valveInputs.capacityRatio = PIvalve.react(ref,pp.parts.gc.p(end));
    % Physical Plant
    pp.timestep(((it-1):it)*timestep,Inputs);
end
toc
figure(1)
subplot(411)
plot(pp.t,pp.parts.gc.record.x(:,1:nCell))
subplot(412)
plot(pp.t,pp.parts.gc.record.x(:,1*nCell+1:2*nCell))
subplot(413)
plot(pp.t,pp.parts.gc.record.x(:,2*nCell+1:3*nCell))
subplot(414)
plot(pp.t,pp.parts.valve.record.x);
disp(['Number of CoolProp bugs were ' num2str(bugnumber)])

% Functions for Physical Plant
function Dx = process(t,x,Inputs,pp)
    % States
    xGC = x(1:pp.parts.gc.nCell*3);
    xValve = x(end);
    % Valve
    DValve = pp.parts.valve.process(t,xValve,Inputs.valveInputs);
    Inputs.gcInputs.DmOutlet = pp.parts.valve.DmState;
    % Gas cooler
    DGC = pp.parts.gc.process(t,xGC,Inputs.gcInputs);
    Inputs.valveInputs.dInlet = pp.parts.gc.d(end);
    Inputs.valveInputs.pInlet = pp.parts.gc.p(end);
    Inputs.valveInputs.hInlet = pp.parts.gc.h(end);
    % Derivative
    Dx = [DGC; DValve];
end
function postProcess(X,pp)
    % States
    XGC = X(:,1:pp.parts.gc.nCell*3);
    XValve = X(:,end);
    % Recording
    pp.parts.gc.record.x = [pp.parts.gc.record.x; XGC];
    pp.parts.valve.record.x = [pp.parts.valve.record.x; XValve];
    pp.parts.gc.p = XGC(end,1:pp.parts.gc.nCell)';
    pp.parts.gc.h = XGC(end,pp.parts.gc.nCell+1:2*pp.parts.gc.nCell)';
    pp.parts.gc.d = XGC(end,2*pp.parts.gc.nCell+1:3*pp.parts.gc.nCell)';
end