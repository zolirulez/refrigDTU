clearvars
global bugnumber 
bugnumber = 0;
% Initialization
givenVolume = 0.0192;
Parameters.InnerTubeDiameter = 0.007;
Parameters.nParallelFlows = 10;
Parameters.OneTubelength = givenVolume/(Parameters.InnerTubeDiameter^2*pi/4);
Parameters.f1 = 0.001; % In case of just a few cells, sensitivity is pretty low
gc = HeatExchanger;
p = [84.8e5 84.8e5-1.57e5];
h = [450e3 298e3];
Tau = 0.1;
ODEoptions = [];
nCell = 2;
gc.initialize(nCell,p,h,Parameters,ODEoptions)
Dm = 0.321;
% HP Valve initialization
HPValve = Valve;
Kv = 0.8;
HPValve.initialize(Kv,Tau,Dm);
Inputs.HPValveInputs.dInlet = gc.d(end);
Inputs.HPValveInputs.pInlet = gc.p(end);
Inputs.HPValveInputs.hInlet = gc.h(end);
Inputs.HPValveInputs.pOutlet = 38e5;
Inputs.HPValveInputs.capacityRatio = 0;
HPValve.enthalpy(Inputs.HPValveInputs.hInlet);
% Gas Cooler input initialization
Inputs.gcInputs.hInlet = h(1);
Inputs.gcInputs.hInlet = 525e3;
Inputs.gcInputs.DQ = -ones(nCell,1)*Dm*(Inputs.gcInputs.hInlet-h(2))/nCell/Parameters.nParallelFlows;
Inputs.gcInputs.DmInlet = Dm;
% Receiver initialization
recVolume = 0.133;
rec = Receiver;
pRec = 35e5;
hRec = 250e3;
DmGas = 0.126;
rec.initialize(pRec,hRec,recVolume,ODEoptions);
rec.separation(); % For the receiver valve
Inputs.recInputs.hInlet = gc.h(end);
Inputs.recInputs.DmInlet = Dm;
Inputs.recInputs.DmGas = DmGas;
Inputs.recInputs.DmLiquid = Inputs.recInputs.DmInlet - Inputs.recInputs.DmGas;
% Receiver Valve initialization
recValve = Valve;
Kv = 2;
recValve.initialize(Kv,Tau,DmGas);
Inputs.recValveInputs.dInlet = rec.dGas;
Inputs.recValveInputs.pInlet = rec.p;
Inputs.recValveInputs.hInlet = rec.hGas;
Inputs.recValveInputs.pOutlet = 30e5;
Inputs.recValveInputs.capacityRatio = 0;
% Controller initialization
PIHPValve = PIController;
PIrecValve = PIController;
K = 1e-6;
Ti = 50;
mn = 0;
mx = 1;
neg = -1;
refHP = 85e5;
refRec = 38e5;
timestep = 5;
initInt = 0;
PIHPValve.initialize(K,Ti,initInt,mn,mx,neg,timestep);
K = 10^-6.5;
Ti = 50;
PIrecValve.initialize(K,Ti,initInt,mn,mx,neg,timestep);



% Initialization of physical plant object
ODEoptions.jancsi = 'jancsi';
pp = PhysicalPlant;
Parts.gc = gc;
Parts.HPValve = HPValve;
Parts.rec = rec;
Parts.recValve = recValve;
Process = @process;
PreProcess = @preProcess; 
PostProcess = @postProcess; 
ODEoptions = [];
Initial = [gc.p; gc.h; gc.d; HPValve.DmState; rec.p; rec.h; rec.d; recValve.DmState];
pp.initialize(Parts,Process,PostProcess,Initial,ODEoptions);
pp.initialInputs(Inputs);

% Simulation
itmax = 100;
tic
for it = 1:itmax
    % Controller
    pp.Inputs.recValveInputs.capacityRatio = PIrecValve.react(refRec,pp.parts.rec.p);
    pp.Inputs.HPValveInputs.capacityRatio =  PIHPValve.react(refHP,pp.parts.gc.p(end));
    % Physical Plant
    pp.timestep(((it-1):it)*timestep);
end
toc
figure(1)
subplot(811)
plot(pp.t,pp.parts.gc.record.x(:,1:nCell))
subplot(812)
plot(pp.t,pp.parts.gc.record.x(:,1*nCell+1:2*nCell))
subplot(813)
plot(pp.t,pp.parts.gc.record.x(:,2*nCell+1:3*nCell))
subplot(814)
plot(pp.t,pp.parts.HPValve.record.x);
subplot(815)
plot(pp.t,pp.parts.rec.record.x(:,1))
subplot(816)
plot(pp.t,pp.parts.rec.record.x(:,2))
subplot(817)
plot(pp.t,pp.parts.rec.record.x(:,3))
subplot(818)
plot(pp.t,pp.parts.recValve.record.x);
disp(['Number of CoolProp bugs were ' num2str(bugnumber)])

% Functions for Physical Plant
function Dx = process(t,x,pp)
% In this function, all the subprocesses are called, and their connections
% are arranged
    % States
    xGC = x(1:pp.parts.gc.nCell*3);
    xHPValve = x(pp.parts.gc.nCell*3+1);
    xRec = x(pp.parts.gc.nCell*3+2:pp.parts.gc.nCell*3+4);
    xRecValve = x(pp.parts.gc.nCell*3+5);
    pp.parts.gc.p = xGC(1:pp.parts.gc.nCell);
    pp.parts.gc.h = xGC(pp.parts.gc.nCell+1:2*pp.parts.gc.nCell);
    pp.parts.gc.d = xGC(2*pp.parts.gc.nCell+1:3*pp.parts.gc.nCell);
    pp.parts.HPValve.DmState = xHPValve;
    pp.parts.rec.p = xRec(1);
    pp.parts.rec.h = xRec(2);
    pp.parts.rec.d = xRec(3);
    pp.parts.recValve.DmState = xRecValve;
    % Static equations
    pp.parts.HPValve.enthalpy(pp.parts.gc.h(end));
    pp.parts.recValve.enthalpy(pp.parts.rec.hGas);
    pp.parts.rec.separation();
    % Inputs
    pp.Inputs.recValveInputs.dInlet = pp.parts.rec.dGas;
    pp.Inputs.recValveInputs.pInlet = pp.parts.rec.p;
    pp.Inputs.recValveInputs.hInlet = pp.parts.rec.hGas;
    pp.Inputs.recInputs.DmGas = pp.parts.recValve.DmState;
    pp.Inputs.recInputs.DmInlet = pp.parts.HPValve.DmState;
    pp.Inputs.recInputs.hInlet = pp.parts.HPValve.hOutlet;
    pp.Inputs.HPValveInputs.dInlet = pp.parts.gc.d(end);
    pp.Inputs.HPValveInputs.pInlet = pp.parts.gc.p(end);
    pp.Inputs.HPValveInputs.hInlet = pp.parts.gc.h(end);
    pp.Inputs.HPValveInputs.pOutlet = pp.parts.rec.p;
    pp.Inputs.gcInputs.DmOutlet = pp.parts.HPValve.DmState;
    % Receiver Valve
    DRecValve = pp.parts.recValve.process(t,xRecValve,pp.Inputs.recValveInputs);
    % Receiver
    DRec = pp.parts.rec.process(t,xRec,pp.Inputs.recInputs);
    % HP Valve
    DHPValve = pp.parts.HPValve.process(t,xHPValve,pp.Inputs.HPValveInputs);
    % Gas cooler
    DGC = pp.parts.gc.process(t,xGC,pp.Inputs.gcInputs);
    % Derivative
    Dx = [DGC; DHPValve; DRec; DRecValve];
end
function postProcess(X,pp)
% In this function, the states and the records are updated
    % States
    XGC = X(:,1:pp.parts.gc.nCell*3);
    XHPValve = X(:,pp.parts.gc.nCell*3+1);
    XRec = X(:,pp.parts.gc.nCell*3+2:pp.parts.gc.nCell*3+4);
    XRecValve = X(:,pp.parts.gc.nCell*3+5);
    % Recording
    pp.parts.gc.record.x = [pp.parts.gc.record.x; XGC];
    pp.parts.HPValve.record.x = [pp.parts.HPValve.record.x; XHPValve];
    pp.parts.rec.record.x = [pp.parts.rec.record.x; XRec];
    pp.parts.recValve.record.x = [pp.parts.recValve.record.x; XRecValve];
    % State setting for next iteration
    pp.parts.gc.p = XGC(end,1:pp.parts.gc.nCell)';
    pp.parts.gc.h = XGC(end,pp.parts.gc.nCell+1:2*pp.parts.gc.nCell)';
    pp.parts.gc.d = XGC(end,2*pp.parts.gc.nCell+1:3*pp.parts.gc.nCell)';
    pp.parts.rec.p = XRec(end,1);
    pp.parts.rec.h = XRec(end,2);
    pp.parts.rec.d = XRec(end,3);
    pp.parts.HPValve.DmState = XHPValve(end,1);
    pp.parts.recValve.DmState = XRecValve(end,1);
end