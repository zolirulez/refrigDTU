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
    % Handles for shorts
    gc = pp.parts.gc;
    HPValve = pp.parts.HPValve;
    rec = pp.parts.rec;
    recValve = pp.parts.recValve;
    gcInputs = pp.Inputs.gcInputs;
    HPValveInputs = pp.Inputs.HPValveInputs;
    recInputs = pp.Inputs.recInputs;
    recValveInputs = pp.Inputs.recValveInputs;
    % States
    xGC = x(1:gc.nCell*3);
    xHPValve = x(gc.nCell*3+1);
    xRec = x(gc.nCell*3+2:gc.nCell*3+4);
    xRecValve = x(gc.nCell*3+5);
    gc.p = xGC(1:gc.nCell);
    gc.h = xGC(gc.nCell+1:2*gc.nCell);
    gc.d = xGC(2*gc.nCell+1:3*gc.nCell);
    HPValve.DmState = xHPValve;
    rec.p = xRec(1);
    rec.h = xRec(2);
    rec.d = xRec(3);
    recValve.DmState = xRecValve;
    % Static equations
    pp.parts.HPValve.enthalpy(gc.h(end));
    pp.parts.recValve.enthalpy(rec.hGas);
    pp.parts.rec.separation();
    % Inputs
    recValveInputs.dInlet = rec.dGas;
    recValveInputs.pInlet = rec.p;
    recValveInputs.hInlet = rec.hGas;
    recInputs.DmGas = recValve.DmState;
    recInputs.DmInlet = HPValve.DmState;
    recInputs.hInlet = HPValve.hOutlet;
    HPValveInputs.dInlet = gc.d(end);
    HPValveInputs.pInlet = gc.p(end);
    HPValveInputs.hInlet = gc.h(end);
    HPValveInputs.pOutlet = rec.p;
    gcInputs.DmOutlet = HPValve.DmState;
    % Receiver Valve
    DRecValve = recValve.process(t,xRecValve,recValveInputs);
    % Receiver
    DRec = rec.process(t,xRec,recInputs);
    % HP Valve
    DHPValve = HPValve.process(t,xHPValve,HPValveInputs);
    % Gas cooler
    DGC = gc.process(t,xGC,gcInputs);
    % Derivative
    Dx = [DGC; DHPValve; DRec; DRecValve];
end
function postProcess(X,pp)
% In this function, the states and the records are updated
    % Handles for shorts
    gc = pp.parts.gc;
    HPValve = pp.parts.HPValve;
    rec = pp.parts.rec;
    recValve = pp.parts.recValve;
    % States
    XGC = X(:,1:gc.nCell*3);
    XHPValve = X(:,gc.nCell*3+1);
    XRec = X(:,gc.nCell*3+2:gc.nCell*3+4);
    XRecValve = X(:,gc.nCell*3+5);
    % Recording
    gc.record.x = [gc.record.x; XGC];
    HPValve.record.x = [HPValve.record.x; XHPValve];
    rec.record.x = [rec.record.x; XRec];
    recValve.record.x = [recValve.record.x; XRecValve];
    % State setting for next iteration
    gc.p = XGC(end,1:gc.nCell)';
    gc.h = XGC(end,gc.nCell+1:2*gc.nCell)';
    gc.d = XGC(end,2*gc.nCell+1:3*gc.nCell)';
    rec.p = XRec(end,1);
    rec.h = XRec(end,2);
    rec.d = XRec(end,3);
    HPValve.DmState = XHPValve(end,1);
    recValve.DmState = XRecValve(end,1);
end