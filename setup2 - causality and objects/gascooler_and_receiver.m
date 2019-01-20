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
% Joint initialization
joint = Joint;
% Cooler initialization
cooler = Evaporator;
coolerVolume = givenVolume*5;
pCooler = 30e5;
hCooler = 400e3;
hOutlet = 433e3;
cooler.initialize(pCooler,hCooler,hOutlet,coolerVolume,ODEoptions);
Inputs.coolerInputs.DmOutlet = 0.151;
Inputs.coolerInputs.DQ = 36e3;
Inputs.coolerInputs.hInlet = 195e3;
% Compressor initialization
comp = Compressor;
IsentropicEfficiency = 0.9;
compDensity = CoolProp.PropsSI('D','H',450e3,'P',30e5,'CO2');
Displacement = Dm/compDensity/20; % TODO
Initial = Dm;
Tau = 1;
comp.initialize(IsentropicEfficiency,Displacement,Tau,Initial);
Inputs.compInputs.dInlet = compDensity;
Inputs.compInputs.pInlet = 30e5;
Inputs.compInputs.pOutlet = 85e5;
Inputs.compInputs.hInlet = 450e3;
comp.enthalpy(Inputs.compInputs.pInlet,Inputs.compInputs.pOutlet,...
    Inputs.compInputs.hInlet);
% Controller initialization
PIHPValve = PIController;
PIrecValve = PIController;
PIcomp = PIController;
K = 1e-6;
Ti = 50;
mn = 0;
mx = 1;
neg = -1;
refHP = 85e5;
refRec = 38e5;
refMT = 30e5;
timestep = 5;
initInt = 0;
PIHPValve.initialize(K,Ti,initInt,mn,mx,neg,timestep);
K = 10^-6.5;
Ti = 50;
PIrecValve.initialize(K,Ti,initInt,mn,mx,neg,timestep);
mx = 5;
Ti = 100;
PIcomp.initialize(K,Ti,initInt,mn,mx,neg,timestep);



% Initialization of physical plant object
ODEoptions.jancsi = 'jancsi';
pp = PhysicalPlant;
Parts.gc = gc;
Parts.HPValve = HPValve;
Parts.rec = rec;
Parts.cooler = cooler;
Parts.comp = comp;
Parts.joint = joint;
Parts.recValve = recValve;
Process = @process;
PreProcess = @preProcess; 
PostProcess = @postProcess; 
ODEoptions = [];
Initial = [gc.p; gc.h; gc.d; HPValve.DmState; rec.p; rec.h; rec.d;...
    recValve.DmState; cooler.p; cooler.h; cooler.d; comp.DmState];
pp.initialize(Parts,Process,PostProcess,Initial,ODEoptions);
pp.initialInputs(Inputs);

% Simulation
itmax = 100;
tic
for it = 1:itmax
    % Controller
    pp.Inputs.recValveInputs.capacityRatio = PIrecValve.react(refRec,pp.parts.rec.p);
    pp.Inputs.HPValveInputs.capacityRatio =  PIHPValve.react(refHP,pp.parts.gc.p(end));
    pp.Inputs.compInputs.frequency =  PIcomp.react(refMT,pp.parts.cooler.p);
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
    comp = pp.parts.comp;
    cooler = pp.parts.cooler;
    joint = pp.parts.joint;
    gcInputs = pp.Inputs.gcInputs;
    HPValveInputs = pp.Inputs.HPValveInputs;
    recInputs = pp.Inputs.recInputs;
    recValveInputs = pp.Inputs.recValveInputs;
    compInputs = pp.Inputs.compInputs;
    coolerInputs = pp.Inputs.coolerInputs;
    % States
    xGC = x(1:gc.nCell*3);
    xHPValve = x(gc.nCell*3+1);
    xRec = x(gc.nCell*3+2:gc.nCell*3+4);
    xRecValve = x(gc.nCell*3+5);
    xCooler = x(gc.nCell*3+6:gc.nCell*3+8);
    xComp = x(gc.nCell*3+9);
    gc.p = xGC(1:gc.nCell);
    gc.h = xGC(gc.nCell+1:2*gc.nCell);
    gc.d = xGC(2*gc.nCell+1:3*gc.nCell);
    HPValve.DmState = xHPValve;
    rec.p = xRec(1);
    rec.h = xRec(2);
    rec.d = xRec(3);
    recValve.DmState = xRecValve;
    cooler.p = xCooler(1);
    cooler.h = xCooler(2);
    cooler.d = xCooler(3);
    comp.DmState = xComp;
    % Static equations
    HPValve.enthalpy(gc.h(end));
    recValve.enthalpy(rec.hGas);
    rec.separation();
    [coolerInputs.DmOutlet, compInputs.hInlet] =...
        joint.noAccummulation([recValve.DmState; 0.0146],[recValve.hOutlet; 500e3],...
        -comp.DmState,cooler.hOutlet); % TODO teh values
    % Inputs
    compInputs.dInlet = CoolProp.PropsSI('D','H',compInputs.hInlet,'P',cooler.p,'CO2');
    compInputs.pInlet = cooler.p;
    recValveInputs.pOutlet = cooler.p;
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
    gcInputs.hInlet = comp.hOutlet;
    gcInputs.DmInlet = comp.DmState;
    % Cooler
    DCooler = cooler.process(t,xCooler,coolerInputs);
    % Compressor
    DComp = comp.process(t,xComp,compInputs);
    % Receiver Valve
    DRecValve = recValve.process(t,xRecValve,recValveInputs);
    % Receiver
    DRec = rec.process(t,xRec,recInputs);
    % HP Valve
    DHPValve = HPValve.process(t,xHPValve,HPValveInputs);
    % Gas cooler
    DGC = gc.process(t,xGC,gcInputs);
    % Derivative
    Dx = [DGC; DHPValve; DRec; DRecValve; DCooler; DComp];
end
function postProcess(X,pp)
% In this function, the states and the records are updated
    % Handles for shorts
    gc = pp.parts.gc;
    HPValve = pp.parts.HPValve;
    rec = pp.parts.rec;
    recValve = pp.parts.recValve;
    comp = pp.parts.comp;
    cooler = pp.parts.cooler;
    % States
    XGC = X(:,1:gc.nCell*3);
    XHPValve = X(:,gc.nCell*3+1);
    XRec = X(:,gc.nCell*3+2:gc.nCell*3+4);
    XRecValve = X(:,gc.nCell*3+5);
    XCooler = X(:,gc.nCell*3+6:gc.nCell*3+8);
    XComp = X(:,gc.nCell*3+9);
    % Recording
    gc.record.x = [gc.record.x; XGC];
    HPValve.record.x = [HPValve.record.x; XHPValve];
    rec.record.x = [rec.record.x; XRec];
    recValve.record.x = [recValve.record.x; XRecValve];
    cooler.record.x = [cooler.record.x; XCooler];
    comp.record.x = [comp.record.x; XComp];
    % State setting for next iteration
    gc.p = XGC(end,1:gc.nCell)';
    gc.h = XGC(end,gc.nCell+1:2*gc.nCell)';
    gc.d = XGC(end,2*gc.nCell+1:3*gc.nCell)';
    rec.p = XRec(end,1);
    rec.h = XRec(end,2);
    rec.d = XRec(end,3);
    cooler.p = XCooler(end,1);
    cooler.h = XCooler(end,2);
    cooler.d = XCooler(end,3);
    HPValve.DmState = XHPValve(end,1);
    recValve.DmState = XRecValve(end,1);
    comp.DmState = XComp(end,1);
end