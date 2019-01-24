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
Initial.Dm = Dm;
Initial.h = 300e3;
HPValve.initialize(Kv,Tau,Initial);
HPValve.hInlet = gc.h(end);
% Gas Cooler input initialization
Inputs.gcDQ = -ones(nCell,1)*Dm*(525e3-300e3)/nCell/Parameters.nParallelFlows;
% Receiver initialization
recVolume = 0.133;
rec = Receiver;
pRec = 35e5;
hRec = 250e3;
DmGas = 0.126;
rec.initialize(pRec,hRec,recVolume,ODEoptions);
rec.separation(); % For the receiver valve
rec.liquid.DmOutlet = Dm - DmGas;
% Receiver Valve initialization
recValve = Valve;
Kv = 2;
Initial.Dm = DmGas/5;
Initial.h = 440e3;
recValve.initialize(Kv,Tau,Initial);
recValve.h = rec.gas.h;
% Joint initialization
MTJoint = Joint;
RecJoint = Joint;
% Cooler initialization
cooler = Evaporator;
coolerVolume = givenVolume*5;
pCooler = 30e5;
hOutletVirtual = 440e3;
hCooler = hOutletVirtual;
cooler.initialize(pCooler,hCooler,hOutletVirtual,coolerVolume,ODEoptions);
% Inputs.coolerInputs.DmOutlet = 0.151;
% Inputs.coolerInputs.DQ = 36e3;
% Inputs.coolerInputs.hInlet = 195e3;
Inputs.coolerDQ = 36e3;
cooler.DmInlet = 0.151;
% Compressor initialization
comp = Compressor;
IsentropicEfficiency = 0.65;
compDensity = CoolProp.PropsSI('D','H',450e3,'P',30e5,'CO2');
Displacement = Dm/compDensity/20; % TODO
Initial.Dm = Dm;
Initial.h = 525e3;
Tau = 1;
comp.initialize(IsentropicEfficiency,Displacement,Tau,Initial);
comp.hInlet = 450e3;
comp.pInlet = 30e5;
comp.pOutlet = 85e5;
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
timestep = 0.1;
initInt = 0;
PIHPValve.initialize(K,Ti,initInt,mn,mx,neg,timestep);
K = 10^-6.5;
Ti = 50;
PIrecValve.initialize(K,Ti,initInt,mn,mx,neg,timestep);
mx = 5;
K = 1e-7;
Ti = 100;
PIcomp.initialize(K,Ti,initInt,mn,mx,neg,timestep);

% Connection object
connect = Connect;

% Initialization of physical plant object
ODEoptions.jancsi = 'jancsi';
pp = PhysicalPlant;
Parts.gc = gc;
Parts.HPValve = HPValve;
Parts.rec = rec;
Parts.cooler = cooler;
Parts.comp = comp;
Parts.MTJoint = MTJoint;
Parts.RecJoint = RecJoint;
Parts.recValve = recValve;
Parts.connect = Connect;
Process = @process;
PreProcess = @preProcess; 
PostProcess = @postProcess; 
ODEoptions = [];
Initial = [gc.p; gc.h; gc.d; HPValve.Dm; rec.p; rec.h; rec.d;...
    recValve.Dm; cooler.p; cooler.h; cooler.d; comp.Dm];
pp.initialize(Parts,Process,PostProcess,Initial,ODEoptions);
pp.initialInputs(Inputs);

% Simulation
itmax = 5;
tic
for it = 1:itmax
    % Controller
    pp.Inputs.recValveCR = PIrecValve.react(refRec,pp.parts.rec.p);
    pp.Inputs.HPValveCR =  PIHPValve.react(refHP,pp.parts.gc.p(end));
    pp.Inputs.compFreq = PIcomp.react(refMT,pp.parts.cooler.p);
    % Physical Plant
    pp.timestep(((it-1):it)*timestep);
end
toc
figure(1)
subplot(621)
plot(pp.t,pp.parts.gc.record.x(:,1:nCell))
subplot(622)
plot(pp.t,pp.parts.gc.record.x(:,1*nCell+1:2*nCell))
subplot(623)
plot(pp.t,pp.parts.gc.record.x(:,2*nCell+1:3*nCell))
subplot(624)
plot(pp.t,pp.parts.HPValve.record.x);
subplot(625)
plot(pp.t,pp.parts.rec.record.x(:,1))
subplot(626)
plot(pp.t,pp.parts.rec.record.x(:,2))
subplot(627)
plot(pp.t,pp.parts.rec.record.x(:,3))
subplot(628)
plot(pp.t,pp.parts.recValve.record.x)
subplot(629)
plot(pp.t,pp.parts.cooler.record.x(:,1))
subplot(6,2,10)
plot(pp.t,pp.parts.cooler.record.x(:,2))
subplot(6,2,11)
plot(pp.t,pp.parts.cooler.record.x(:,3))
subplot(6,2,12)
plot(pp.t,pp.parts.comp.record.x)
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
    MTJoint = pp.parts.MTJoint;
    RecJoint = pp.parts.RecJoint;
    connect = pp.parts.connect;
    % Disturbances
    gcDQ = pp.Inputs.gcDQ;
    coolerDQ = pp.Inputs.coolerDQ;
    recValveCR = pp.Inputs.recValveCR;
    HPValveCR = pp.Inputs.HPValveCR;
    compFreq = pp.Inputs.compFreq;
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
    HPValve.Dm = xHPValve;
    rec.p = xRec(1);
    rec.h = xRec(2);
    rec.d = xRec(3);
    recValve.Dm = xRecValve;
    cooler.p = xCooler(1);
    cooler.h = xCooler(2);
    cooler.d = xCooler(3);
    comp.Dm = xComp;
    % Static equations
    MTJoint.noAccummulation(cooler.p,...
        [recValve.Dm; 0.046],[recValve.h; 530e3],-comp.Dm,cooler.h);
    RecJoint.noAccummulation(rec.p,...
        -0.046,rec.liquid.h,-cooler.DmInlet,rec.liquid.h);
    % Connections 
    connect.inlet(comp,MTJoint,{'p','h','d'});
    connect.outlet(comp,gc,{'p'});
    connect.outlet(recValve,MTJoint,{'p'});
    connect.inlet(recValve,rec.gas,{'p','h','d'});
    connect.outlet(rec,recValve,{'Dm'},'gas');
    connect.inlet(rec,HPValve,{'Dm','h'});
    connect.outlet(rec,RecJoint,{'Dm'},'liquid');
    connect.outlet(HPValve,rec,{'p'});
    connect.inlet(HPValve,gc,{'p','h','d'});
    connect.inlet(gc,comp,{'Dm','h'});
    connect.outlet(gc,HPValve,{'Dm'});
    connect.outlet(cooler,MTJoint,{'Dm'});
    connect.inlet(cooler,RecJoint,{'h'});
    % Cooler
    DCooler = cooler.process(t,xCooler,coolerDQ);
    % Compressor
    DComp = comp.process(t,xComp,compFreq);
    % Receiver Valve
    DRecValve = recValve.process(t,xRecValve,recValveCR);
    % Receiver
    DRec = rec.process(t,xRec);
    % HP Valve
    DHPValve = HPValve.process(t,xHPValve,HPValveCR);
    % Gas cooler
    DGC = gc.process(t,xGC,gcDQ);
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
    HPValve.Dm = XHPValve(end,1);
    recValve.Dm = XRecValve(end,1);
    comp.Dm = XComp(end,1);
end