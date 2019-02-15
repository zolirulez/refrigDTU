clearvars
global bugnumber 
bugnumber = 0;
% Initialization
% Gas Cooler
givenVolume = 0.0192;
Parameters.InnerTubeDiameter = 0.007;
Parameters.nParallelFlows = 5;
Parameters.OneTubelength = givenVolume/(Parameters.InnerTubeDiameter^2*pi/4);
Parameters.f1 = 0.005; % In case of just a few cells, sensitivity is pretty low
% Gas Cooler, thermal
Parameters.NominalVolumeFlow = 3.33;
Parameters.ConductionRatio = 0.2;
Parameters.Qnominal = 74300;
Parameters.Tpi = 273+107;
Parameters.Tpo = 273+33;
Parameters.Tsi = 273+30;
Parameters.Tso = 273+40;
gc = HeatExchanger;
p = [90e5 87e5];
h = [450e3 320e3];
Tau = 1;
ODEoptions = [];
nCell = 4;
gc.initialize(nCell,p,h,Parameters)
% HP Valve initialization
HPValve = Valve;
Kv = 0.8;
Dm = 0.321;
Initial.Dm = Dm;
Initial.h = 298e3;
HPValve.initialize(Kv,Tau,Initial);
HPValve.hInlet = gc.h(end);
% Receiver initialization
recVolume = 0.133;
rec = Receiver;
pRec = 38e5;
hRec = 280e3;
DmGas = 0.126;
rec.initialize(pRec,hRec,recVolume);
rec.separation(); % For the receiver valve
rec.liquid.DmOutlet = Dm - DmGas;
% Receiver Valve initialization
recValve = Valve;
Kv = 2;
Initial.Dm = DmGas;
Initial.h = 440e3;
recValve.initialize(Kv,Tau,Initial);
recValve.hInlet = rec.gas.h;
% Joint initialization
MTJoint = Joint;
RecJoint = Joint;
% Cooler initialization
cooler = Evaporator;
coolerVolume = givenVolume*5;
pCooler = 26.5e5;
hOutletVirtual = 440e3;
hCooler = hOutletVirtual;
cooler.initialize(pCooler,hCooler,hOutletVirtual,coolerVolume);
Inputs.coolerDQ = 36e3;
cooler.DmInlet = 0.151;
% Compressor initialization
comp = Compressor;
IsentropicEfficiency = 0.65;
compDensity = CoolProp.PropsSI('D','H',450e3,'P',26.5e5,'CO2');
Displacement = (9.6+6.5+4.3)/1450/60;
Initial.Dm = Dm;
Initial.h = 525e3;
Tau = 1;
comp.initialize(IsentropicEfficiency,Displacement,Tau,Initial);
comp.hInlet = 400e3;
comp.pInlet = 26.5e5;
comp.pOutlet = 85e5;
% Fan initialization
fan = Fan;
Initial.DV = 3.33;
Initial.T = 273.15+30;
maxVolumeFlow = 6.66;
Tau = 1;
fan.initialize(maxVolumeFlow,Tau,Initial);
% Controller initialization
PIHPValve = PIController;
PIrecValve = PIController;
PIcomp = PIController;
PIfan = PIController;
K = 4e-6;
Ti = 5;
mn = 0;
mx = 1;
neg = -1;
refHP = 85e5;
refRec = 38e5;
refMT = 26.5e5;
refFan = 273.15 + 33;
timestep = 0.5;
initInt = 0.25;
PIHPValve.initialize(K,Ti,initInt,mn,mx,neg,timestep);
K = 2e-6;
Ti = 5;
initInt = 0.25;
PIrecValve.initialize(K,Ti,initInt,mn,mx,neg,timestep);
mx = 5;
K = 2e-6;
Ti = 10;
initInt = 0.25;
PIcomp.initialize(K,Ti,initInt,mn,mx,neg,timestep);
K = 1e-1;
Ti = 5;
initInt = 1;
PIfan.initialize(K,Ti,initInt,mn,mx,neg,timestep);

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
Parts.fan = fan;
Process = @process;
PreProcess = @preProcess; 
PostProcess = @postProcess; 
ODEoptions = [];
Initial = [gc.p; gc.h; gc.d; HPValve.Dm; rec.p; rec.h; rec.d;...
    recValve.Dm; cooler.p; cooler.h; cooler.d; comp.Dm; fan.DV];
pp.initialize(Parts,Process,PostProcess,Initial,ODEoptions);
pp.initialInputs(Inputs);

% Simulation
itmax = 5000;
tic
for it = 1:itmax
    % Controller
    pp.Inputs.recValveCR = PIrecValve.react(refRec,pp.parts.rec.p);
    pp.Inputs.HPValveCR =  PIHPValve.react(refHP,pp.parts.gc.p(end));
    pp.Inputs.compFreq = PIcomp.react(refMT,pp.parts.cooler.p);
    pp.Inputs.fanCR = PIfan.react(refFan+2,pp.parts.gc.T(end));
    % Physical Plant
    pp.timestep(((it-1):it)*timestep);
    % Plotting
    if ~rem(it,10)
        pause(0.1)
        subplot(521)
        plot(pp.t,pp.parts.gc.record.x(:,1:nCell))
        subplot(522)
        plot(pp.t,pp.parts.gc.record.x(:,1*nCell+1:2*nCell))
        subplot(523)
        plot(pp.t,pp.parts.HPValve.record.x);
        subplot(524)
        plot(pp.t,pp.parts.rec.record.x(:,1))
        subplot(525)
        plot(pp.t,pp.parts.rec.record.x(:,2))
        subplot(526)
        plot(pp.t,pp.parts.recValve.record.x)
        subplot(527)
        plot(pp.t,pp.parts.cooler.record.x(:,1))
        subplot(528)
        plot(pp.t,pp.parts.cooler.record.x(:,2))
        subplot(529)
        plot(pp.t,pp.parts.comp.record.x)
        subplot(5,2,10)
        plot(pp.t,pp.parts.fan.record.x)
    end
end
toc
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
    fan = pp.parts.fan;
    % Disturbances
    coolerDQ = pp.Inputs.coolerDQ;
    recValveCR = pp.Inputs.recValveCR;
    HPValveCR = pp.Inputs.HPValveCR;
    compFreq = pp.Inputs.compFreq;
    fanCR = pp.Inputs.fanCR;
    % States
    xGC = x(1:gc.nCell*3);
    xHPValve = x(gc.nCell*3+1);
    xRec = x(gc.nCell*3+2:gc.nCell*3+4);
    xRecValve = x(gc.nCell*3+5);
    xCooler = x(gc.nCell*3+6:gc.nCell*3+8);
    xComp = x(gc.nCell*3+9);
    xFan = x(gc.nCell*3+10);
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
    fan.DV = xFan;
    % Static equations: inner propagation of states
    rec.separation;
    % Connections: joints
    LTDm = 11000/(450e3-rec.liquid.h); % Boundary condition
    recValve.enthalpy;
    MTJoint.noAccummulation(cooler.p,...
        [recValve.Dm; LTDm],[recValve.h; 530e3],-comp.Dm,cooler.h);
    comp.enthalpy;
    RecJoint.noAccummulation(rec.p,...
        -LTDm,rec.liquid.h,-cooler.DmInlet,rec.liquid.h);
    % Connections: paired components
    connect.inlet(comp,MTJoint,{'p','h','d'});
    connect.outlet(comp,gc,{'p'});
    connect.outlet(recValve,MTJoint,{'p'});
    connect.inlet(recValve,rec.gas,{'p','h','d'});
    connect.outlet(rec,recValve,{'Dm'},'gas');
    connect.outlet(rec,RecJoint,{'Dm'},'liquid');
    connect.outlet(HPValve,rec,{'p'});
    connect.inlet(HPValve,gc,{'p','h','d'});
    connect.inlet(gc,comp,{'Dm','h'});
    connect.outlet(gc,HPValve,{'Dm'});
    connect.inlet(gc,fan,{'T','DV'});
    connect.outlet(cooler,MTJoint,{'Dm'});
    connect.inlet(cooler,RecJoint,{'h'});
    % Static equations: inner propagations of connected inlets
    HPValve.enthalpy;
    % Connections: paired components
    connect.inlet(rec,HPValve,{'Dm','h'});    
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
    DGC = gc.process(t,xGC);
    % Fan
    DFan = fan.process(t,xFan,fanCR);
    % Derivative
    Dx = [DGC; DHPValve; DRec; DRecValve; DCooler; DComp; DFan];
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
    MTJoint = pp.parts.MTJoint;
    RecJoint = pp.parts.RecJoint;
    connect = pp.parts.connect;
    fan = pp.parts.fan;
    % States
    XGC = X(:,1:gc.nCell*3);
    XHPValve = X(:,gc.nCell*3+1);
    XRec = X(:,gc.nCell*3+2:gc.nCell*3+4);
    XRecValve = X(:,gc.nCell*3+5);
    XCooler = X(:,gc.nCell*3+6:gc.nCell*3+8);
    XComp = X(:,gc.nCell*3+9);
    XFan = X(:,gc.nCell*3+10);
    % Recording
    gc.record.x = [gc.record.x; XGC];
    HPValve.record.x = [HPValve.record.x; XHPValve];
    rec.record.x = [rec.record.x; XRec];
    recValve.record.x = [recValve.record.x; XRecValve];
    cooler.record.x = [cooler.record.x; XCooler];
    comp.record.x = [comp.record.x; XComp];
    fan.record.x = [fan.record.x; XFan];
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
    fan.DV = XFan(end,1);
    % Static equations: inner propagation of states
    rec.separation;
    % Connections: joints
    LTDm = 11000/(450e3-rec.liquid.h); % Boundary condition
    recValve.enthalpy;
    MTJoint.noAccummulation(cooler.p,...
        [recValve.Dm; LTDm],[recValve.h; 530e3],-comp.Dm,cooler.h);
    comp.enthalpy;
    RecJoint.noAccummulation(rec.p,...
        -LTDm,rec.liquid.h,-cooler.DmInlet,rec.liquid.h);
    % Connections: paired components
    connect.inlet(comp,MTJoint,{'p','h','d'});
    connect.outlet(comp,gc,{'p'});
    connect.outlet(recValve,MTJoint,{'p'});
    connect.inlet(recValve,rec.gas,{'p','h','d'});
    connect.outlet(rec,recValve,{'Dm'},'gas');
    connect.outlet(rec,RecJoint,{'Dm'},'liquid');
    connect.outlet(HPValve,rec,{'p'});
    connect.inlet(HPValve,gc,{'p','h','d'});
    connect.inlet(gc,comp,{'Dm','h'});
    connect.outlet(gc,HPValve,{'Dm'});
    connect.inlet(gc,fan,{'T','DV'});
    connect.outlet(cooler,MTJoint,{'Dm'});
    connect.inlet(cooler,RecJoint,{'h'});
    % Static equations: inner propagations of connected inlets
    HPValve.enthalpy;
    % Connections: paired components
    connect.inlet(rec,HPValve,{'Dm','h'});   
end