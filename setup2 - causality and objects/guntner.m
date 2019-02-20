clearvars
global bugnumber 
bugnumber = 0;
% Initialization
% Gas Cooler
givenVolume = 0.0192;
Parameters.InnerTubeDiameter = 0.007;
Parameters.nParallelFlows = 5;
Parameters.OneTubelength = givenVolume/(Parameters.InnerTubeDiameter^2*pi/4);
Parameters.f1 = 0.0001; % In case of just a few cells, sensitivity is pretty low
Parameters.TauDm = 1;
% Gas Cooler, thermal
Parameters.NominalVolumeFlow = 3.33;
Parameters.ConductionRatio = 0.2;
Parameters.Qnominal = 74300;
Parameters.Tpi = 273+107;
Parameters.Tpo = 273+33;
Parameters.Tsi = 273+30;
Parameters.Tso = 273+40;
Parameters.TauTa = 5;
gc = HeatExchanger;
p = [88e5 85e5];
h = [450e3 305e3];
Ta = [32 50]+273;
iDm = 0.321;
Tau = 1;
ODEoptions = [];
nCell = 10;
gc.initialize(nCell,p,h,iDm,Ta,Parameters)
% HP Valve initialization
HPValve = Valve;
Kv = 0.8;
Dm = 0.321;
Initial.Dm = Dm;
Initial.h = 298e3;
HPValve.initialize(Kv,Tau,Initial);
HPValve.hInlet = gc.h(end);
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
K = 5e-7;
Ti = 5;
mn = 0;
mx = 1;
neg = -1;
refHP = 85e5;
refFan = 273.15 + 33;
timestep = 2;
initInt = 0.92;
PIHPValve.initialize(K,Ti,initInt,mn,mx,neg,timestep);
mx = 1;
K = 5e-2;
Ti = 5;
initInt = 0.52;
PIfan.initialize(K,Ti,initInt,mn,mx,neg,timestep);

% Connection object
connect = Connect;

% Initialization of physical plant object
ODEoptions.jancsi = 'jancsi';
pp = PhysicalPlant;
Parts.gc = gc;
Parts.HPValve = HPValve;
Parts.connect = Connect;
Parts.fan = fan;
Process = @process;
PreProcess = @preProcess; 
PostProcess = @postProcess; 
ODEoptions = [];
Initial = [gc.p; gc.h; gc.d; gc.Ta; gc.iDm; HPValve.Dm; fan.DV];
pp.initialize(Parts,Process,PostProcess,Initial,ODEoptions);
Inputs = [];
pp.initialInputs(Inputs);

% Simulation
itmax = 5000;
tic
for it = 1:itmax
    % Controller
    pp.Inputs.HPValveCR =  PIHPValve.react(refHP,pp.parts.gc.p(end));
    pp.Inputs.fanCR = PIfan.react(refFan+1,pp.parts.gc.T(end));
    % Physical Plant
    pp.timestep(((it-1):it)*timestep);
    % Plotting
    if ~rem(it,1)
        pause(0.1)
        subplot(321)
        plot(pp.t,pp.parts.gc.record.x(:,1:nCell))
        subplot(322)
        plot(pp.t,pp.parts.gc.record.x(:,1*nCell+1:2*nCell))
        subplot(323)
        plot(pp.t,pp.parts.gc.record.x(:,2*nCell+1:3*nCell))
        subplot(324)
        plot(pp.t,pp.parts.HPValve.record.x);
        subplot(325)
        plot(pp.t,pp.parts.fan.record.x);
        subplot(326)
        plot(pp.t,(gc.hInlet-gc.record.x(:,2*nCell)).*(HPValve.record.x+gc.DmInlet)/2);
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
    connect = pp.parts.connect;
    fan = pp.parts.fan;
    % Disturbances
    HPValveCR = pp.Inputs.HPValveCR;
    fanCR = pp.Inputs.fanCR;
    % States
    xGC = x(1:gc.nCell*5-1);
    xHPValve = x(gc.nCell*5-1+1);
    xFan = x(gc.nCell*5-1+2);
    gc.p = xGC(1:gc.nCell);
    gc.h = xGC(gc.nCell+1:2*gc.nCell);
    gc.d = xGC(2*gc.nCell+1:3*gc.nCell);
    gc.Ta = xGC(3*gc.nCell+1:4*gc.nCell);
    gc.iDm = xGC(4*gc.nCell+1:5*gc.nCell-1);
    HPValve.Dm = xHPValve;
    fan.DV = xFan;
    % Boundary Conditions
    gc.DmInlet = 0.321;
    gc.hInlet = 530e3;
    HPValve.pOutlet = 38e5;
    % Connections: paired components
    connect.inlet(HPValve,gc,{'p','h','d'});
    connect.outlet(gc,HPValve,{'Dm'});
    connect.inlet(gc,fan,{'T','DV'});
    % HP Valve
    DHPValve = HPValve.process(t,xHPValve,HPValveCR);
    % Gas cooler
    DGC = gc.process(t,xGC);
    % Fan
    DFan = fan.process(t,xFan,fanCR);
    % Derivative
    Dx = [DGC; DHPValve; DFan];
end
function postProcess(X,pp)
% In this function, the states and the records are updated
    % Handles for shorts
    gc = pp.parts.gc;
    HPValve = pp.parts.HPValve;
    connect = pp.parts.connect;
    fan = pp.parts.fan;
    % States
    XGC = X(:,1:gc.nCell*5-1);
    XHPValve = X(:,gc.nCell*5-1+1);
    XFan = X(:,gc.nCell*5-1+2);
    % Recording
    gc.record.x = [gc.record.x; XGC];
    HPValve.record.x = [HPValve.record.x; XHPValve];
    fan.record.x = [fan.record.x; XFan];
    % State setting for next iteration
    gc.p = XGC(end,1:gc.nCell)';
    gc.h = XGC(end,gc.nCell+1:2*gc.nCell)';
    gc.d = XGC(end,2*gc.nCell+1:3*gc.nCell)';
    gc.d = XGC(end,3*gc.nCell+1:4*gc.nCell)';
    gc.iDm = XGC(end,4*gc.nCell+1:5*gc.nCell-1)';
    HPValve.Dm = XHPValve(end,1);
    fan.DV = XFan(end,1);
    % Connections: paired components
    connect.inlet(HPValve,gc,{'p','h','d'});
    connect.outlet(gc,HPValve,{'Dm'});
    connect.inlet(gc,fan,{'T','DV'});
end