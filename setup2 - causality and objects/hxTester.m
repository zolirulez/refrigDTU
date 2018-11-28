clearvars
givenVolume = 0.0192;
Parameters.InnerTubeDiameter = 0.0192;
Parameters.nParallelFlows = 20;
Parameters.OneTubelength = givenVolume/(Parameters.InnerTubeDiameter^2*pi/4);
Parameters.f1 = 0.01;
Parameters
hx = HeatExchanger;
p = [84.8e5 84.8e5-1.57e5];
h = [525e3 298e3];
tau = [0.001 0.01 0.01];
ODEoptions = [];
nCell = 10;
hx.initialize(nCell,p,h,Parameters,tau,ODEoptions)
Dm = 0.3186;
inputs.hInlet = 525e3;
DQ = -ones(nCell,1)*Dm*(inputs.hInlet-h(2))/nCell/Parameters.nParallelFlows;
inputs.DmInlet = Dm;
inputs.DmOutlet = Dm;
tic
timestep = 1;
for it = 1:500
    it
    inputs.DQ = DQ + 0*ones(nCell,1)*sin(2*pi*it*timestep)*1000;
    hx.timestep([(it-1)*timestep it*timestep],inputs);
end
toc
subplot(311)
plot(hx.record.t,hx.record.x(:,1:nCell))
subplot(312)
plot(hx.record.t,hx.record.x(:,1*nCell+1:2*nCell))
subplot(313)
plot(hx.record.t,hx.record.x(:,2*nCell+1:3*nCell))