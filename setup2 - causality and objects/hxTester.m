clearvars
givenVolume = 0.0192;
Parameters.InnerTubeDiameter = 0.0192;
Parameters.nParallelFlows = 20;
Parameters.OneTubelength = givenVolume/(Parameters.InnerTubeDiameter^2*pi/4);
Parameters.f1 = 0.008;
Parameters
hx = HeatExchanger;
p = [90e5 89e5];
h = [550e3 300e3];
tau = [0.1 0.01 0.01];
ODEoptions = [];
nCell = 10;
hx.initialize(10,p,h,Parameters,tau,ODEoptions)
Dm = 0.45;
inputs.hInlet = 545e3;
inputs.DQ = -ones(nCell,1)*Dm*(inputs.hInlet-h(2))/nCell/Parameters.nParallelFlows;
inputs.DmInlet = Dm;
inputs.DmOutlet = Dm;
tic
timestep = 0.001;
for it = 1:100
    hx.timestep([(it-1)*timestep it*timestep],inputs);
end
toc
subplot(311)
plot(hx.record.t,hx.record.x(:,1:10))
subplot(312)
plot(hx.record.t,hx.record.x(:,11:20))
subplot(313)
plot(hx.record.t,hx.record.x(:,21:30))