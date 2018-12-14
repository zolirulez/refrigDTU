clearvars
global bugnumber partials ratios derivatives
partials = [];
bugnumber = 0;
derivatives = [];
givenVolume = 0.0192;
Parameters.InnerTubeDiameter = 0.0192;
Parameters.nParallelFlows = 20;
Parameters.OneTubelength = givenVolume/(Parameters.InnerTubeDiameter^2*pi/4);
Parameters.f1 = 0.01; % In case of just a few cells, sensitivity is pretty low
Parameters
hx = HeatExchanger;
p = [85e5 85e5-1.57e5];
h = [525e3 298e3];
tau = [0.1 0.01 0.01];
ODEoptions = [];
nCell = 4;
hx.initialize(nCell,p,h,Parameters,tau,ODEoptions)
Dm = 0.3186;
inputs.hInlet = 525e3;
DQ = -ones(nCell,1)*Dm*(inputs.hInlet-h(2))/nCell/Parameters.nParallelFlows;
inputs.DmInlet = Dm;
inputs.DmOutlet = Dm;
tic
timestep = 1;
itmax = 500;
t_Dm = timestep:timestep:itmax*timestep;
Dm = zeros(hx.nCell+1,length(t_Dm));
for it = 1:itmax
    inputs.DQ = DQ + 0*ones(nCell,1)*sin(2*pi*it/itmax*10)*100;
    hx.timestep([(it-1)*timestep it*timestep],inputs);
    Dm(:,it) = hx.Dm;
end
toc
figure(1)
subplot(411)
plot(hx.record.t,hx.record.x(:,1:nCell))
subplot(412)
plot(hx.record.t,hx.record.x(:,1*nCell+1:2*nCell))
subplot(413)
plot(hx.record.t,hx.record.x(:,2*nCell+1:3*nCell))
subplot(414)
plot(t_Dm,Dm);
disp(['Number of CoolProp bugs were ' num2str(bugnumber)])
figure(2)
plot(partials')
figure(4)
plot(derivatives(1:4,:)')
figure(5)
plot(derivatives(5:8,:)')