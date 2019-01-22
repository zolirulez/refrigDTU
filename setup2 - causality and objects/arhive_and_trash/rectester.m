clearvars
global bugnumber partials 
partials = [];
bugnumber = 0;
Volume = 0.133;
rec = Receiver;
p = 38e5;
h = 300e3;
ODEoptions = [];
rec.initialize(p,h,Volume,ODEoptions);
inputs.DmInlet = 0.321;
inputs.DmGas = 0.126;
inputs.DmLiquid = inputs.DmInlet - inputs.DmGas;
inputs.hInlet = 294.75e3;
tic
timestep = 1;
itmax = 500;
for it = 1:itmax
    rec.timestep([(it-1)*timestep it*timestep],inputs);
end
toc
figure(1)
subplot(311)
plot(rec.record.t,rec.record.x(:,1))
subplot(312)
plot(rec.record.t,rec.record.x(:,2))
subplot(313)
plot(rec.record.t,rec.record.x(:,3))
disp(['Number of CoolProp bugs were ' num2str(bugnumber)])
figure(2)
plot(partials')