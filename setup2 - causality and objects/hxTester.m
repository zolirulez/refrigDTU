hx = HeatExchanger;
hx.initialize(10,[90e5 89e5],[550e3 300e3],0.0192,[0.1 0.1 0.01],[])
inputs.DQ = -[1 1 1 1 1 1 1 1 1 1]'*0.05*245e3/10;
inputs.DmInlet = 0.05;
inputs.DmOutlet = 0.05;
inputs.hInlet = 545e3;
tic
hx.timestep([0 10],inputs);
toc
subplot(311)
plot(hx.record.t,hx.record.x(:,31:40))
subplot(312)
plot(hx.record.t,hx.record.x(:,41:50))
subplot(313)
plot(hx.record.t,hx.record.x(:,51:60))