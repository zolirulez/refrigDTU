clearvars
close all
global bugnumber partials 
partials = [];
bugnumber = 0;
Volume = 0.0192*5;
ev = Evaporator;
p = 30e5;
h = 400e3;
hOutlet = 433e3;
ODEoptions = [];
ev.initialize(p,h,hOutlet,Volume,ODEoptions);
inputs.DmOutlet = 0.151;
inputs.DQ = 36e3;
tic
timestep = 1;
itmax = 500;
for it = 1:itmax
    if it < itmax/2
        inputs.hInlet = 186e3;
    else
        inputs.hInlet = -inputs.DQ/inputs.DmOutlet + hOutlet;
    end
    ev.timestep([(it-1)*timestep it*timestep],inputs);
end
toc
figure(1)
subplot(311)
plot(ev.record.t,ev.record.x(:,1))
subplot(312)
plot(ev.record.t,ev.record.x(:,2))
subplot(313)
plot(ev.record.t,ev.record.x(:,3))
disp(['Number of CoolProp bugs were ' num2str(bugnumber)])
figure(2)
plot(partials')