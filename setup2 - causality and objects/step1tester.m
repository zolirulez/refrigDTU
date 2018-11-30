% Sketch for the whole simulation
timestep = 0.1;
tic
for it = 1:1000
    % Pi controllers provide frequency and capacityRatio
    freq = PIcomp(ev.p(end));
    cR = PIvalve(gx.p(end));
    % Compressor
    [Dm, h] = comp.process(freq,ev.d(end),ev.p(end),gc.p(1),ev.h(end));
    gcInputs.hInlet = h;
    gcInputs.DmInlet = Dm;
    evInputs.DmOutlet = Dm;
    % Valve
    [Dm, h] = valve.process(cR,dInlet,pInlet,pOutlet,hInlet);
    evInputs.hInlet = h;
    evInputs.DmInlet = Dm;
    gcInputs.DmOutlet = Dm;
    % Gas cooler
    gc.timestep([(it-1)*timestep it*timestep],gcInputs);
    % Evaporator
    ev.timestep([(it-1)*timestep it*timestep],evInputs);
end
toc
