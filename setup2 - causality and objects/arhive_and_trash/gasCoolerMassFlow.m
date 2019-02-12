clearvars
givenVolume = 0.0192;
InnerTubeDiameter = 0.0192;
nParallelFlows = 5;
OneTubelength = givenVolume/(InnerTubeDiameter^2*pi/4);

T = 107.6+273.15;
d = CoolProp.PropsSI('D','T',T,'P',84.8e5,'CO2');
Dm = d*7.87/3600;
atmero = 0.0192;

% 1 and 2 indexes mean liquid min and vapour max
mu1 = 156e-6;
mu2 = 12.5e-6;
Dm1 = Dm/nParallelFlows
Dm2 = Dm/nParallelFlows
deltap = 1.57e5;
L = OneTubelength;

for i = 1:100
    % Reynolds
    Re1 = 4*Dm1/pi/atmero/mu1;
    Re2 = 4*Dm2/pi/atmero/mu2;
    % Gnielinki
    f1 = 0.5/(0.79*log(Re1)-1.64)^2;
    f2 = 0.5/(0.79*log(Re2)-1.64)^2;
    % mass flow in one tube
    Dm1 = sqrt(deltap/f1*pi*pi*atmero^5*d/16/L);
    Dm2 = sqrt(deltap/f2*pi*pi*atmero^5*d/16/L);
end
Dm1
Dm2