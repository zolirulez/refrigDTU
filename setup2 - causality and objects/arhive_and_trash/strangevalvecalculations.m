clearvars
% Calculations for HP valve

pRec = 40e5;
dp = [linspace(0.1,2,1000) linspace(2,50e5,100)];
mdotSmooth = 0.0005;
effectiveFlowAreaTypical = 0.2e-6;
f = 2;

d = zeros(1,length(dp));
mdotdens = zeros(1,length(dp));
dpsmooth = zeros(1,length(dp));
for i=1:length(dp)
    d(i) = CoolProp.PropsSI('D','P',pRec+dp(i),'H',300e3,'CO2');
    fd = f*d(i);
    dpsmooth(i) = mdotSmooth/effectiveFlowAreaTypical/fd;
    a = -1/(4*(dpsmooth(i)*fd)^(5/2));
    b = 5/(4*sqrt(dpsmooth(i)*fd));
    if (dp(i) > dpsmooth)
        mdotdens(i) = sqrt(fd*dp(i));
    else
        mdotdens(i) = (a*fd*dp(i)*fd*dp(i) + b)*fd*dp(i);
    end
end
figure(1)
subplot(311)
loglog(dp,d)
title(['Receiver pressure: ' num2str(pRec)])
xlabel('dp')
ylabel('Density')
subplot(312)
semilogx(dp,dpsmooth)
xlabel('dp')
ylabel('dp_{smooth}')
subplot(313)
loglog(dp,mdotdens)
xlabel('dp')
ylabel('Mass flow rate density [mdot/A]')