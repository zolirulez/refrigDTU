figure(1)
clf
pvector = linspace(14e5,85e5,4);
for it3 = 1:4
    p = pvector(it3);
    h = linspace(212e3,525e3,100);
    Dd = linspace(-100,100,10);
    d = zeros(100,10);
    Dp = zeros(100,10);
    Dh = zeros(100,10);
    for it2 = 1:10
        for it1 = 1:100
            d(it1,it2) = CoolProp.PropsSI('D','H',h(it1),'P',p,'CO2');
            DhDd = CoolProp.PropsSI('d(H)/d(D)|P','H',h(it1),'P',p,'CO2');
            DhDp = CoolProp.PropsSI('d(H)/d(P)|D','H',h(it1),'P',p,'CO2');
            Dp(it1,it2) = 1/(d(it1)*DhDp-1)*(d(it1)*DhDd+p/d(it1))*Dd(it2); % no excitation
            Dh(it1,it2) = DhDp*Dp(it1,it2)+DhDd*Dd(it2);
        end
        subplot(4,2,it3*2-1)
        hold on
        plot(h,Dp(:,it2))
        hold off
        subplot(4,2,it3*2)
        hold on
        plot(h,Dh(:,it2))
        hold off
    end
end
xlabel('enthalpy')