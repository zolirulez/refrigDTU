classdef Receiver   < Tank
    properties
        % Outlet states
        gas         % Structure with fields p,h,d,DmOutlet
        liquid      % Structure with fields p,h,d,DmOutlet
        DmInlet
        hInlet
        DmOutlet
    end
    methods
        function separation(rec)
            rec.gas.h = CoolProp.PropsSI('H','P',rec.p,'Q',1,'CO2');
            rec.liquid.h = CoolProp.PropsSI('H','P',rec.p,'Q',0,'CO2');
            rec.gas.d = CoolProp.PropsSI('D','P',rec.p,'H',rec.gas.h,'CO2');
            rec.liquid.d = CoolProp.PropsSI('D','P',rec.p,'H',rec.liquid.h,'CO2');
            rec.gas.p = rec.p;
            rec.liquid.p = rec.p;
            if rec.gas.h < rec.h
                disp('Error: receiver is empty of liquid')
            elseif rec.liquid.h > rec.h
                disp('Error: receiver is full of liquid')
            end
        end
        function Dx = process(rec,t,x)
            % Time
            rec.t = t;
            % States
            rec.p = x(1,1);
            rec.h = x(2,1);
            rec.d = x(3,1);
            % Inputs
            rec.DmOutlet = rec.gas.DmOutlet + rec.liquid.DmOutlet;
            % Process
            rec.massAccummulation();
            rec.separation();
            rec.excitation([rec.DmInlet; -rec.gas.DmOutlet; -rec.liquid.DmOutlet],...
                [rec.hInlet; rec.gas.h; rec.liquid.h],0);
            rec.potentialAccummulation();
            Dx = [rec.Dp; rec.Dh; rec.Dd];
        end
        function initialize(rec,p,h,Volume,ODEoptions)
            % Function help: provide initial pressure, enthalpy, and the
            %   volume, and calculates separations

            initialize@Tank(rec,p,h,Volume,ODEoptions);
            rec.separation();
        end
    end
end