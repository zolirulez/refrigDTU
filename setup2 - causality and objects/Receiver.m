classdef Receiver   < Tank
    properties
        % Outlet states
        hGas
        hLiquid
        dGas
        dLiquid
    end
    methods
        function separation(rec)
            rec.hGas = CoolProp.PropsSI('H','P',rec.p,'Q',1,'CO2');
            rec.hLiquid = CoolProp.PropsSI('H','P',rec.p,'Q',0,'CO2');
            rec.dGas = CoolProp.PropsSI('D','P',rec.p,'H',rec.hGas,'CO2');
            rec.dLiquid = CoolProp.PropsSI('D','P',rec.p,'H',rec.hLiquid,'CO2');
            if rec.hGas < rec.h
                disp('Error: receiver is empty of liquid')
            elseif rec.hLiquid > rec.h
                disp('Error: receiver is full of liquid')
            end
        end
        function excitation(rec,hInlet,DmInlet,DmGas,DmLiquid)
            rec.Dpsi = [DmInlet -DmGas -DmLiquid]*...
                [hInlet; rec.hGas; rec.hLiquid]/rec.Volume;
        end
        function timestep(rec,t,inputs)
            % Function help: 
            
            x = [rec.p; rec.h; rec.d;];
            [t, x] = ode15s(@rec.process,[t(1) t(2)],x,rec.ODEoptions,inputs);
            rec.record.t = [rec.record.t; t];
            rec.record.x = [rec.record.x; x];
            rec.p = x(end,1)';
            rec.h = x(end,2)';
            rec.d = x(end,3)';
        end
        function Dx = process(rec,t,x,Inputs)
            % Time
            rec.t = t;
            % States
            rec.p = x(1,1);
            rec.h = x(2,1);
            rec.d = x(3,1);
            % Inputs
            DmInlet = Inputs.DmInlet;
            DmLiquid = Inputs.DmLiquid;
            DmGas = Inputs.DmGas;
            hInlet = Inputs.hInlet;
            % Process
            rec.massAccummulation(DmInlet,DmGas+DmLiquid);
            rec.excitation(hInlet,DmInlet,DmGas,DmLiquid);
            rec.separation();
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