classdef Receiver   < matlab.mixin.Copyable % TODO
    properties
        % Parameters
        Volume      % Total VLE volume 
        % Variables
        Dm                      % Mass flow rate
        p                       % Pressure
        h                       % Enthalpy
        d                       % Density
        Dp                      % Pressure derivative
        Dh                      % Enthalpy derivative
        Dd                      % Density derivative
        record                  % Record of results
        ODEoptions              % Options of ODE solver
        t                       % Last time instant
        % Partial derivatives
        Dd_Dh
        Dd_Dp
        % Outlet states
        hGas
        hLiquid
        dGas
        dLiquid
    end
    methods
        function massAccummulation(rec,DmInlet,DmGas,DmLiquid)
            rec.Dd = (DmInlet - DmGas - DmLiquid)/rec.Volume;
        end
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
        function potentialAccummulation(rec,hInlet,DmInlet,DmGas,DmLiquid)
            Dpsi = [DmInlet -DmGas -DmLiquid]*...
                [hInlet; rec.hGas; rec.hLiquid]/rec.Volume;
            try
                Dd_Dp = CoolProp.PropsSI('d(D)/d(P)|H','H',rec.h,'P',rec.p,'CO2');
                Dd_Dh = CoolProp.PropsSI('d(D)/d(H)|P','H',rec.h,'P',rec.p,'CO2');
                rec.Dd_Dp = Dd_Dp;
                rec.Dd_Dh = Dd_Dh;
                global partials
                partials = [partials [Dd_Dp; Dd_Dh]];
            catch
                global bugnumber
                bugnumber = bugnumber+1;
            end
            DpDh_vector = [-1 rec.d; rec.Dd_Dp rec.Dd_Dh]\[Dpsi; rec.Dd];
            rec.Dp = DpDh_vector(1,1);
            rec.Dh = DpDh_vector(2,1);
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
            rec.massAccummulation(DmInlet,DmGas,DmLiquid);
            rec.separation();
            rec.potentialAccummulation(hInlet,DmInlet,DmGas,DmLiquid);
            Dx = [rec.Dp; rec.Dh; rec.Dd];
        end
        function initialize(rec,p,h,Volume,ODEoptions)
            % Function help: provide initial pressure, enthalpy, and the
            %   volume

            rec.Volume = Volume;
            rec.p = p;
            rec.h = h;
            rec.ph2d;
            rec.ODEoptions = ODEoptions;
            Record.t = [];
            Record.x = [];
            rec.record = Record;
        end
        function reinitialize(rec,p,h)
            rec.p = p;
            rec.h = h;
            rec.ph2d;
            rec.record.t = [];
            rec.record.x = [];
        end
        function ph2d(rec)
            rec.d = CoolProp.PropsSI('D','P',rec.p,'H',rec.h,'CO2');
        end
        function dh2p(rec)
            rec.p = CoolProp.PropsSI('P','D',rec.d,'H',rec.h,'CO2');
        end
    end
end