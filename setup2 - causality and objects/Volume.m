classdef Volume < matlab.mixin.Copyable 
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
        DPsi                    % Excitation
        volord                  % volord of results
        ODEoptions              % Options of ODE solver
        t                       % Last time instant
        % Partial derivatives
        Dd_Dh
        Dd_Dp
    end
    methods
        function massAccummulation(vol,DmInlet,DmOutlet)
            vol.Dd = (DmInlet - DmOutlet)/vol.Volume;
        end
        function excitation(vol,hInlet,DmInlet,DmOutlet)
            vol.Dpsi = [DmInlet -DmOutlet]*[hInlet; vol.h]/vol.Volume;
        end
        function potentialAccummulation(vol)
            try
                Dd_Dp = CoolProp.PropsSI('d(D)/d(P)|H','H',vol.h,'P',vol.p,'CO2');
                Dd_Dh = CoolProp.PropsSI('d(D)/d(H)|P','H',vol.h,'P',vol.p,'CO2');
                vol.Dd_Dp = Dd_Dp;
                vol.Dd_Dh = Dd_Dh;
                global partials
                partials = [partials [Dd_Dp; Dd_Dh]];
            catch
                global bugnumber
                bugnumber = bugnumber+1;
            end
            DpDh_vector = [-1 vol.d; vol.Dd_Dp vol.Dd_Dh]\[vol.Dpsi; vol.Dd];
            vol.Dp = DpDh_vector(1,1);
            vol.Dh = DpDh_vector(2,1);
        end
        function timestep(vol,t,inputs)
            % Function help: 
            
            x = [vol.p; vol.h; vol.d;];
            [t, x] = ode15s(@vol.process,[t(1) t(2)],x,vol.ODEoptions,inputs);
            vol.volord.t = [vol.volord.t; t];
            vol.volord.x = [vol.volord.x; x];
            vol.p = x(end,1)';
            vol.h = x(end,2)';
            vol.d = x(end,3)';
        end
        function Dx = process(vol,t,x,Inputs)
            % Time
            vol.t = t;
            % States
            vol.p = x(1,1);
            vol.h = x(2,1);
            vol.d = x(3,1);
            % Inputs
            DmInlet = Inputs.DmInlet;
            DmLiquid = Inputs.DmLiquid;
            DmGas = Inputs.DmGas;
            hInlet = Inputs.hInlet;
            % Process
            vol.massAccummulation(DmInlet,DmGas,DmLiquid);
            vol.excitation();
            vol.potentialAccummulation(hInlet,DmInlet,DmGas,DmLiquid);
            Dx = [vol.Dp; vol.Dh; vol.Dd];
        end
        function initialize(vol,p,h,Volume,ODEoptions)
            % Function help: provide initial pressure, enthalpy, and the
            %   volume

            vol.Volume = Volume;
            vol.p = p;
            vol.h = h;
            vol.ph2d;
            vol.ODEoptions = ODEoptions;
            Record.t = [];
            Record.x = [];
            vol.record = Record;
        end
        function reinitialize(vol,p,h)
            vol.p = p;
            vol.h = h;
            vol.ph2d;
            vol.record.t = [];
            vol.record.x = [];
        end
        function ph2d(vol)
            vol.d = CoolProp.PropsSI('D','P',vol.p,'H',vol.h,'CO2');
        end
        function dh2p(vol)
            vol.p = CoolProp.PropsSI('P','D',vol.d,'H',vol.h,'CO2');
        end
    end
end