classdef Tank < matlab.mixin.Copyable 
    properties
        % Parameters
        Volume      % Total VLE volume 
        % Variables
        p                       % Pressure
        h                       % Enthalpy
        d                       % Density
        Dp                      % Pressure derivative
        Dh                      % Enthalpy derivative
        Dd                      % Density derivative
        Dpsi                    % Excitation
        record                  % Record of results
        ODEoptions              % Options of ODE solver
        t                       % Last time instant
        % Partial derivatives
        Dd_Dh
        Dd_Dp
    end
    methods
        function massAccummulation(tank)
            tank.Dd = (tank.DmInlet - tank.DmOutlet)/tank.Volume;
        end
        function excitation(tank,Dm,h,DQ)
            tank.Dpsi = (Dm'*h + DQ)/tank.Volume;
        end
        function potentialAccummulation(tank)
            try
                Dd_Dp = CoolProp.PropsSI('d(D)/d(P)|H','H',tank.h,'P',tank.p,'CO2');
                Dd_Dh = CoolProp.PropsSI('d(D)/d(H)|P','H',tank.h,'P',tank.p,'CO2');
                tank.Dd_Dp = Dd_Dp;
                tank.Dd_Dh = Dd_Dh;
                global partials
                partials = [partials [Dd_Dp; Dd_Dh]];
            catch
                global bugnumber
                bugnumber = bugnumber+1;
            end
            DpDh_vector = [-1 tank.d; tank.Dd_Dp tank.Dd_Dh]\[tank.Dpsi; tank.Dd];
            tank.Dp = DpDh_vector(1,1);
            tank.Dh = DpDh_vector(2,1);
        end
        function timestep(tank,t,inputs)
            % Function help: 
            
            x = [tank.p; tank.h; tank.d];
            [t, x] = ode15s(@tank.process,[t(1) t(2)],x,tank.ODEoptions,inputs);
            tank.record.t = [tank.record.t; t];
            tank.record.x = [tank.record.x; x];
            tank.p = x(end,1)';
            tank.h = x(end,2)';
            tank.d = x(end,3)';
        end
        function Dx = process(tank,t,x,DQ)
            % Time
            tank.t = t;
            % States
            tank.p = x(1,1);
            tank.h = x(2,1);
            tank.d = x(3,1);
            % Process
            tank.massAccummulation();
            tank.excitation([tank.DmInlet; -tank.DmOutlet],...
                [tank.hInlet; tank.h*ones(size(tank.DmOutlet))],DQ);
            tank.potentialAccummulation();
            Dx = [tank.Dp; tank.Dh; tank.Dd];
        end
        function initialize(tank,p,h,Volume,ODEoptions)
            % Function help: provide initial pressure, enthalpy, and the
            %   volume

            tank.Volume = Volume;
            tank.p = p;
            tank.h = h;
            tank.ph2d;
            tank.ODEoptions = ODEoptions;
            Record.t = [];
            Record.x = [];
            tank.record = Record;
        end
        function reinitialize(tank,p,h)
            tank.p = p;
            tank.h = h;
            tank.ph2d;
            tank.record.t = [];
            tank.record.x = [];
        end
        function ph2d(tank)
            tank.d = CoolProp.PropsSI('D','P',tank.p,'H',tank.h,'CO2');
        end
        function dh2p(tank)
            tank.p = CoolProp.PropsSI('P','D',tank.d,'H',tank.h,'CO2');
        end
    end
end