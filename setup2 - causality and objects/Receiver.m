classdef Receiver   < matlab.mixin.Copyable % TODO
    properties
        % Parameters
        Volume      % Total VLE volume 
        % Variables
        Dm                      % Mass flow rate
        p                       % Pressure
        h                       % Enthalpy
        d                       % Density
        pState                  % Pressure state for smoothing
        pStateTimeConstant      % First order time constant for pState
        Dp                      % Pressure derivative
        Dh                      % Enthalpy derivative
        Dd                      % Density derivative
        DpState                 % Pressure state derivative
        record                  % Record of results
        ODEoptions              % Options of ODE solver
        t                       % Last time instant
        % Partial derivatives
        Dd_Dh
        Dd_Dp
    end
    methods
        function massAccummulation(hx)
            hx.Dd = -diff(hx.Dm)/hx.OneTubeCellVolume;
        end
        function separation(hx)
        end
        function potentialAccummulation(hx,hInlet)
            Dpsi = -diff(hx.Dm .*[hInlet; hx.h(1:end)])/hx.OneTubeCellVolume;
            DpDh_vector = zeros(2,hx.nCell);
            try
                Dd_Dp = CoolProp.PropsSI('d(D)/d(P)|H','H',hx.h(it),'P',hx.p(it),'CO2');
                Dd_Dh = CoolProp.PropsSI('d(D)/d(H)|P','H',hx.h(it),'P',hx.p(it),'CO2');
                hx.Dd_Dp = Dd_Dp;
                hx.Dd_Dh = Dd_Dh;
                global partials
                partials = [partials [Dd_Dp; Dd_Dh]];
            catch
                global bugnumber
                bugnumber = bugnumber+1;
            end
            DpDh_vector(:,it) = [-1 hx.d(it); hx.Dd_Dp hx.Dd_Dh]\...
                [Dpsi(it); hx.Dd(it)];
            hx.Dp = DpDh_vector(1,:)';
            hx.Dh = DpDh_vector(2,:)';
        end
        function timestep(hx,t,inputs)
            % Function help: 
            
            x = [hx.p; hx.h; hx.d; hx.pState];
            [t, x] = ode15s(@hx.process,[t(1) t(2)],x,hx.ODEoptions,inputs);
            hx.record.t = [hx.record.t; t];
            hx.record.x = [hx.record.x; x];
            hx.p = x(end,1:hx.nCell)';
            hx.h = x(end,hx.nCell+1:2*hx.nCell)';
            hx.d = x(end,2*hx.nCell+1:3*hx.nCell)';
        end
        function Dx = process(hx,t,x,inputs)
            % Time
            hx.t = t;
            % States
            hx.p = x(1:hx.nCell,1);
            hx.h = x(hx.nCell+1:2*hx.nCell,1);
            hx.d = x(2*hx.nCell+1:3*hx.nCell,1);
            %hx.pState = x(3*hx.nCell+1:4*hx.nCell,1);
            % hx.dh2p;
            % Inputs
            DmInlet = inputs.DmInlet;
            DmOutlet = inputs.DmOutlet;
            hInlet = inputs.hInlet;
            DQ = inputs.DQ;
            % Process
            massAccummulation(hx);
            potentialAccummulation(hx,hInlet,DQ);
            hx.separation;
            Dx = [hx.Dp; hx.Dh; hx.Dd];%[hx.DpState; hx.Dh; hx.Dd; hx.DpState];
        end
        function initialize(hx,nCell,p,h,Parameters,TimeConstants,ODEoptions)
            % Function help: provide two-element vectors for pressure and 
            %   enthalpy inlets and outlets, and provide volume.

            hx.nCell = nCell;
            hx.InnerTubeDiameter = Parameters.InnerTubeDiameter;
            hx.OneTubeLength = Parameters.OneTubelength;
            hx.nParallelFlows = Parameters.nParallelFlows;
            hx.f1 = Parameters.f1;
            hx.OneTubeTotalVolume =...
                hx.InnerTubeDiameter^2*pi/4*hx.OneTubeLength;
            hx.OneTubeTotalResistance = sqrt(16*hx.f1*hx.OneTubeLength/...
                (pi^2*hx.InnerTubeDiameter^5));
            hx.discretize();
            hx.p = linspace(p(1),p(2),nCell)';
            hx.pState = hx.p;
            hx.pStateTimeConstant = TimeConstants(1);
            hx.h = linspace(h(1),h(2),nCell)';
            hx.ph2d;
            hx.ODEoptions = ODEoptions;
            Record.t = [];
            Record.x = [];
            hx.record = Record;
        end
        function reinitialize(hx,p,h)
            hx.p = linspace(p(1),p(2),hx.nCell)';
            hx.pState = hx.p;
            hx.pStateTimeConstant = TimeConstants(1);
            hx.h = linspace(h(1),h(2),hx.nCell)';
            hx.ph2d;
            hx.record.t = [];
            hx.record.x = [];
        end
        function ph2d(hx)
            hx.d = zeros(hx.nCell,1);
            for it = 1:hx.nCell
                hx.d(it) = CoolProp.PropsSI('D','P',hx.p(it),'H',hx.h(it),'CO2');
            end
        end
        function dh2p(hx)
            hx.p = zeros(hx.nCell,1);
            for it = 1:hx.nCell
                hx.p(it) = CoolProp.PropsSI('P','D',hx.d(it),'H',hx.h(it),'CO2');
            end
        end
    end
end