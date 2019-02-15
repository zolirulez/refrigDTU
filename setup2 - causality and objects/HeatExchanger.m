classdef HeatExchanger < matlab.mixin.Copyable
    properties
        % Parameters
        OneTubeTotalVolume      % Total VLE volume of one tube
        OneTubeCellVolume       % Volume of a cell
        InnerTubeDiameter
        OneTubeLength           % Length for the tube for one flow
        nParallelFlows          % Number of parallel flows
        f1                      % Friction factor
        OneTubeTotalResistance  % Total resistance of one tube
        OneCellResistance       % Resistance of one piece of n parallel tubes
        % Discretization parameter
        nCell
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
        % Boundaries
        DmInlet
        DmOutlet
        hInlet
        % Heat Exchange
        T
        Ta
        NominalThermalResistance
        NominalVolumeFlow
        ConductionSlope
        OneCellConductionSlope
        NaturalConduction             
        OneCellNaturalConduction
        oneCellThermalResistance
        DVInlet
        TInlet
        cp                   
        da         
        Qnominal
        dTlog
        DQ
    end
    methods
        function massflow(hx)
            deltap = hx.p(1:end-1)-hx.p(2:end);
            inducedDm = sign(deltap).*...
                sqrt(abs(deltap.*hx.d(1:end-1)))/hx.OneCellResistance;
%             inducedDm = -diff(hx.h-hx.p./hx.d)./hx.OneCellResistance;

            Dm1  = [hx.DmInlet; inducedDm; hx.DmOutlet];
            Dm2  = linspace(hx.DmInlet,hx.DmOutlet,hx.nCell+1)';
            hx.Dm = (Dm1*0.5 + Dm2*0.5);
            hx.Dm(1) = hx.DmInlet;
            hx.Dm(end) = hx.DmOutlet;
        end
        function massAccummulation(hx)
            hx.Dd = -diff(hx.Dm)/hx.OneTubeCellVolume;
        end
        function potentialAccummulation(hx)
             for it = 1:hx.nCell
                try
                    hx.ph2T;
                catch
                    global bugnumber
                    bugnumber = bugnumber+1;
                end
             end
             hx.oneCellThermalResistance =...
                 1/(hx.OneCellNaturalConduction+hx.OneCellConductionSlope*hx.DVInlet)...
                 + 1/(hx.da*hx.DVInlet*hx.cp);
             weightFactor = hx.oneCellThermalResistance*hx.da*hx.DVInlet*hx.cp;
             hx.Ta = 1/(weightFactor + 1)*flip(hx.T) +...
                 1/(1/weightFactor + 1)*[hx.TInlet; hx.Ta(1:end-1)];
             hx.DQ = (flip(hx.Ta) - hx.T)/hx.oneCellThermalResistance;
             Dpsi = (-diff(hx.Dm .*[hx.hInlet; hx.h(1:end)]) + hx.DQ)/hx.OneTubeCellVolume;
             DpDh_vector = zeros(2,hx.nCell);
             for it = 1:hx.nCell
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
                Dd = hx.Dd; %sign(hx.Dd).*sqrt(abs(hx.Dd));
                DpDh_vector(:,it) = [-1 hx.d(it); hx.Dd_Dp hx.Dd_Dh]\...
                    [Dpsi(it)-hx.p(it)/hx.d(it)*Dd(it); Dd(it)]; % TODO
            end
            hx.Dp = DpDh_vector(1,:)';
            hx.Dh = DpDh_vector(2,:)';
%             for it = 1:hx.nCell
%                 hx.Dp(it) = sign(hx.Dp(it))*min(abs(hx.Dp(it)),1e7);
%                 hx.Dh(it) = sign(hx.Dh(it))*min(abs(hx.Dh(it)),1e5);
%             end
        end
        function timestep(hx,t,inputs)
            % Function help: simulates the heat exchanger from time t1 to
            %   time t2. Input is a structure that has the following
            %   fields: DmInlet, DmOutlet and the heatflow
            
            x = [hx.p; hx.h; hx.d; hx.pState];
            [t, x] = ode15s(@hx.process,[t(1) t(2)],x,hx.ODEoptions,inputs);
            hx.record.t = [hx.record.t; t];
            hx.record.x = [hx.record.x; x];
            hx.p = x(end,1:hx.nCell)';
            hx.h = x(end,hx.nCell+1:2*hx.nCell)';
            hx.d = x(end,2*hx.nCell+1:3*hx.nCell)';
        end
        function Dx = process(hx,t,x)
            % Time
            hx.t = t;
            % States
            hx.p = x(1:hx.nCell,1);
            hx.h = x(hx.nCell+1:2*hx.nCell,1);
            hx.d = x(2*hx.nCell+1:3*hx.nCell,1);
            % Process
            massflow(hx);
            massAccummulation(hx);
            potentialAccummulation(hx);
            Dx = [hx.Dp; hx.Dh; hx.Dd];
        end
        function initialize(hx,nCell,p,h,Parameters)
            % Function help: provide two-element vectors for pressure and 
            %   enthalpy inlets and outlets, and provide volume.

            hx.nCell = nCell;
            % Resistance
            hx.InnerTubeDiameter = Parameters.InnerTubeDiameter;
            hx.OneTubeLength = Parameters.OneTubelength;
            hx.nParallelFlows = Parameters.nParallelFlows;
            hx.f1 = Parameters.f1;
            hx.OneTubeTotalVolume =...
                hx.InnerTubeDiameter^2*pi/4*hx.OneTubeLength;
            hx.OneTubeTotalResistance = sqrt(16*hx.f1*hx.OneTubeLength/...
                (pi^2*hx.InnerTubeDiameter^5));
            % States
            hx.p = linspace(p(1),p(2),nCell)';
            hx.h = linspace(h(1),h(2),nCell)';
            hx.ph2d;
            hx.ph2T;
            % Thermal properties
            hx.thermalInitialization(Parameters);
            % Discretizing
            hx.discretize();
            % Simulation
            Record.t = [];
            Record.x = [];
            hx.record = Record;
        end
        function thermalInitialization(hx,Parameters)
            % Thermal properties
            hx.Qnominal = Parameters.Qnominal;
            dTi = Parameters.Tpi - Parameters.Tso;
            dTo = Parameters.Tpo - Parameters.Tsi;
            hx.dTlog = (dTo - dTi)/log(dTo/dTi); 
            CorrectionFactor = 2;
            hx.NominalThermalResistance = hx.dTlog/hx.Qnominal/CorrectionFactor;
            hx.NominalVolumeFlow = Parameters.NominalVolumeFlow;
            hx.ConductionSlope = (1 - Parameters.ConductionRatio)/...
                hx.NominalThermalResistance/hx.NominalVolumeFlow;
            hx.cp = 1000;
            hx.da = 1.25;
            hx.NaturalConduction =...
                Parameters.ConductionRatio/hx.NominalThermalResistance;
        end
        function reinitialize(hx,p,h)
            hx.p = linspace(p(1),p(2),hx.nCell)';
            hx.h = linspace(h(1),h(2),hx.nCell)';
            hx.ph2d;
            hx.ph2T;
            hx.record.t = [];
            hx.record.x = [];
        end
        function discretize(hx)
            hx.OneTubeCellVolume = hx.OneTubeTotalVolume/hx.nCell;
            hx.OneCellResistance = sqrt(hx.OneTubeTotalResistance^2/...
                hx.nCell/hx.nParallelFlows);
            hx.OneCellConductionSlope = hx.ConductionSlope/hx.nCell;
            hx.OneCellNaturalConduction = hx.NaturalConduction/hx.nCell;
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
        function ph2T(hx)
            hx.T = zeros(hx.nCell,1);
            for it = 1:hx.nCell
                hx.T(it) = CoolProp.PropsSI('T','P',hx.p(it),'H',hx.h(it),'CO2');
            end
        end
    end
end