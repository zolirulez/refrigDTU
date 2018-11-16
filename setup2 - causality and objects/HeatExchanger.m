classdef HeatExchanger < handle
    properties
        OneTubeTotalVolume      % Total VLE volume of one tube
        OneTubeCellVolume       % Volume of a cell
        InnerTubeDiameter
        OneTubeLength           % Length for the tube for one flow
        nParallelFlows          % Number of parallel flows
        f1                      % Friction factor
        OneTubeTotalResistance  % Total resistance of one tube
        OneTubeCellResistance   % Resistance of one piece of one tube
%         stiffnessMatrix
%         invStiffnessMatrix
        nCell
        Dm 
        p
        h
        d
        pState
        hState
        dState
        pStateTimeConstant
        hStateTimeConstant
        dStateTimeConstant
        Dp
        Dh
        Dd
        DL
        DpState
        DhState
        DdState
        record
        ODEoptions
        t
    end
    methods
        function discretize(hx)
            hx.OneTubeCellVolume = hx.OneTubeTotalVolume/hx.nCell;
            hx.OneTubeCellResistance = sqrt(hx.OneTubeTotalResistance^2/hx.nCell);
%             *** If FEM was needed: ***
% 
%             hx.stiffnessMatrix = zeros(nCell);
%             for it = 1:nCell-1
%                 hx.stiffnessMatrix(it:it+1,it:it+1) = ...
%                     hx.stiffnessMatrix(it:it+1,it:it+1) + [1 -1; -1 1];
%             end
%             hx.invStiffnessMatrix = inv(hx.stiffnessMatrix);
        end
        function massflow(hx,DmInlet,DmOutlet)
            deltap = hx.p(1:end-1)-hx.p(2:end);
            inducedDm = sign(deltap).*...
                sqrt(abs(deltap.*(0.5*(hx.d(1:end-1)+hx.d(2:end)))))...
                /hx.OneTubeCellResistance;
            hx.Dm  = [DmInlet/hx.nParallelFlows; inducedDm;...
                DmOutlet/hx.nParallelFlows];
        end
        function massAccummulation(hx)
            hx.Dd = -diff(hx.Dm)/hx.OneTubeCellVolume;
%             hx.Dd = (hx.dState-hx.d)/hx.dStateTimeConstant;
        end
        function potentialAccummulation(hx,hInlet,DQ)
            Dpsi = (-diff(hx.Dm .*...
                [hInlet; 0.5*(hx.h(1:end-1)+hx.h(2:end)); hx.h(end)])... 298e3 put this back TODO
                + DQ)/hx.OneTubeCellVolume;
            DpDh_vector = zeros(2,hx.nCell);
            for it = 1:hx.nCell
                Dd_Dp = CoolProp.PropsSI('d(D)/d(P)|H','H',hx.h(it),'P',hx.p(it),'CO2');
                Dd_Dh = CoolProp.PropsSI('d(D)/d(H)|P','H',hx.h(it),'P',hx.p(it),'CO2');
                DpDh_vector(:,it) = [-1 hx.d(it); Dd_Dp Dd_Dh]\[Dpsi(it); hx.Dd(it)];
            end
            hx.DpState = DpDh_vector(1,:)';
            hx.Dh = DpDh_vector(2,:)';
            hx.Dp = (hx.pState-hx.p)/hx.pStateTimeConstant;
%             hx.Dh = (hx.hState-hx.h)/hx.hStateTimeConstant;
        end
        function timestep(hx,t,inputs)
            % Function help: simulates the heat exchanger from time t1 to
            %   time t2. Input is a structure that has the following
            %   fields: DmInlet, DmOutlet and the heatflow
            
            x = [hx.p; hx.h; hx.d; hx.pState]; %; hx.pState; hx.hState; hx.dState];
            [t, x] = ode15s(@hx.process,[t(1) t(2)],x,hx.ODEoptions,inputs);
            hx.record.t = [hx.record.t; t];
            hx.record.x = [hx.record.x; x];
            hx.p = x(end,1:hx.nCell)';
            hx.h = x(end,hx.nCell+1:2*hx.nCell)';
            hx.d = x(end,2*hx.nCell+1:3*hx.nCell)';
            hx.pState = x(end,3*hx.nCell+1:4*hx.nCell)';
%             hx.hState = x(1,4*hx.nCell+1:5*hx.nCell)';
%             hx.dState = x(1,5*hx.nCell+1:6*hx.nCell)';
        end
        function Dx = process(hx,t,x,inputs)
            % Time
            hx.t = t;
            % States
            hx.p = x(1:hx.nCell,1);
            hx.h = x(hx.nCell+1:2*hx.nCell,1);
            hx.d = x(2*hx.nCell+1:3*hx.nCell,1);
            hx.pState = x(3*hx.nCell+1:4*hx.nCell,1);
%             hx.hState = x(4*hx.nCell+1:5*hx.nCell,1);
%             hx.dState = x(5*hx.nCell+1:6*hx.nCell,1);
            % hx.ph2d;
            % Inputs
            DmInlet = inputs.DmInlet;
            DmOutlet = inputs.DmOutlet;
            hInlet = inputs.hInlet;
            DQ = inputs.DQ;
            % Process
            massflow(hx,DmInlet,DmOutlet);
            massAccummulation(hx);
            potentialAccummulation(hx,hInlet,DQ);
%             Dx = [hx.Dp; hx.Dh; hx.Dd; hx.DpState; hx.DhState; hx.DdState];
            Dx = [hx.DpState; hx.Dh; hx.Dd; hx.DpState];
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
            hx.hState = hx.h;
            hx.hStateTimeConstant = TimeConstants(2);
            hx.ph2d;
            hx.dState = hx.d;
            hx.dStateTimeConstant = TimeConstants(3);
            hx.ODEoptions = ODEoptions;
            record.t = [];
            record.x = [];
            hx.record = record;
        end
        function ph2d(hx)
            hx.d = zeros(hx.nCell,1);
            for it = 1:hx.nCell
                hx.d(it) = CoolProp.PropsSI('D','P',hx.p(it),'H',hx.h(it),'CO2');
            end
        end
    end
end