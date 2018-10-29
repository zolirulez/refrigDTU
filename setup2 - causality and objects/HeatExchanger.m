classdef HeatExchanger < handle
    properties
        totalVolume
        cellVolume
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
        DpState
        DhState
        DdState
        record
        ODEoptions
        t
    end
    methods
        function discretize(hx)
            hx.cellVolume = hx.totalVolume/hx.nCell;
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
            averageResistance = abs(hx.p(1)-hx.p(end)) /...
                ((DmInlet+DmOutlet)/2)/(hx.nCell-1);
            hx.Dm  = [DmInlet;...
                (hx.p(1:end-1)-hx.p(2:end))/averageResistance;...
                DmOutlet];
%             figure(2)
%             plot(1:hx.nCell,hx.Dm)
%             title(num2str(averageResistance))
%             pause(0.5)
        end
        function massAccummulation(hx)
            hx.DdState = -diff(hx.Dm)/hx.cellVolume;
            hx.Dd = (hx.dState-hx.d)/hx.dStateTimeConstant;
        end
        function potentialAccummulation(hx,hInlet,DQ)
            Dpsi = (-diff(hx.Dm .*...
                [hInlet; 0.5*(hx.h(1:end-1)+hx.h(2:end)); hx.h(end)])...
                + DQ)/hx.cellVolume;
            DpDh_vector = zeros(2,hx.nCell);
            for it = 1:hx.nCell
                dddp = CoolProp.PropsSI('d(D)/d(P)|H','H',hx.h(it),'P',hx.p(it),'CO2');
                dddh = CoolProp.PropsSI('d(D)/d(H)|P','H',hx.h(it),'P',hx.p(it),'CO2');
                DpDh_vector(:,it) = [-1 hx.d(it); dddp dddh]\[Dpsi(it); hx.Dd(it)];
            end
            hx.DpState = DpDh_vector(1,:)';
            hx.DhState = DpDh_vector(2,:)';
            hx.Dp = (hx.pState-hx.p)/hx.pStateTimeConstant;
            hx.Dh = (hx.hState-hx.h)/hx.hStateTimeConstant;
        end
        function timestep(hx,t,inputs)
            % Function help: simulates the heat exchanger from time t1 to
            %   time t2. Input is a structure that has the following
            %   fields: DmInlet, DmOutlet and the heatflow
            
            x = [hx.p; hx.h; hx.d; hx.pState; hx.hState; hx.dState];
            [t, x] = ode15s(@hx.process,[t(1) t(2)],x,hx.ODEoptions,inputs);
            hx.record.t = [hx.record.t; t];
            hx.record.x = [hx.record.x; x];
            hx.p = x(1,1:hx.nCell)';
            hx.h = x(1,hx.nCell+1:2*hx.nCell)';
            hx.d = x(1,2*hx.nCell+1:3*hx.nCell)';
            hx.pState = x(1,3*hx.nCell+1:4*hx.nCell)';
            hx.hState = x(1,4*hx.nCell+1:5*hx.nCell)';
            hx.dState = x(1,5*hx.nCell+1:6*hx.nCell)';
        end
        function Dx = process(hx,t,x,inputs)
            % Time
            hx.t = t;
            % States
            hx.p = x(1:hx.nCell,1);
            hx.h = x(hx.nCell+1:2*hx.nCell,1);
            hx.d = x(2*hx.nCell+1:3*hx.nCell,1);
            hx.pState = x(3*hx.nCell+1:4*hx.nCell,1);
            hx.hState = x(4*hx.nCell+1:5*hx.nCell,1);
            hx.dState = x(5*hx.nCell+1:6*hx.nCell,1);
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
            Dx = [hx.Dp; hx.Dh; hx.Dd; hx.DpState; hx.DhState; hx.DdState];
        end
        function initialize(hx,nCell,p,h,Volume,TimeConstants,ODEoptions)
            % Function help: provide two-element vectors for pressure and 
            %   enthalpy inlets and outlets, and provide volume.
                 
            hx.nCell = nCell;
            hx.totalVolume = Volume;
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