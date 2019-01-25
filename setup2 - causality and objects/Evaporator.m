classdef Evaporator < Tank
    properties
        % Parameters
        hOutletVirtual
        DmInlet
        DmOutlet
        hInlet
    end
    methods
        function boundaryCondition(ev,DQ)
            % Valve operation for necessary heat transfer (no dynamics)
            ev.DmInlet = -DQ/(ev.hInlet - ev.hOutletVirtual);
        end
        function Dx = process(ev,t,x,DQ)
            % Time
            ev.t = t;
            % States
            ev.p = x(1,1);
            ev.h = x(2,1);
            ev.d = x(3,1);
            % Process
            ev.boundaryCondition(DQ);
            ev.massAccummulation();
            ev.excitation([ev.DmInlet; -ev.DmOutlet],[ev.hOutletVirtual; ev.h],0);
            ev.potentialAccummulation();
            Dx = [ev.Dp; ev.Dh; ev.Dd];
        end
        function initialize(ev,p,h,hOutletVirtual,Volume,ODEoptions)
            % Function help: provide initial pressure, enthalpy, and the
            %   volume

            initialize@Tank(ev,p,h,Volume,ODEoptions);
            ev.hOutletVirtual = hOutletVirtual;
        end
    end
end