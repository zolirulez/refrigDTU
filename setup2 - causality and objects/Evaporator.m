classdef Evaporator < Tank
    properties
        % Parameters
        hOutlet
        DmInlet
    end
    methods
        function boundaryCondition(ev,DQ,hInlet)
            % Valve operation for necessary heat transfer (no dynamics)
            ev.DmInlet = -DQ/(hInlet - ev.hOutlet);
        end
        function Dx = process(ev,t,x,Inputs)
            % Time
            ev.t = t;
            % States
            ev.p = x(1,1);
            ev.h = x(2,1);
            ev.d = x(3,1);
            % Inputs
%             DQ = Inputs.DQ;
%             DmOutlet = Inputs.DmOutlet;
%             hInlet = Inputs.hInlet;
            % Process
            ev.boundaryCondition(DQ,hInlet);
            ev.massAccummulation(ev.DmInlet,DmOutlet);
            ev.excitation([hInlet; ev.hOutlet],[ev.DmInlet; -DmOutlet],[],DQ);
            ev.potentialAccummulation();
            Dx = [ev.Dp; ev.Dh; ev.Dd];
        end
        function initialize(ev,p,h,hOutlet,Volume,ODEoptions)
            % Function help: provide initial pressure, enthalpy, and the
            %   volume

            initialize@Tank(ev,p,h,Volume,ODEoptions);
            ev.hOutlet = hOutlet;
        end
    end
end