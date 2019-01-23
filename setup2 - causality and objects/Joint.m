classdef Joint < matlab.mixin.Copyable
    properties
        p
        h
        d
        Dm
    end
    methods
        function [Dm2, h1] = noAccummulation(joint,p,Dm,h,Dm1,h2)
            % No index denotes fully defined inflow. Index 1 denotes
            %   inflow with unknown enthalpy, index 2 denotes inflow
            %   with unknown mass flow. Outflows are denoted by negative
            %   signs.
            joint.p = p;
            if nargin > 3
                Dm2 = -sum([Dm; Dm1]); 
                h1 = -[Dm; Dm2]'*[h; h2]/Dm1;
            else
                % In this case index 1 and 2 are the same
                Dm2 = -sum(Dm);
                h1 = -Dm'*h/Dm2;
            end
            joint.h = h1;
            joint.d = CoolProp.PropsSI('D','H',h1,'P',joint.p,'CO2');
            joint.Dm = Dm2;
        end
    end
end