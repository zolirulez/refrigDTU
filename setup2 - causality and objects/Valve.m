classdef Valve
   properties
      Kv
   end
   methods
      function mdot = flow(valve,capacityRatio,dInlet,pInlet,pOutlet)
         % Kv value is in non-SI units!
         mdot = 8.7841e-06*capacityRatio*valve.Kv*sqrt(dInlet*(pInlet-pOutlet));
      end
      function hOutlet = enthalpy(valve,hInlet)
          hOutlet = hInlet;
      end
      function [mdot, hOutlet] = process(valve,capacityRatio,dInlet,pInlet,pOutlet,hInlet)
         mdot = flow(valve,capacityRatio,dInlet,pInlet,pOutlet);
         hOutlet = enthalpy(valve,hInlet);
      end
   end
end