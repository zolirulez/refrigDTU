classdef Valve < handle
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
      function [Dm, hOutlet] = process(valve,capacityRatio,Inputs)
          dInlet = Inputs.dInlet;
          pInlet = Inputs.pInlet;
          pOutlet = Inputs.pOutlet;
          hInlet = Inputs.hInlet;
          Dm = flow(valve,capacityRatio,dInlet,pInlet,pOutlet);
          hOutlet = enthalpy(valve,hInlet);
      end
      function initialize(valve,Kv)
          valve.Kv = Kv;      
      end
   end
end